#--------------------------------------#
# PDAC IHC processing pipeline 2-------#
# @Author: Haoyang Mi -----------------#
# Date: September 23th 2021 -----------#



library(ggplot2); library(dplyr); library(reshape2); library(RColorBrewer); library(matrixStats)
library(ggvoronoi); library(igraph);library(concaveman);library(sf); library(spatstat);
library(survival); library(survminer); 
library(foreach); library(doParallel); 
library(sp); library(largeVis); library(dbscan)
library(alphahull); library(caramellar)
library(ComplexHeatmap); library(circlize)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


setwd('..')

source("./Codes/Function.r")

#-----------------------------------------------#
#------------ Data Preparations ----------------#
#-----------------------------------------------#

#-- read PDAC
PDACdata <- readRDS('./exportedFiles/PDACdata_USETHIS.rds')
#-- long versus short survival
cohort <- read.csv('./exportedFiles/sample_key_10.23.csv', encoding = 'UTF-8-BE') %>%
  filter(Tumor.type == 'PDAC') %>%
  filter(Rx.type == 'naive') %>%
  filter(OS_days != 'NA') 
colnames(cohort)[1] <- 'sample'
colnames(cohort)[7] <- 'DFS_days'

# split to long versus short

thresh <- quantile(cohort$OS_days, c(0.5))


cohort$label <-  1* ifelse(cohort$OS_days > thresh, 1, 0)
#
long.group <- cohort %>%
  filter(label == 1) %>%
  select(sample) %>%
  as.matrix()


short.group <- cohort %>%
  filter(label == 0) %>%
  select(sample) %>%
  as.matrix()

# prepare file for GCN
long.group.indexed <- cbind.data.frame(long.group, 1) %>%
  `colnames<-` (c('sample', 'group'))
short.group.indexed <- cbind.data.frame(short.group, 0) %>%
  `colnames<-` (c('sample', 'group'))

survival.group <- rbind.data.frame(long.group.indexed, short.group.indexed)


km_fit <- survfit(Surv(OS_days, death) ~ label, data=cohort)

ggsurvplot(
  fit = km_fit, 
  xlab = "Days", 
  ylab = "Overall survival probability",
  risk.table = TRUE,
  pval = TRUE)


surv_pvalue(km_fit, data = cohort, method = 'survdiff')
#write.csv(survival.group, 'survival_group.csv')
#-----------------------------------------------#
#----------- Tertiary Lymphoid Structure -------#
#-----------------------------------------------#



# create some empty vectors
TLSarea_all_region_sample <- c()
for(Sample in unique(PDACdata$sample)){
  
  #Sample <- 5203
  
  # iterate through each subregion
  TLSarea_all_region <- c()
  
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  # get all ROIs
  ROIs <- unique(sampleDat$ROIs)
  
  
  for(ROI in ROIs){    
    #ROI <- 'ROI04'
    # single-cell dataset of ROI
    dat <- sampleDat %>%
      filter(ROIs == ROI)
    
    
    # construct bounding box
    Rect <- data.frame(cbind(x = c(min(dat$X), min(dat$X), max(dat$X), max(dat$X)), x = c(min(dat$Y), max(dat$Y), max(dat$Y), min(dat$Y))))
    
    p <- ggplot(dat,aes(X, Y)) +
      geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +
      scale_fill_manual(values = c('Tumor' = '#f8f8f8', 'Macrophages' = '#f8f8f8', 'Bcells' = '#82b4e5',
                                   'CD4 Tcells' = '#f8f8f8', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#f8f8f8'))  +
      coord_fixed(ratio = 1) 
    p
    
    # identify tertiary lymphoid structure
    
    memberships <- dat %>%
      filter(class == 'Bcells' |
               class == 'CD8 Tcells')
    
    if(nrow(memberships) >= 100){
      
      # combine membership cells
      #test <- hdbscan(memberships[, c('X', 'Y')], minPts = 100)
      
      # plot(memberships[ ,c('X', 'Y')], col=test$cluster + 1, pch=20)
      
      #op <- Network(memberships, prox_thresh = 100)
      
      
      Dual_EdgeList <- op[[1]]
      Dual_NodeList <- op[[2]]
      
      
      g <- graph_from_data_frame(vertices = Dual_NodeList[, c('nodes', 'x', 'y')], d= Dual_EdgeList[, c('from', 'to')], directed = FALSE) #%>%
      
      clusters <- split(names(V(g)), components(g)$membership) # all connected component
      
      plot(g,  vertex.color=rainbow(85, alpha=0.6)[components(g)$membership], vertex.label = NA, vertex.size = 2)
      
      
      
      TLSid <- Filter(function(x) {length(which(memberships[x, 'class'] == 'Bcells')) > 100}, clusters) # most connected components
      
      # for each TLS, compute TLS area, if any
      if(length(TLSid) > 0){
        TLSarea_all <- c()
        
        for(cid in names(TLSid)){
          
          # using the cluster id to get the components
          comp <- TLSid[[as.character(cid)]]
          
          # get cell data with the corresponding node data
          nodesTLS <- memberships[comp,]
          
          # use the concave hull area to approximate the TLS area
          TLSarea <- st_area(st_polygon(list(concaveman(as.matrix(nodesTLS[, c('X', 'Y')]), concavity = 2)))) / 1000000
          plot(st_polygon(list(concaveman(as.matrix(nodesTLS[, c('X', 'Y')]), concavity = 2))))
          
          #combine area
          TLSarea_all <- rbind(TLSarea_all, cbind(cid, TLSarea))
        }
        # add region data
        TLSarea_all_region <- rbind(TLSarea_all_region, cbind(TLSarea_all, ROI)) 
      }
      #--- end here -----#
    }
    
  }
  
  if(length(TLSarea_all_region) >0){
    TLSarea_all_region_sample <- rbind(TLSarea_all_region_sample, cbind(TLSarea_all_region, Sample))
  }
}










#-------------- HDBSCAN method to find clusters ----------------#
# create some empty vectors
clus.info <- c()
B.clus.info <- c()
CD8T.clus.info <- c()

# network list
network.list <- list()
nid <- 1
for(Sample in unique(PDACdata$sample)){
  
  #Sample <- 5203
  
  # iterate through each subregion
  
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  # get all ROIs
  ROIs <- unique(sampleDat$ROIs)
  
  
  for(ROI in ROIs){    
    # ROI <- 'ROI04'
    # single-cell dataset of ROI
    dat <- sampleDat %>%
      filter(ROIs == ROI)
    
    
    # construct bounding box
    Rect <- data.frame(cbind(x = c(min(dat$X), min(dat$X), max(dat$X), max(dat$X)), x = c(min(dat$Y), max(dat$Y), max(dat$Y), min(dat$Y))))
    
    #p <- ggplot(dat,aes(X, Y)) +
    #  geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +
    #  scale_fill_manual(values = c('Tumor' = '#f8f8f8', 'Macrophages' = '#f8f8f8', 'Bcells' = '#82b4e5',
    #                               'CD4 Tcells' = '#f8f8f8', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#f8f8f8'))  +
    #  coord_fixed(ratio = 1) 
    #p
    
    # identify tertiary lymphoid structure
    
    memberships <- dat %>%
      filter(class == 'Bcells' |
               class == 'CD8 Tcells')
    
    if(nrow(memberships) >= 100){
      
      # combine membership cells
      
      # HDBSCAN algorithm by DBSCAN package
      clusters <- dbscan::hdbscan(memberships[, c('X', 'Y')], minPts = 50)
      
      
      # HDBSCAN algorithm by largeVis package
      #minPts = 30
      #K = 4
      #dat <- as.matrix(memberships[, c('X', 'Y')])
      #cd8vis <- largeVis(t(dat), dim=2, K = K, tree_threshold = 100, max_iter = 5,sgd_batches = 1, threads = 1)
      #clusters <- largeVis::hdbscan(cd8vis, K=K, minPts = minPts, verbose = FALSE, threads = 2)
      
      #write.csv(memberships[, c('X', 'Y')], 'data.csv')
      cls <- cbind.data.frame(clusters$cluster, clusters$outlier_scores) %>%
        `colnames<-` (c('cluster', 'GLOSH')) %>%
        cbind.data.frame(memberships) %>%
        filter(GLOSH < 0.7) %>%
        #na.omit()
        filter(cluster != 0)
      #plot(cls[, 'X'], cls[, 'Y'], col = cls$cluster)
      
      
      # tests <- cbind.data.frame(cls.raw$cluster, cls.raw$outlier_scores) %>%
      # cbind.data.frame(memberships[ ,c('X', 'Y', 'class')]) %>%
      # `colnames<-` (c('cluster', 'GLOSH', 'X', 'Y', 'class')) 
      
      
      
      #plot(tests[, 'X'], tests[, 'Y'], col = tests$cluster)
      #colors <- mapply(function(col, i) adjustcolor(col, alpha.f = tests$GLOSH[i]), 
      #                 palette()[tests$cluster+1], seq_along(tests$cluster))
      #points(tests[tests$cluster == 2, 'X'], tests[tests$cluster == 2, 'Y'], col=colors, pch=20)
      
      #cls <- cls.raw$cluster %>%
      #  cbind.data.frame(cls.raw$outlier_scores, memberships[ ,c('X', 'Y', 'class')]) %>%
      #  `colnames<-` (c('cluster', 'GLOSH', 'X', 'Y', 'class')) %>%
      #  filter(cluster != 0) %>%
      #  filter(GLOSH < 0.7)
      
      
      
      #cls.id <- length(which(!is.na(unique(clusters$cluster))))
      cls.id <- unique(cls$cluster)
      
      #plot(cls$X, cls$Y, col = cls$cluster)
      
      if(!is.empty(cls.id)){
        for(cl in cls.id){
          
          #cl <- 1
          # number of members within this cluster
          #n.member <- clusdata %>%
          #  filter(cluster == cl) %>%
          #  nrow()
          
          # cl <- 2
          # cluster data
          cls.data <- cls %>%
            filter(cluster == cl)
          
          if(nrow(cls.data[cls.data$class == 'Bcells', ]) > 100){
            
            # voronoi tessellation
            cls.data$id <- seq_len(nrow(cls.data))
            adj.voronoi <- voronoi_adjacency(cls.data, id~ X + Y)$Adjacencies %>%
              apply(., 1, function(data) which(data == "TRUE"))
            
            voronoi.network <- data.frame(matrix(nrow = 0, ncol = 0))
            for(id in seq_len(length(adj.voronoi))){
              #
              voronoi.network <- rbind.data.frame(voronoi.network, cbind.data.frame(id, adj.voronoi[[id]]))
              
            }
            
            
            voronoi.network <- voronoi.network %>%
              `colnames<-` (c('source', 'target')) %>%
              mutate(from = pmin(source, target),
                     to = pmax(source, target)) %>%
              distinct(from, to)
            
            
            # save node and link data to list
            
            network.list[[nid]] <- list('node' = cls.data, 'link' = voronoi.network, 'Sample' = Sample, 'ROI' = ROI, 'cid' = cl)
            
            nid <- nid + 1
            
            # plot data frame graph
            
            g <- graph_from_data_frame(vertices = cls.data[, c('id', 'X', 'Y')], d= voronoi.network[, c('from', 'to')], directed = FALSE) #%>%
            
            
            p <- ggplot() +
              geom_segment(aes(x = cls.data[voronoi.network$from, 'X'], y = cls.data[voronoi.network$from, 'Y'],
                               xend = cls.data[voronoi.network$to, 'X'], yend = cls.data[voronoi.network$to, 'Y']),
                           color = 'grey') +
              geom_point(data = cls.data, shape = 21, aes(X, Y, fill = class), size = 5) +
              scale_fill_manual(values = c('Bcells' = '#82b4e5', 'CD8 Tcells' = '#4875b0')) +
              theme_void() +
              coord_fixed(ratio = 1) +
              theme(legend.position = 'none')
            #plot(p)
            #ggsave(p, file=paste0('Figures/TLS/TLS_networks/', Sample, "_", ROI, "_", cl, ".png"), width = 10, height = 10, units = "in", dpi = 300)
            
            # proximity method
            # op <- Network(cls.data, prox_thresh = 300)
            # 
            # Dual_EdgeList <- op[[1]]
            # Dual_NodeList <- op[[2]]
            # Dual_NodeList$class <- cls.data[Dual_NodeList$nodes, 'class']
            # 
            # g <- graph_from_data_frame(vertices = Dual_NodeList[, c('nodes', 'x', 'y')], d= Dual_EdgeList[, c('from', 'to')], directed = FALSE) #%>%
            # 
            # 
            # p <- ggplot() +
            #   geom_segment(aes(x = Dual_NodeList[Dual_EdgeList$from, 'x'], y = Dual_NodeList[Dual_EdgeList$from, 'y'],
            #                    xend = Dual_NodeList[Dual_EdgeList$to, 'x'], yend = Dual_NodeList[Dual_EdgeList$to, 'y']),
            #                color = 'grey') +
            #   geom_point(data = Dual_NodeList, aes(x, y, color = class), size = 4) +
            #   scale_color_manual(values = c('Bcells' = '#82b4e5', 'CD8 Tcells' = '#4875b0')) +
            #   theme_void() +
            #   coord_fixed(ratio = 1) +
            #   theme(legend.position = 'none')
            # plot(p)
            # 
            # ggsave(p, file=paste0('Figures/TLS/TLS_networks/', Sample, "_", ROI, "_", cl, ".png"), width = 10, height = 10, units = "in", dpi = 300)
            # 
            # 
            # 
            
            
            #plot(g,  vertex.color=rainbow(85, alpha=0.6)[components(g)$membership], vertex.label = NA, vertex.size = 2)
            #plot(st_polygon(list(concaveman(as.matrix(cls.data[, c('X', 'Y')]), concavity = 1.5))))
            
            # calculate area
            clus.morph <- getEllipse(cls.data)
            
            
            #TLSarea <- st_area(st_polygon(list(concaveman(as.matrix(cls.data[, c('X', 'Y')]), concavity = 4)))) / 1000000
            
            # calculate cell density per cluster
            B_density <- nrow(cls.data[cls.data$class == 'Bcells',]) / clus.morph$Aconv
            CD8T_density <- nrow(cls.data[cls.data$class == 'CD8 Tcells',]) / clus.morph$Aconv
            
            
            # compute ratio of each cell type
            IRPvars = c("PD1.", "PDL1.", 'KI67.', 'EOMES.', 'GZMB.', "IL10.", 'ICOS.')
            
            exprData <- cls.data %>%
              summarise_each_(funs(sum(.==1,na.rm=TRUE)), IRPvars)
            exprData <- exprData/sum(exprData)
            
            clus.info <- rbind.data.frame(clus.info, cbind.data.frame(cl, B_density, CD8T_density, clus.morph, exprData, ROI, Sample))
            
            # B cell specific
            exprData <- cls.data %>%
              filter(class == 'Bcells') %>%
              summarise_each_(funs(sum(.==1,na.rm=TRUE)), IRPvars)
            exprData <- exprData/sum(exprData)
            
            B.clus.info <- rbind.data.frame(B.clus.info, cbind.data.frame(cl, B_density, CD8T_density, clus.morph, exprData, ROI, Sample))
            
            # CD8 T cell specific
            exprData <- cls.data %>%
              filter(class == 'CD8 Tcells') %>%
              summarise_each_(funs(sum(.==1,na.rm=TRUE)), IRPvars)
            exprData <- exprData/sum(exprData)
            CD8T.clus.info <- rbind.data.frame(CD8T.clus.info, cbind.data.frame(cl, B_density, CD8T_density, clus.morph, exprData, ROI, Sample))
            
            
          }
          
          
        }
        
      }
      
      #--- end here -----#
    }
  }
  
  #if(length(TLSarea_all_region) >0){
  #  TLSarea_all_region_sample <- rbind(TLSarea_all_region_sample, cbind(TLSarea_all_region, Sample))
  #}
}


#Sample <- 5203

#ROI <- 'ROI04'


#ROI <- 'ROI05'
# single-cell dataset of ROI
dat <- PDACdata %>%
  filter(sample == Sample) %>%
  filter(ROIs == ROI)


# construct bounding boxr
Rect <- data.frame(cbind(x = c(min(dat$X), min(dat$X), max(dat$X), max(dat$X)), x = c(min(dat$Y), max(dat$Y), max(dat$Y), min(dat$Y))))

p <- ggplot(dat,aes(X, Y)) +
  geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +
  scale_fill_manual(values = c('Tumor' = '#f8f8f8', 'Macrophages' = '#f8f8f8', 'Bcells' = '#82b4e5',
                               'CD4 Tcells' = '#f8f8f8', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#f8f8f8'))  +
  #theme_void() +
  theme(legend.position = 'none') +
  coord_fixed(ratio = 1) 
p
ggsave(p, file=paste0('Figures/TLS_5203_ROI04.png'), width = 15, height = 10, units = "in", dpi = 300)





#---------------------------------#
# Voronoi graph for ROIs with TLS #
#---------------------------------#



for(sp in unique(TLSarea_all$Sample)){
  
  sp <- 2148
  
  # iterate through each subregion
  
  # get all ROIs
  ROIs <- unique(TLSarea_all[TLSarea_all$Sample == sp, 'ROI'])
  
  
  for(ROI in ROIs){    
    ROI <- 'ROI08'
    # single-cell dataset of ROI
    
    label <- unique(TLSarea_all[TLSarea_all$Sample == sp, 'group'])
    
    dat <- PDACdata %>%
      filter(sample == sp) %>%
      filter(ROIs == ROI)
    
    
    # construct bounding box
    Rect <- data.frame(cbind(x = c(min(dat$X), min(dat$X), max(dat$X), max(dat$X)), x = c(min(dat$Y), max(dat$Y), max(dat$Y), min(dat$Y))))
    
    
    # only B cells and CD8 T cells
    p <- ggplot(dat,aes(X, Y)) +
      theme_void() +
      geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +
      scale_fill_manual(values = c('Tumor' = '#f8f8f8', 'Macrophages' = '#f8f8f8', 'Bcells' = '#82b4e5',
                                   'CD4 Tcells' = '#f8f8f8', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#f8f8f8'))  +
      theme(legend.position = 'none') +
      coord_fixed(ratio = 1) 
    p
    ggsave(p, file=paste0('Figures/TLS/noTumor/', label, '_', sp, '_', ROI, '.png'), width = 10, height = 10, units = "in", dpi = 300)
    
    # B cells, CD8 T cells, and Tumor cells
    p <- ggplot(dat,aes(X, Y)) +
      theme_void() +
      geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +
      scale_fill_manual(values = c('Tumor' = '#fb8648', 'Macrophages' = '#f8f8f8', 'Bcells' = '#82b4e5',
                                   'CD4 Tcells' = '#f8f8f8', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#f8f8f8'))  +
      coord_fixed(ratio = 1) +
      theme(legend.position = 'none')
    p
    ggsave(p, file=paste0('Figures/TLS/Tumor/', label, '_', sp, '_', ROI, '.png'), width = 10, height = 10, units = "in", dpi = 300)
    
    # identify tertiary lymphoid structure
  }
}






#------------------------------------------#
#------- Figure plot ---------#
#------------------------------------------#


clus.info$group <- 'long'
clus.info[clus.info$Sample %in% short.group, 'group'] <- 'short'
clus.info$group <- factor(clus.info$group, levels = c('short', 'long'))
#write.csv(clus.info, 'clus_info_long_short.csv')
wilcox.test(clus.info[clus.info$group == 'long', 'Aconv'], clus.info[clus.info$group == 'short', 'Aconv'])

# I: TLS/B cell follicle size
p <- ggplot(clus.info, aes(x = group, y = Aconv, color = group)) +
  theme_classic() +
  stat_boxplot( aes(group, Aconv), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = group), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 2) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 15, annotations = 'ns') +
  ylim(0, 16) +
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  scale_fill_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none') 

p
ggsave(p, file=paste0('Figures/TLS_area_compare.png'), width = 3, height = 4, units = "in", dpi = 300)


# II: B cell density per cluster
wilcox.test(clus.info[clus.info$group == 'long', 'B_density'], clus.info[clus.info$group == 'short', 'B_density'])

p <- ggplot(clus.info, aes(x = group, y = B_density, color = group)) +
  theme_classic() +
  stat_boxplot( aes(group, B_density), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = group), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 2) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 2400, annotations = 'ns') +
  ylim(0, 2700) +
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  scale_fill_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none') 

p
ggsave(p, file=paste0('Figures/TLS_Bdensity_compare.png'), width = 3, height = 4, units = "in", dpi = 300)




# III: T cell density per cluster
wilcox.test(clus.info[clus.info$group == 'long', 'CD8T_density'], clus.info[clus.info$group == 'short', 'CD8T_density'])


p <- ggplot(clus.info, aes(x = group, y = CD8T_density, color = group)) +
  theme_classic() +
  stat_boxplot( aes(group, CD8T_density), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = group), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 2) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 1600, annotations = '**') +
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  scale_fill_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none') +
  ylim(0, 1800)

p
ggsave(p, file=paste0('Figures/TLS_CD8Tdensity_compare.png'), width = 3, height = 4, units = "in", dpi = 300)

# IV: Eccentricity
wilcox.test(clus.info[clus.info$group == 'long', 'eccentricity'], clus.info[clus.info$group == 'short', 'eccentricity'])

p <- ggplot(clus.info, aes(x = group, y = eccentricity, color = group)) +
  theme_classic() +
  stat_boxplot( aes(group, eccentricity), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = group), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 2) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 1, annotations = 'ns') +
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  scale_fill_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none') +
  ylim(0, 1.2)

p
ggsave(p, file=paste0('Figures/TLS_CD8Tdensity_compare.png'), width = 3, height = 4, units = "in", dpi = 300)



# V: Circularity
wilcox.test(clus.info[clus.info$group == 'long', 'circularity'], clus.info[clus.info$group == 'short', 'circularity'])

p <- ggplot(clus.info, aes(x = group, y = circularity, color = group)) +
  theme_classic() +
  stat_boxplot( aes(group, circularity), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = group), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 2) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 0.9, annotations = 'ns') +
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  scale_fill_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none') +
  ylim(0, 1)

p
ggsave(p, file=paste0('Figures/TLS_CD8Tdensity_compare.png'), width = 3, height = 4, units = "in", dpi = 300)


# VI: Convexity
wilcox.test(clus.info[clus.info$group == 'long', 'convexity'], clus.info[clus.info$group == 'short', 'convexity'])

p <- ggplot(clus.info, aes(x = group, y = convexity, color = group)) +
  theme_classic() +
  stat_boxplot( aes(group, convexity), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = group), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 2) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 1, annotations = 'ns') +
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  scale_fill_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none') +
  ylim(0.4, 1.1)

p
ggsave(p, file=paste0('Figures/TLS_CD8Tdensity_compare.png'), width = 3, height = 4, units = "in", dpi = 300)



#--------------------------------------#
#----------- Network analysis ---------#
#--------------------------------------#







#----------------------------------------------#
#---- Component marker expression analysis ----#
#----------------------------------------------#

B_CD8T.clus.info <- cbind.data.frame(CD8T.clus.info, B.clus.info)
B_CD8T.clus.info[B_CD8T.clus.info == "NaN"] <- 0
B_CD8T.clus.info$group <- 1 * (B_CD8T.clus.info$Sample %in% long.group)

B_CD8T.clus.info[B_CD8T.clus.info$group == '1', 'group'] <- 'long'
B_CD8T.clus.info[B_CD8T.clus.info$group == '0', 'group'] <- 'short'


# rename to add appendix
colnames(B_CD8T.clus.info)[9:15] <- paste0(colnames(B_CD8T.clus.info)[9:15], '_CD8T')
colnames(B_CD8T.clus.info)[26:32] <- paste0(colnames(B_CD8T.clus.info)[26:32], '_B')


png('Figures/Clustering_ComplexHeatmap.jpeg', width=22, height=35, units = 'cm', res = 300) 


foo2 = as.matrix(B_CD8T.clus.info[, 19:20])
row_ha = rowAnnotation('Density, mm' = row_anno_barplot(foo2, axis = TRUE, axis_param = list(side = "top"),
                                                        gp = gpar(fill = c("#82b4e5", "#4875b0"))), width = unit(4, "cm"))

lgd = Legend(at = c("CD8 T", "B cell"), title = "Cell type", legend_gp = gpar(fill = c("#4875b0", "#82b4e5")))


q <- Heatmap(as.matrix(B_CD8T.clus.info[, c(9:15, 26:32)]), name="ratio",
             col = colorRamp2(c(0, 1), c("white", "#377EB8")),
             width = unit(8, "cm"), 
             height = unit(20, "cm"), 
             row_km = 2, column_km = 3,
             border = TRUE, 
             row_dend_width = unit(4, 'cm'),
             column_dend_height = unit(4, 'cm'),
             row_title = NULL,
             column_title = NULL,
             heatmap_legend_param = list(legend_height = unit(4, "cm")),
             left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c('#4986c3', 'black')),
                                                              labels = c("Cluster 1", "Cluster 2"), 
                                                              labels_gp = gpar(col = "white", fontsize = 10))),             #cell_fun = TextFunc(pvalMatrices.vec, numdat = F)
) +
  Heatmap(B_CD8T.clus.info$label, name = "Cluster shape", width = unit(5, "mm"), col = c("circular" = "#a2b1d2", "elongated" = "#76b666", 'else' = '#f0a47c')) +
  Heatmap(B_CD8T.clus.info$group, name = "Survival gorup", width = unit(5, "mm"), col = c("short" = "#c22828", "long" = "#325698")) +
  row_ha 
q

draw(q, heatmap_legend_list = list(lgd))
dev.off()


#----------------------------------------------------------#
# assign cluster (cluster 1 or 2) to each T/B cell cluster #
#----------------------------------------------------------#


clusterlist = row_order(q)


chiseq.clus.info <- B_CD8T.clus.info

chiseq.clus.info$cluster <- 'c1'
chiseq.clus.info[clusterlist[[2]], 'cluster'] <- 'c2'

test <- chiseq.clus.info[, c('group', 'cluster')] %>%
  table()
test
chisq.test(test)



test1 <- chiseq.clus.info[, c('label', 'cluster')] %>%
  table()
test1
chisq.test(test1)
# 
c1 <- chiseq.clus.info[chiseq.clus.info$cluster == 'c1',]
c2 <- chiseq.clus.info[chiseq.clus.info$cluster == 'c2',]
wilcox.test(c1$B_density, c2$B_density)

mean(c1$B_density)
mean(c2$B_density)
chiseq.clus.info_2 <- chiseq.clus.info[, c(1, 2, 3, 16, 17, 36)]
p <- ggplot(chiseq.clus.info_2, aes(x = cluster, y = B_density, color = cluster)) +
  theme_classic() +
  stat_boxplot( aes(cluster, B_density), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = cluster), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 3) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("c1", "c2")), 
              map_signif_level=TRUE, size = 0.5, textsize = 11, color = 'black', 
              y_position = 2300, annotations = '**') +
  #ylim(0, 16) +
  scale_color_manual(values = c('c1' = '#4986c3', 'c2' = 'black')) +
  scale_fill_manual(values = c('c1' = '#4986c3', 'c2' = 'black')) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 30),
        legend.position = 'none') 

p
ggsave(p, file=paste0('Figures/B_density_cluster1.png'), width = 6, height = 12, units = "in", dpi = 300)

wilcox.test(c1$CD8T_density, c2$CD8T_density)


p <- ggplot(chiseq.clus.info_2, aes(x = cluster, y = CD8T_density, color = cluster)) +
  theme_classic() +
  stat_boxplot( aes(cluster, B_density), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = cluster), alpha = 0.2, outlier.shape =  NA) +
  geom_jitter(size = 3) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("c1", "c2")), 
              map_signif_level=TRUE, size = 0.5, textsize = 11, color = 'black', 
              y_position = 2400, annotations = '****') +
  #ylim(0, 16) +
  scale_color_manual(values = c('c1' = '#4986c3', 'c2' = 'black')) +
  scale_fill_manual(values = c('c1' = '#4986c3', 'c2' = 'black')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 30),
        legend.position = 'none') 

p
ggsave(p, file=paste0('Figures/CD8T_density_cluster1.png'), width = 6, height = 12, units = "in", dpi = 300)

#----------------------------------------#
#--------- Network properties -----------#
#----------------------------------------#

bet_CD8_all <- data.frame(matrix(nrow = 0, ncol = 0))
bet_B_all <- data.frame(matrix(nrow = 0, ncol = 0))

clos_CD8_all <- data.frame(matrix(nrow = 0, ncol = 0))
clos_B_all <- data.frame(matrix(nrow = 0, ncol = 0))

deg_CD8_all <- data.frame(matrix(nrow = 0, ncol = 0))
deg_B_all <- data.frame(matrix(nrow = 0, ncol = 0))

trans_CD8_all <- data.frame(matrix(nrow = 0, ncol = 0))
trans_B_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(nid in seq_len(length(network.list))){
  
  # cluster data
  cls.data <- network.list[[nid]]$node
  
  # cluster ID
  cl <- network.list[[nid]]$cid
  
  # voronoi links
  voronoi.network <- network.list[[nid]]$link
  
  # Sample
  Sample <- network.list[[nid]]$Sample
  
  # ROI
  ROI <- network.list[[nid]]$ROI
  
  # B cell ids and T cell ids
  Bid <- which(cls.data$class == 'Bcells')
  Tid <- which(cls.data$class == 'CD8 Tcells')
  
  
  #----- Construct igraph object ------#
  
  g <- graph_from_data_frame(vertices = cls.data[, c('id', 'X', 'Y')], d = voronoi.network[, c('from', 'to')], directed = FALSE) #%>%
  
  mst.g <- mst(g)
  plot(mst.g, vertex.label = NA, vertex.size = 2)
  #---- calculate properties ------#
  
  
  # expression levels
  expr_CD8 <- cls.data[Tid, 32:38]
  expr_B <- cls.data[Bid, 32:38]
  
  # betweeness
  bet_CD8 <-  betweenness(mst.g, v = V(mst.g), directed = FALSE)[Tid]
  bet_B <-  betweenness(mst.g, v = V(mst.g), directed = FALSE)[Bid]
  
  # closeness
  clos_CD8 <- closeness(mst.g, v = V(mst.g))[Tid]
  clos_B <- closeness(mst.g, v = V(mst.g))[Bid]
  
  # degree
  deg_CD8 <- igraph::degree(mst.g, v = V(mst.g))[Tid]
  deg_B <- igraph::degree(mst.g, v = V(mst.g))[Bid]
  
  # transitivity (clustering coefficient)
  trans_CD8 <- transitivity(mst.g, type = 'local')[Tid]
  trans_B <- transitivity(mst.g, type = 'local')[Bid]
  
  #---- combine altogether ------#
  bet_CD8_all <- rbind.data.frame(bet_CD8_all, cbind.data.frame(bet_CD8, expr_CD8, ROI, Sample, cl))
  bet_B_all <- rbind.data.frame(bet_B_all, cbind.data.frame(bet_B, expr_B, ROI, Sample, cl))
  
  clos_CD8_all <- rbind.data.frame(clos_CD8_all, cbind.data.frame(clos_CD8, expr_CD8, ROI, Sample, cl))
  clos_B_all <- rbind.data.frame(clos_B_all, cbind.data.frame(clos_B, ROI, expr_B, Sample, cl))
  
  deg_CD8_all <- rbind.data.frame(deg_CD8_all, cbind.data.frame(deg_CD8, expr_CD8, ROI, Sample, cl))
  deg_B_all <- rbind.data.frame(deg_B_all, cbind.data.frame(deg_B, ROI, expr_B, Sample, cl))
  
  trans_CD8_all <- rbind.data.frame(trans_CD8_all, cbind.data.frame(trans_CD8, expr_CD8, ROI, Sample, cl))
  trans_B_all <- rbind.data.frame(trans_B_all, cbind.data.frame(trans_B, expr_B, ROI, Sample, cl))
  
}



# merge centrality metrics to cluster id

bet_CD8_all_clusid <- merge(bet_CD8_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))
bet_B_all_clusid <- merge(bet_B_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))

trans_CD8_all_clusid <- merge(trans_CD8_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))
trans_B_all_clusid <- merge(trans_B_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))

clos_CD8_all_clusid <- merge(clos_CD8_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))
clos_B_all_clusid <- merge(clos_B_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))

deg_CD8_all_clusid <- merge(deg_CD8_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))
deg_B_all_clusid <- merge(deg_B_all, chiseq.clus.info_2, by = c('cl', 'Sample', 'ROI'))

bet_CD8_all_clus1 <- bet_CD8_all_clusid %>%
  filter(cluster == 'c1')

bet_CD8_all_clus2 <- bet_CD8_all_clusid %>%
  filter(cluster == 'c2')


bet_B_all_clus1 <- bet_B_all_clusid %>%
  filter(cluster == 'c1')

bet_B_all_clus2 <- bet_B_all_clusid %>%
  filter(cluster == 'c2')



clos_CD8_all_clus1 <- clos_CD8_all_clusid %>%
  filter(cluster == 'c1')

clos_CD8_all_clus2 <- clos_CD8_all_clusid %>%
  filter(cluster == 'c2')


clos_B_all_clus1 <- clos_B_all_clusid %>%
  filter(cluster == 'c1')

clos_B_all_clus2 <- clos_B_all_clusid %>%
  filter(cluster == 'c2')


deg_CD8_all_clus1 <- deg_CD8_all_clusid %>%
  filter(cluster == 'c1')

deg_CD8_all_clus2 <- deg_CD8_all_clusid %>%
  filter(cluster == 'c2')


deg_B_all_clus1 <- deg_B_all_clusid %>%
  filter(cluster == 'c1')

deg_B_all_clus2 <- deg_B_all_clusid %>%
  filter(cluster == 'c2')


# clustering coefficient
trans_CD8_all_clus1 <- trans_CD8_all_clusid %>%
  filter(cluster == 'c1')

trans_CD8_all_clus2 <- trans_CD8_all_clusid %>%
  filter(cluster == 'c2')


trans_B_all_clus1 <- trans_B_all_clusid %>%
  filter(cluster == 'c1')

trans_B_all_clus2 <- trans_B_all_clusid %>%
  filter(cluster == 'c2')


p <- ggplot() +
  geom_density(data = clos_CD8_all_clus1[clos_CD8_all_clus1$IL10. == 1, ], aes(clos_CD8, ..scaled..), color = '#71c8d7', size = 2) + # blue
  geom_density(data = clos_CD8_all_clus1[clos_CD8_all_clus1$IL10. == 0, ], aes(clos_CD8, ..scaled..), color = '#d16c6c', size = 2) +  # red
  theme_classic() +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 26)) +
  ylab('Density') +
  xlab('closeness')
p
ggsave(p, file=paste0('Figures/closeness_CD8_cluster1.png'), width = 10, height = 6, units = "in", dpi = 300)

p <- ggplot() +
  geom_density(data = clos_B_all_clus2[clos_B_all_clus2$EOMES == 1, ], aes(clos_B, ..scaled..), color = '#71c8d7', size = 2) +
  geom_density(data = clos_B_all_clus2[clos_B_all_clus2$EOMES. == 0, ], aes(clos_B, ..scaled..), color = '#d16c6c', size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 26)) +
  ylab('Density') +
  xlab('closeness')
p
ggsave(p, file=paste0('Figures/closeness_Bcell_cluster2.png'), width = 10, height = 6, units = "in", dpi = 300)




#-------- Boxplot for closeness comparison -----------#
p <- ggplot(clos_B_all_clus2, aes(x = as.factor(EOMES.), y = clos_B, color = as.factor(EOMES.))) +
  theme_classic() +
  stat_boxplot( aes(as.factor(EOMES.), clos_B), 
                geom='errorbar', linetype=1, width=0.3, size =1)+  #whiskers
  geom_boxplot(size = 1, aes(fill = as.factor(EOMES.)), alpha = 0.2, outlier.shape =  NA) +
  #geom_jitter(size = 3) +
  #geom_violin(size = 1) +
  geom_signif(comparisons = list(c("0", "1")), 
              map_signif_level=TRUE, size = 0.5, textsize = 11, color = 'black', 
              y_position = 0.001, annotations = '****') +
  #ylim(0, 16) +
  scale_color_manual(values = c('c1' = '#4986c3', 'c2' = 'black')) +
  scale_fill_manual(values = c('c1' = '#4986c3', 'c2' = 'black')) +
  
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 30),
        legend.position = 'none') 

p
ggsave(p, file=paste0('Figures/CD8T_density_cluster1.png'), width = 6, height = 12, units = "in", dpi = 300)

##
mean(chiseq.clus.info[chiseq.clus.info$cluster == 'c2', 2] / chiseq.clus.info[chiseq.clus.info$cluster == 'c2', 3])
