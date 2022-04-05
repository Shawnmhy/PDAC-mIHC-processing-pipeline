#--------------------------------------#
# PDAC IHC processing pipeline validation-------#
# @Author: Haoyang Mi -----------------#
# Date: December 1st 2021 -----------#



library(ggplot2); library(dplyr); library(reshape2); library(RColorBrewer); library(matrixStats)
library(ggvoronoi); library(igraph);library(concaveman);library(sf); library(spatstat);
library(survival); library(survminer); 
library(foreach); library(doParallel); 
library(sp); library(dbscan)
library(alphahull); library(caramellar)
library(ComplexHeatmap); library(circlize)
library(Rmisc)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd('..')
source("./Codes/Function.r")

# how many ROIs per sample

cohort <- read.csv('./exportedFiles/sample_key_10.23.csv', encoding = 'UTF-8-BE') %>%
  filter(Tumor.type == 'PDAC') %>%
  filter(Rx.type == 'naive') #%>%
  #filter(OS_days != 'NA')
colnames(cohort)[1] <- 'sample'
colnames(cohort)[7] <- 'DFS_days'


# long versus short survival group
thresh <- 619
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






sampleID <- cohort$sample

# read single-cell dataset
allFiles <- list.files('./Dataset_2021/ClassifiedCSV_10012021T8') 



# create some empty vectors
dat_regions_samples <- c()
for(sample in sampleID){
  
  #sample <- 5218
  subsets <- grep(sample, allFiles, value = TRUE)
  
  # number of regions
  nregion <- length(subsets)
  
  # iterate through each subregion
  dat_regions <- c()
  
  for(subregion in subsets){
    
    #subregion <- subsets[8]
    # single-cell dataset of subregion
    dat <- read.csv(paste0('./Dataset_2021/ClassifiedCSV_10012021T8/', subregion), row.names = 1) %>%
      filter(class != 'Noise')
    
    
    # construct bounding box
    
    ROIs <- rev(strsplit(subregion, "\\.|\\,|_")[[1]])[2]
    
    dat_regions <- rbind.data.frame(dat_regions, cbind.data.frame(dat, ROIs))
  }
  
  dat_regions_samples <- rbind.data.frame(dat_regions_samples, cbind.data.frame(dat_regions, sample))
}

colnames(dat_regions_samples)[30:36] <- c('PD1.', 'PDL1.', 'KI67.', 'EOMES.', 'GZMB.', 'IL10.', 'ICOS.')
colnames(dat_regions_samples)[27:28] <- c('X', 'Y')

#saveRDS(dat_regions_samples, 'PDACdata_USETHIS_validation.rds')





# voronoi tessellation
# create some empty vectors
for(sample in sampleID){
  
  #sample <- 2009
  subsets <- grep(sample, allFiles, value = TRUE)
  
  # number of regions
  nregion <- length(subsets)
  
  # iterate through each subregion
  TLSarea_all_region <- c()
  
  for(subregion in subsets){
    
    #ROIs
    ROIs <- rev(strsplit(subregion, "\\.|\\,|_")[[1]])[2]
    
    #subregion <- subsets[6]
    # single-cell dataset of subregion
    dat <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class != 'Noise')
    
    
    # construct bounding box
    Rect <- data.frame(cbind(x = c(min(dat$X), min(dat$X), max(dat$X), max(dat$X)), x = c(min(dat$Y), max(dat$Y), max(dat$Y), min(dat$Y))))
    
    p <- ggplot(dat,aes(X, Y)) +
      theme_void() +
      theme(legend.position = 'NA') +
      geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +
      scale_fill_manual(values = c('Tumor' = '#fb8648', 'Macrophages' = '#c074eb', 'Bcells' = '#82b4e5',
                                   'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#ebe459'))  +
      
      coord_fixed(ratio = 1) 
    p
    ggsave(p, file=paste0("Figures/Voronoi_Tessellation_Validation/", sample,'_', ROIs, '.png'), width = 8, height = 6, units = "in", dpi = 300)
    
    
    
  }
}




#---------------------------------------#
#---- CD8+ T cell effective score ------#
#---------------------------------------#
PDACdata <- readRDS('./exportedFiles/PDACdata_USETHIS_validation.rds')

eff.score_all <- data.frame(matrix(nrow = 0, ncol = 0))

CD8_all <- PDACdata[PDACdata$class == 'CD8 Tcells' ,] %>%
  filter(GZMB. == 1)
quantile(CD8_all$GZMB, probs = 0.75)
#quantile(CD8_all$PD1)

for(Sample in unique(PDACdata$sample)){
  
  #Sample <- 2089
  # sample data
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  # get all ROIs
  ROIs <- unique(sampleDat$ROIs)
  
  # CD4 -> Macrophage
  
  
  for(ROI in ROIs){
    #ROI <- 'ROI01'
    # roi data
    roiDat <- sampleDat %>%
      filter(ROIs == ROI)
    
    
    # get single-cell data for PD1+ CD4 T, IL10+ Macrophage, and Tumor cells
    CD8T <- roiDat %>%
      filter(class == 'CD8 Tcells') %>%
      filter(GZMB. == 1)
    
    
    IL10_Mpg <- roiDat %>%
      filter(class == 'Macrophages') %>%
      filter(IL10. == 1) %>%
      filter(GZMB. == 0) %>%
      filter(PD1. == 0) %>%
      filter(PDL1. == 0) %>%
      filter(ICOS. == 0) %>%
      filter(KI67. == 0) %>%
      filter(EOMES. == 0)
    
    CD4T <- roiDat %>%
      filter(class == 'CD4 Tcells') %>%
      filter(PD1. == 1) %>%
      filter(GZMB. == 0) %>%
      filter(IL10. == 0) %>%
      filter(PDL1. == 0) %>%
      filter(ICOS. == 0) %>%
      filter(KI67. == 0) %>%
      filter(EOMES. == 0)
    #filter(IL10. == 1)
    
    # distance to nearest IL10+ macrophage
    if(nrow(CD8T) >0 & nrow(CD4T) > 0 & nrow(IL10_Mpg) > 0){
      dist2Mac <- flexclust::dist2(IL10_Mpg[, c('X', 'Y')], CD8T[, c('X', 'Y')]) %>%
        as.matrix() %>%
        rowMins()
      
      
      # distance to nearest Tumor cell
      
      dist2Tumor <- flexclust::dist2(IL10_Mpg[, c('X', 'Y')], CD4T[, c('X', 'Y')]) %>%
        as.matrix() %>%
        rowMins()
      
      # Effective score
      eff.score <- dist2Tumor / (dist2Tumor + dist2Mac)
      
      # combine all scores
      
      eff.score_all <- rbind.data.frame(eff.score_all, cbind.data.frame(eff.score, ROI, Sample))
      
    }
    
  }
  print(Sample)
}


eff.scores <- eff.score_all %>%
  #group_by(Sample) %>%
  #dplyr::summarise(mean = mean(eff.score)) %>%
  mutate(class = 1 * (Sample %in% long.group))


eff.scores[eff.scores$class == '1', 'class'] <- 'long'
eff.scores[eff.scores$class == '0', 'class'] <- 'short'

tgc <- summarySE(eff.scores, measurevar= "eff.score", groupvars=c("class"))
tgc$class <- factor(tgc$class, level = c('short', 'long'))

pd <- position_dodge(0.1) # move them .05 to the left and right
wilcox.test(eff.scores[eff.scores$class == 'short', 'eff.score'], eff.scores[eff.scores$class == 'long', 'eff.score'])

p <- ggplot(tgc, aes(x= class, y = eff.score, colour=class, group=class)) + 
  theme_bw() +
  geom_errorbar(aes(ymin=eff.score - 10*se, ymax=eff.score + 10*se, color = class), width=.05) +
  geom_point(aes(color = class), fill = 'white', position=pd, size=3, shape=21) + # 21 is filled circle
  
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  ylab(expression(IL-10^'+' ~ macrophage~risk~score)) +
  expand_limits(y=0) +                        # Expand y range
  ylim(0, 1) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 0.9, annotations = '****'
  ) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.position = 'none') 
p
ggsave(p, file=paste0("Figures/Macrophage_Risk_Score_GZMB.jpg"), width = 5, height = 6, units = "in", dpi = 300)



mean(eff.scores[eff.scores$class == 'short', 'eff.score'])
mean(eff.scores[eff.scores$class == 'long', 'eff.score'])






#-------------------------------------------#
#---- Graph of links to CD8 and CD4 T ------#
#-------------------------------------------#

for(Sample in unique(PDACdata$sample)){
  
  #Sample <- 5213
  # sample data
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  # get all ROIs
  ROIs <- unique(sampleDat$ROIs)
  
  
  for(ROI in ROIs){
    #ROI <- 'ROI06'
    # roi data
    roiDat <- sampleDat %>%
      filter(ROIs == ROI)
    
    
    
    CD8_linkData <- roiDat %>%
      filter(class == 'CD8 Tcells') %>%
      filter(GZMB. == 1)
    
    
    
    M_linkData <- roiDat %>%
      filter(class == 'Macrophages') %>%
      filter(IL10. == 1) %>%
      filter(GZMB. == 0) %>%
      filter(PD1. == 0) %>%
      filter(PDL1. == 0) %>%
      filter(ICOS. == 0) %>%
      filter(KI67. == 0) %>%
      filter(EOMES. == 0)
    
    CD4_linkData <- roiDat %>%
      filter(class == 'CD4 Tcells') %>%
      filter(PD1. == 1) %>%
      filter(GZMB. == 0) %>%
      filter(IL10. == 0) %>%
      filter(PDL1. == 0) %>%
      filter(ICOS. == 0) %>%
      filter(KI67. == 0) %>%
      filter(EOMES. == 0)
    
    
    #testFD <- rbind.data.frame(CD8_linkData, M_linkData, CD4_linkData)
    
    #test <- testFD %>%
    #  group_by(ROIs, sample, class) %>%
    #  tally()
    
    if(nrow(CD8_linkData) >0 & nrow(CD4_linkData) > 0 & nrow(M_linkData) > 0){
      # compute links
      dist1 <- flexclust::dist2(M_linkData[, c('X', 'Y')], CD8_linkData[, c('X', 'Y')]) %>%
        as.matrix()
      
      min.id1 <- apply(dist1, 1, FUN = which.min)
      
      segment_CD8 <- data.frame(x = M_linkData$X, y = M_linkData$Y, xend = CD8_linkData[min.id1, 'X'], yend = CD8_linkData[min.id1, 'Y'])
      
      
      dist2 <- flexclust::dist2(M_linkData[, c('X', 'Y')], CD4_linkData[, c('X', 'Y')]) %>%
        as.matrix()
      
      min.id2 <- apply(dist2, 1, FUN = which.min)
      
      
      segment_CD4 <- data.frame(x = M_linkData$X, y = M_linkData$Y, xend = CD4_linkData[min.id2, 'X'], yend = CD4_linkData[min.id2, 'Y'])
      
      
      
      
      #------------------------------------------------#
      #------ Plot the link map long versus short -----#
      #------------------------------------------------#
      
      
      p <- ggplot() +
        geom_segment(data = segment_CD8, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_segment(data = segment_CD4, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_point(data = M_linkData, aes(X, Y, color = class), size = 4) +
        geom_point(data = CD8_linkData[min.id1,], aes(X, Y, color = class), size = 3) +
        geom_point(data = CD4_linkData[min.id2,], aes(X, Y, color = class), size = 3) +
        theme_void() +
        scale_color_manual(values = c('Macrophages' = '#c074eb', 'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0')) +
        theme(#axis.text = element_blank(),
          axis.title = element_blank(),
          #axis.ticks = element_blank(),
          legend.position = 'none') 
      p
      
      ggsave(p, file=paste0("Figures/validation_Links/", Sample, "_", ROI, "_Links_.jpg"), width = 6, height = 6, units = "in", dpi = 300)
      
    }
    
  }
}




#---------------------------------------#
#---- Tertiary Lymphoid Structure ------#
#---------------------------------------#


#-------------- HDBSCAN method to find clusters ----------------#
# create some empty vectors
PDACdata <- readRDS('./exportedFiles/PDACdata_USETHIS_validation.rds')

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
            plot(p)
            ggsave(p, file=paste0('Figures/TLS/TLS_networks_validation/', Sample, "_", ROI, "_", cl, ".png"), width = 10, height = 10, units = "in", dpi = 300)
            
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
#write.csv(clus.info, 'clus.info_validation.csv')







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


png('Figures/Clustering_ComplexHeatmap_validation.jpeg', width=22, height=35, units = 'cm', res = 300) 


foo2 = as.matrix(B_CD8T.clus.info[, 19:20])
row_ha = rowAnnotation('Density, mm' = row_anno_barplot(foo2, axis = TRUE, axis_param = list(side = "top"),
                                                        gp = gpar(fill = c("#82b4e5", "#4875b0"))), width = unit(4, "cm"))

lgd = Legend(at = c("CD8 T", "B cell"), title = "Cell type", legend_gp = gpar(fill = c("#4875b0", "#82b4e5")))


q <- Heatmap(as.matrix(B_CD8T.clus.info[, c(9:15, 26:32)]), name="ratio",
             col = colorRamp2(c(0, 1), c("white", "#377EB8")),
             width = unit(8, "cm"), 
             height = unit(20, "cm"), 
             row_km = 1, column_km = 1,
             border = TRUE, 
             row_dend_width = unit(4, 'cm'),
             column_dend_height = unit(4, 'cm'),
             row_title = NULL,
             column_title = NULL,
             heatmap_legend_param = list(legend_height = unit(4, "cm")),
             cluster_row_slices = FALSE
             #left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c('#4986c3', 'black')),
            #                                                  labels = c("Cluster 1", "Cluster 2"), 
            #                                                  labels_gp = gpar(col = "white", fontsize = 10))),             #cell_fun = TextFunc(pvalMatrices.vec, numdat = F)
) +
  Heatmap(B_CD8T.clus.info$label, name = "Cluster shape", width = unit(5, "mm"), col = c("circular" = "#a2b1d2", "elongated" = "#76b666", 'else' = '#f0a47c')) +
  #Heatmap(B_CD8T.clus.info$group, name = "Survival gorup", width = unit(5, "mm"), col = c("short" = "#c22828", "long" = "#325698")) +
  row_ha 
q

draw(q, heatmap_legend_list = list(lgd))
dev.off()




#-------------------------------------------#




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






#---------------------------------------------#
# Long- versus short term survival -----------#
#---------------------------------------------#

cohort <- read.csv('./exportedFiles/sample_key_10.23.csv', encoding = 'UTF-8-BE') %>%
  filter(Tumor.type == 'PDAC') %>%
  filter(Rx.type == 'naive') %>%
  filter(OS_days != 'NA') 
colnames(cohort)[1] <- 'sample'
colnames(cohort)[7] <- 'DFS_days'


cohort_valid <- read.csv('./exportedFiles/sample_key_10.23.csv', encoding = 'UTF-8-BE') %>%
  filter(Tumor.type == 'PDAC') %>%
  filter(Rx.type == 'NeoAdj') %>%
  filter(OS_days != 'NA') 
colnames(cohort_valid)[1] <- 'sample'
colnames(cohort_valid)[7] <- 'DFS_days'

# split to long versus short

thresh <- quantile(cohort$OS_days, c(0.5))

cohort$label <-  1* ifelse(cohort$OS_days > thresh, 1, 0)


# combine experimental dataset and validation dataset
#write.csv(cohort, 'cohort_short_long.csv')
km_fit <- survfit(Surv(DFS_days, death) ~ label, data=cohort)

ggsurvplot(
  fit = km_fit, 
  xlab = "Days", 
  ylab = "Overall survival probability",
  risk.table = TRUE,
  pval = TRUE)



#---------------------------------------------#
# Correlation analysis: Cell Density----------#
#---------------------------------------------#
PDACdata_valid <- readRDS('./exportedFiles/PDACdata_USETHIS_validation.rds')
PDACdata_explore <- readRDS('./exportedFiles/PDACdata_USETHIS.rds')


# Cell Density

# create some empty vectors
density_subregion_sample <- c()
each_density_subregion_sample <- c()

allFiles <- list.files('./Dataset_2021/ClassifiedCSV_10012021T8') 


for(sample in long.group){
  
  #sample <- 2009
  subsets <- grep(sample, allFiles, value = TRUE)
  
  # number of regions
  nregion <- length(subsets)
  
  # iterate through each subregion
  density_subregion <- c()
  each_density_subregion <- c()
  for(subregion in subsets){
    
    #ROIs
    ROIs <- rev(strsplit(subregion, "\\.|\\,|_")[[1]])[2]
    
    #subregion <- subsets[6]
    # single-cell dataset of subregion, immune cell coutns
    dat <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class != 'Noise') %>%
      filter(class != 'Other Cells') %>%
      filter(class != 'Tumor')
    
    
    
    
    # diffrent kind of cell density
    B.cells <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'Bcells') %>%
      nrow()
    
    Macrophages <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'Macrophages') %>%
      nrow()
    
    CD8 <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'CD8 Tcells') %>%
      nrow()
    
    CD4 <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'CD4 Tcells') %>%
      nrow()
    
    CD45_Other <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'CD45 Other') %>%
      nrow() 
    
    Tumor <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'Tumor') %>%
      nrow() 
    # construct bounding box
    
    Rect_area <- (  max(dat$X) - min(dat$X) * 0.5 ) * (max(dat$Y) - min(dat$Y) * 0.5) / 1000000
    
    # combine all densities
    density_subregion <- rbind.data.frame(density_subregion, cbind.data.frame( nrow(dat) / Rect_area, ROIs))
    
    # combine all densities for each cell type
    #dens.cells <- no.cells %>%
    #  mutate(density = n / Rect_area)
    
    each_density_subregion <- rbind.data.frame(each_density_subregion, cbind.data.frame(Tumor / Rect_area, B.cells / Rect_area, Macrophages / Rect_area, CD8 / Rect_area,CD4 / Rect_area, CD45_Other / Rect_area, ROIs))
    
  }
  
  density_subregion_sample <- rbind.data.frame(density_subregion_sample, cbind.data.frame(density_subregion, sample))
  each_density_subregion_sample <- rbind.data.frame(each_density_subregion_sample, cbind.data.frame(each_density_subregion, sample))
  
}



#---------- Marker Expression Validation -------------#

funcvars = c("PD1", "PDL1", 'Ki67', 'EOMES',  "IL10", 'ICOS', 'GZMB')


exprShort <- PDACdata_valid %>%
  filter(class != 'Other Cells' | class != 'Tumor') %>%
  select(29:36) %>%
  group_by(class) %>% 
  `colnames<-`(c('class', 'PD1', 'PDL1', 'Ki67', 'EOMES', 'GZMB', 'IL10', 'ICOS')) %>%
  summarise_each_(funs(sum(.==1,na.rm=TRUE)), funcvars) %>% 
  setNames(c(names(.)[1], paste0(funcvars))) %>%
  filter(class != 'Other Cells' & class != 'Tumor' & class != 'CD45 Other')



#----- plot, stacked bar plots
exprLong.stacked <- melt(exprShort)
p <- ggplot(exprLong.stacked, aes(fill=variable, y=value, x=as.factor(class))) + 
  scale_fill_manual(values = c('PDL1' = '#d16c6c', 'PD1' = '#71c8d7', 'IL10' = '#4986c3',
                               'Ki67' = '#2c955a', 'EOMES' = '#Bcb8b0', 'ICOS' = '#c6bae4', 'GZMB' = '#eab676'))  +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 20),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)
        #legend.position = 'NA'
  ) +
  scale_x_discrete(labels= c('B cells', 'CD4+ T', 'CD8+ T', 'Mpg'))

#scale_y_discrete(position = "right")
p

ggsave(p, file=paste0("Figures/Stacked_Functional_Markers_all_valid.jpg"), width = 6, height = 6, units = "in", dpi = 300)



#-- Validation

density_subregion_sample_valid <- c()
each_density_subregion_sample_valid <- c()

for(sample in cohort_valid$sample){
  
  #sample <- 2009
  subsets <- grep(sample, allFiles, value = TRUE)
  
  # number of regions
  nregion <- length(subsets)
  
  # iterate through each subregion
  density_subregion <- c()
  each_density_subregion <- c()
  for(subregion in subsets){
    
    #ROIs
    ROIs <- rev(strsplit(subregion, "\\.|\\,|_")[[1]])[2]
    
    #subregion <- subsets[6]
    # single-cell dataset of subregion, immune cell coutns
    dat <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class != 'Noise') %>%
      filter(class != 'Other Cells') %>%
      filter(class != 'Tumor')
    
    
    
    
    # diffrent kind of cell density
    B.cells <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'Bcells') %>%
      nrow()
    
    Macrophages <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'Macrophages') %>%
      nrow()
    
    CD8 <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'CD8 Tcells') %>%
      nrow()
    
    CD4 <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'CD4 Tcells') %>%
      nrow()
    
    CD45_Other <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'CD45 Other') %>%
      nrow() 
    
    Tumor <- read.csv(paste0('./Dataset_2021/Classified2021_CSV/', subregion), row.names = 1) %>%
      filter(class == 'Tumor') %>%
      nrow() 
    # construct bounding box
    
    Rect_area <- (  max(dat$X) - min(dat$X) * 0.5 ) * (max(dat$Y) - min(dat$Y) * 0.5) / 1000000
    
    # combine all densities
    density_subregion <- rbind.data.frame(density_subregion, cbind.data.frame( nrow(dat) / Rect_area, ROIs))
    
    # combine all densities for each cell type
    #dens.cells <- no.cells %>%
    #  mutate(density = n / Rect_area)
    
    each_density_subregion <- rbind.data.frame(each_density_subregion, cbind.data.frame(Tumor / Rect_area, B.cells / Rect_area, Macrophages / Rect_area, CD8 / Rect_area,CD4 / Rect_area, CD45_Other / Rect_area, ROIs))
    
  }
  
  density_subregion_sample_valid <- rbind.data.frame(density_subregion_sample_valid, cbind.data.frame(density_subregion, sample))
  each_density_subregion_sample_valid <- rbind.data.frame(each_density_subregion_sample_valid, cbind.data.frame(each_density_subregion, sample))
  
}






#---------------------------------------#
#--------- Density Profile -------------#
#---------------------------------------#

density_sample_mean <- each_density_subregion_sample_valid %>%
  group_by(sample) %>%
  dplyr::summarize(meanDensity_CD4 = mean(`CD4/Rect_area`),
                   meanDensity_CD8 = mean(`CD8/Rect_area`),
                   meanDensity_B = mean(`B.cells/Rect_area`),
                   meanDensity_Macrophage = mean(`Macrophages/Rect_area`))


write.csv(density_sample_mean, 'validation_cohort_cell_density.csv')

