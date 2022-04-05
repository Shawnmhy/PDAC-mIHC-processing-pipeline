#--------------------------------------#
# PDAC IHC processing pipeline 1-------#
# @Author: Haoyang Mi -----------------#
# Date: Aug 9th 2021 ------------------#

library(ggplot2); library(dplyr); library(reshape2); library(RColorBrewer); library(matrixStats)
library(ggvoronoi); library(tripack);library(igraph);library(concaveman);library(sf); library(spatstat); library(caramellar)
library(survival); library(survminer); library(ComplexHeatmap); library(circlize)
library(pdist); library(SNFtool); library(flexclust); library(gplots); library(gdata)
library(foreach); library(doParallel); library(caret); library(umap);
library(rstatix); library(ggpubr)
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
  dplyr::filter(Tumor.type == 'PDAC') %>%
  dplyr::filter(Rx.type == 'naive') %>%
  dplyr::filter(OS_days != 'NA') 
colnames(cohort)[1] <- 'sample'
colnames(cohort)[7] <- 'DFS_days'

# split to long versus short

thresh <- quantile(cohort$OS_days, c(0.5))
cohort$label <-  1* ifelse(cohort$OS_days > thresh, 1, 0)

#
long.group <- cohort %>%
  dplyr::filter(label == 1) %>%
  select(sample) %>%
  as.matrix()


short.group <- cohort %>%
  dplyr::filter(label == 0) %>%
  select(sample) %>%
  as.matrix()



#-----------------------------------------------------#
#------------ Expression Data Summary ----------------#
#-----------------------------------------------------#


#--- LONG TERM SURVIVAL GROUP
# Case 1: PD1
funcvars = c("PD1", "PDL1", 'Ki67', 'EOMES',  "IL10", 'ICOS')

exprLong <- PDACdata %>%
  dplyr::filter(sample %in% long.group) %>%
  select(29:36) %>%
  group_by(class) %>% 
  `colnames<-`(c('class', 'PD1', 'PDL1', 'Ki67', 'EOMES', 'GZMB', 'IL10', 'ICOS')) %>%
  summarise_each_(funs(sum(.==1,na.rm=TRUE)), funcvars) %>% 
  setNames(c(names(.)[1], paste0(funcvars))) %>%
  dplyr::filter(class != 'Other Cells' & class != 'Tumor' & class != 'CD45 Other')


exprShort <- PDACdata %>%
  dplyr::filter(sample %in% short.group) %>%
  dplyr::filter(class != 'Other Cells' | class != 'Tumor') %>%
  select(29:36) %>%
  group_by(class) %>% 
  `colnames<-`(c('class', 'PD1', 'PDL1', 'Ki67', 'EOMES', 'GZMB', 'IL10', 'ICOS')) %>%
  summarise_each_(funs(sum(.==1,na.rm=TRUE)), funcvars) %>% 
  setNames(c(names(.)[1], paste0(funcvars))) %>%
  dplyr::filter(class != 'Other Cells' & class != 'Tumor' & class != 'CD45 Other')



#----- plot, stacked bar plots
exprLong.stacked <- melt(exprLong)
p <- ggplot(exprLong.stacked, aes(fill=variable, y=value, x=as.factor(class))) + 
  scale_fill_manual(values = c('PDL1' = '#d16c6c', 'PD1' = '#71c8d7', 'IL10' = '#4986c3',
                               'Ki67' = '#2c955a', 'EOMES' = '#Bcb8b0', 'ICOS' = '#c6bae4'))  +
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
ggsave(p, file=paste0("Figures/Stacked_Functional_Markers_Long.jpg"), width = 6, height = 6, units = "in", dpi = 300)



exprShort.stacked <- melt(exprShort)

p <- ggplot(exprShort.stacked, aes(fill=variable, y=value, x=as.factor(class))) + 
  scale_fill_manual(values = c('PDL1' = '#d16c6c', 'PD1' = '#71c8d7', 'IL10' = '#4986c3',
                               'Ki67' = '#2c955a', 'EOMES' = '#Bcb8b0', 'ICOS' = '#c6bae4'))  +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 20),
        axis.ticks = element_blank(),
        legend.position = 'NA'
  ) +
  scale_x_discrete(labels= c('B cells', 'CD4+ T', 'CD8+ T', 'Mpg'))
#scale_y_discrete(position = "right")
p
ggsave(p, file=paste0("Figures/Stacked_Functional_Markers_Short.jpg"), width = 5, height = 6, units = "in", dpi = 300)






#-------------------------------------------------------------------------#
# Compare the overall expression of markers on different cell types ------#
#-------------------------------------------------------------------------#


boxdf <- PDACdata %>%
  select(c('PD1', 'PDL1', 'IL10', 'ICOS', 'KI67', 'EOMES', 'class', 'sample')) %>%
  dplyr::filter(class != 'Other Cells') %>%
  dplyr::filter(class != 'Tumor') %>%
  dplyr::filter(class != 'CD45 Other')
  
boxdf$label <-  1* ifelse(boxdf$sample %in% long.group, 1, 0)
boxdf[boxdf$label == '1', 'label'] <- 'long-term'
boxdf[boxdf$label == '0', 'label'] <- 'short-term'

# remove sample
boxdf <- boxdf %>%
  select(-'sample')


boxdf <- melt(boxdf, id.vars = c('class', 'label'))
boxdf$class <- as.factor(boxdf$class)


#anno_df <- compare_means(value ~ label, group.by = "class", data = boxdf_melt, p.adjust.method = "holm") %>%
#  mutate(y_pos = 40, p.adj = format.pval(p.adj, digits = 2))
boxdf[boxdf$label == '1', 'label'] <- 'long-term'
boxdf[boxdf$label == '0', 'label'] <- 'short-term'

stat.test <- boxdf %>%
  group_by(variable, class) %>%
  wilcox_test(value ~ label) %>%
  adjust_pvalue(method = 'fdr') %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = 'label')

jpeg('./Figures/Functional_Marker_Profile_Outcome.jpeg', units = 'in', width = 12, height =16, res = 300)
ggboxplot(boxdf, x = 'label', y = 'value', palette = 'lancet', color = 'label', facet.by = c('variable', 'class')) +
  #facet_grid(class ~ variable) +
  geom_boxplot(data=boxdf, mapping=aes(x=label, y= value, color=label)) +
  stat_pvalue_manual(stat.test, label = 'p.adj.signif', y.position = max(boxdf$value) + 0.05, label.size = 6) +
  ylim(0, 1) +
  ylab('Expression') +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = 16)) 
dev.off()






#----------------------------#
funcvars = c("PD1", "PDL1", 'Ki67', 'EOMES', 'GZMB', "IL10", 'ICOS')

exprAll <- PDACdata %>%
  select(29:36) %>%
  group_by(class) %>% 
  `colnames<-`(c('class', 'PD1', 'PDL1', 'Ki67', 'EOMES', 'GZMB', 'IL10', 'ICOS')) %>%
  summarise_each_(funs(sum(.==1,na.rm=TRUE)), funcvars) %>% 
  setNames(c(names(.)[1], paste0(funcvars))) %>%
  filter(class != 'Other Cells' & class != 'Tumor' & class != 'CD45 Other')



exprAll.stacked <- melt(exprAll)

p <- ggplot(exprAll.stacked, aes(fill=variable, y=value, x=as.factor(class))) + 
  scale_fill_manual(values = c('PDL1' = '#d16c6c', 'PD1' = '#71c8d7', 'IL10' = '#4986c3', 'GZMB' = '#eab676',
                               'Ki67' = '#2c955a', 'EOMES' = '#Bcb8b0', 'ICOS' = '#c6bae4'))  +
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
ggsave(p, file=paste0("Figures/Stacked_Functional_Markers_all.jpg"), width = 6, height = 6, units = "in", dpi = 300)




#------------------------------------------#
#-------- Simulation analysis -------------# 
#------------------------------------------#
IRPvars <- c("PD1.", "PDL1.", 'KI67.', 'EOMES.', 'GZMB.',  "IL10.", 'ICOS.')

es.Mpg_All <- data.frame(matrix(nrow = 0, ncol = 0))
es.CD4T_All <- data.frame(matrix(nrow = 0, ncol = 0))

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
    
    CD4T <- roiDat %>%
      filter(class == 'CD4 Tcells')
    
    Mpg <- roiDat %>%
      filter(class == 'Macrophages')
    
    Rect <- data.frame(cbind(X = c(min(roiDat$X), max(roiDat$X), max(roiDat$X), min(roiDat$X)), Y = c(min(roiDat$Y), min(roiDat$Y), max(roiDat$Y), max(roiDat$Y))))
    
    
    if(nrow(CD4T)  > 10 & nrow(Mpg) > 10){
      
      
      
      ### Case 1: CD4 T cell to Macrophage
      distMatrix <- flexclust::dist2(CD4T[, c('X', 'Y')], Mpg[, c('X', 'Y')]) %>%
        as.matrix()
      
      # These CD4 T cells has at least 1 adjacent Macrophages
      adj.CD4T.id <- which(rowMins(distMatrix) < 100)
      
      # Get the single-cell data for these CD4 T cells
      adj.CD4T <- CD4T[adj.CD4T.id, ]
      
      if(!(is.empty(adj.CD4T.id))){
        
        # Real ratio
        realRatio_CD4T <- adj.CD4T %>%
          summarise_each_(funs(sum(.==1,na.rm=TRUE)), IRPvars) %>%
          `/` (nrow(adj.CD4T))
        
        
      }
      
      
      ### Case 2: Macrophage to CD4 T
      
      distMatrix <- flexclust::dist2(Mpg[, c('X', 'Y')], CD4T[, c('X', 'Y')]) %>%
        as.matrix()
      # These macrophages have at least 1 adjacent CD4T
      adj.Mpg.id <- which(rowMins(distMatrix) < 100)
      
      # Get the single-cell data for these Macrophages
      adj.Mpg <- Mpg[adj.Mpg.id, ]
      
      
      
      
      if(!(is.empty(adj.Mpg.id))){
        # Real ratio
        realRatio_Mpg <- adj.Mpg %>%
          summarise_each_(funs(sum(.==1,na.rm=TRUE)), IRPvars) %>%
          `/` (nrow(adj.Mpg))
        
      }
      
      
      
      # simulated ratios
      simRatios <- ratio.sim(CD4T, Mpg, Rect, 500)
      
      # enrichment score for CD4 T cells clusters
      es.CD4T_all <- data.frame(matrix(nrow = 1, ncol = 0))
      for(id in seq_len(7)){
        
        #id <- 1
        # enrichment score
        es.CD4T <- (length(which(realRatio_CD4T[,id] > simRatios[[1]][,id])) + 500 - nrow(simRatios[[1]]))/ 500
        
        es.CD4T_all <- cbind(es.CD4T_all, es.CD4T)
      }
      
      es.CD4T_all <- es.CD4T_all %>%
        cbind.data.frame(ROI, Sample) %>%
        `colnames<-` (c(colnames(realRatio_CD4T), 'ROI', 'Sample'))
      
      es.CD4T_All <- rbind.data.frame(es.CD4T_All, es.CD4T_all)
      
      # enrichment score for CD4 T cells clusters
      es.Mpg_all <- data.frame(matrix(nrow = 1, ncol = 0))
      for(id in seq_len(7)){
        
        #id <- 1
        
        # enrichment score
        es.Mpg <- (length(which(realRatio_Mpg[,id] > simRatios[[2]][,id])) + 500 - nrow(simRatios[[2]]))/ 500
        
        es.Mpg_all <- cbind.data.frame(es.Mpg_all, es.Mpg)
      }
      
      es.Mpg_all <- es.Mpg_all %>%
        cbind.data.frame(ROI, Sample) %>%
        `colnames<-` (c(colnames(realRatio_Mpg), 'ROI', 'Sample'))
      
      es.Mpg_All <- rbind.data.frame(es.Mpg_All, es.Mpg_all)
      
      print(Sample)
    }      
  }
}





mean(es.Mpg_All$PD1.)
mean(es.Mpg_All$PDL1.)
mean(es.Mpg_All$KI67.)
mean(es.Mpg_All$EOMES.)
mean(es.Mpg_All$IL10.)
mean(es.Mpg_All$ICOS.)
mean(es.Mpg_All$GZMB.)

boxplot(es.Mpg_All$PD1., es.Mpg_All$PDL1., es.Mpg_All$KI67., es.Mpg_All$EOMES.,  es.Mpg_All$IL10., es.Mpg_All$ICOS.)
####
mean(es.CD4T_All$PD1.)

mean(es.CD4T_All$PDL1.)
mean(es.CD4T_All$KI67.)
mean(es.CD4T_All$EOMES.)
mean(es.CD4T_All$IL10.)
mean(es.CD4T_All$ICOS.)
mean(es.CD4T_All$PD1.)



es.CD4T_All_BOX <- melt(es.CD4T_All[, 1:7])
es.Mpg_All_BOX <- melt(es.Mpg_All[, 1:7])

es.CD4T_All_BOX$variable <- factor(es.CD4T_All_BOX$variable, levels = c('PDL1.', 'PD1.', 'IL10.', 'KI67.', 'EOMES.', 'GZMB.','ICOS.'))

es.Mpg_All_BOX$variable <- factor(es.Mpg_All_BOX$variable, levels = c('PDL1.', 'PD1.', 'IL10.', 'KI67.', 'EOMES.', 'GZMB.', 'ICOS.'))
p <- ggplot(es.CD4T_All_BOX, aes(x = variable, y = value, color = variable)) +
  theme_bw() +
  stat_boxplot( aes(variable, value), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot(size = 1) +
  #geom_violin(size = 1) +
  #geom_path(data = df_for_line, aes(x = Region, y = mean_y, group = 1), size = 1, color = 'black') +
  #geom_point(data = df_for_line, aes(x = Region, y = mean_y, group = 1), size = 4) +
  scale_color_manual(values = c('PDL1.' = '#d16c6c', 'PD1.' = '#71c8d7', 'IL10.' = '#4986c3', 'GZMB.' = '#eab676',
                                'KI67.' = '#2c955a', 'EOMES.' = '#Bcb8b0', 'ICOS.' = '#c6bae4'))  +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = 'none') 

p

ggsave(p, file=paste0("Figures/Simulation_test_Mpg.jpg"), width = 4, height = 6, units = "in", dpi = 300)




#---------------------------------------#
#---- CD8+ T cell effective score ------#
#---------------------------------------#
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



# statistics
min(eff.scores[eff.scores$class == 'long', 'eff.score'])
mean(eff.scores[eff.scores$class == 'long', 'eff.score'])
mean(eff.scores[eff.scores$class == 'short', 'eff.score'])


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
  
  Sample <- 5213
  # sample data
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  # get all ROIs
  ROIs <- unique(sampleDat$ROIs)
  
  # short or long
  label <- 'short'
  if(Sample %in% long.group){
    label <- 'long'
  }
  
  for(ROI in ROIs){
    ROI <- 'ROI06'
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
        #theme_void() +
        scale_color_manual(values = c('Macrophages' = '#c074eb', 'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0')) +
        theme(#axis.text = element_blank(),
          axis.title = element_blank(),
          #axis.ticks = element_blank(),
          legend.position = 'none') 
      p
      
      ggsave(p, file=paste0("Figures/Links/", label, '_', Sample, "_", ROI, "_Links_.jpg"), width = 6, height = 6, units = "in", dpi = 300)
      
    }
   
  }
}



p <- ggplot() +
  geom_segment(data = segment_CD8, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_segment(data = segment_CD4, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = M_linkData, aes(X, Y, color = class), size = 4) +
  geom_point(data = CD8_linkData[min.id1,], aes(X, Y, color = class), size = 3) +
  geom_point(data = CD4_linkData[min.id2,], aes(X, Y, color = class), size = 3) +
  theme_bw() +
  scale_color_manual(values = c('Macrophages' = '#c074eb', 'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0')) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 
p

ggsave(p, file=paste0("Figures/Links_short.jpg"), width = 6, height = 6, units = "in", dpi = 300)
