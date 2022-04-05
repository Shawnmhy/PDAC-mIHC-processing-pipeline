# PDAC Analysis Part 0
# @Author: Haoyang Mi



library(ggplot2)
library(stringr)
library(ggvoronoi)
library(RANN)
library(dplyr); library(tidyr); library(reshape2)
library(rjson)
library(matrixStats)
library(umap)
library(spatstat)
library(tripack); library(pracma)
library(Matrix); library(ggrepel)
library(foreach); library(doParallel)
library(survival); library(survminer); library(ComplexHeatmap); library(circlize)
library(doParallel); library(foreach)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd('..')
source("./Codes/Function.r")

# how many ROIs per sample

PDACdata <- data.frame(matrix(nrow = 0, ncol = 0))
cohort <- read.csv('./exportedFiles/sample_key_10.23.csv', encoding = 'UTF-8-BE') %>%
  filter(Tumor.type == 'PDAC') %>%
  filter(Rx.type == 'naive') %>%
  filter(OS_days != 'NA')
colnames(cohort)[1] <- 'sample'
colnames(cohort)[7] <- 'DFS_days'

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

saveRDS(dat_regions_samples, 'PDACdata_USETHIS.rds')

#----------Heatmaps---------------#
# scale using quantile

cofactor <- 5
dat_regions_samples[,1:24] <- asinh(dat_regions_samples[,1:24]/cofactor)

rng <- colQuantiles(as.matrix(dat_regions_samples[,1:24]), probs = c(0.01, 0.99))


dat_regions_samples0 <- t((t(dat_regions_samples[,1:24]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
dat_regions_samples0[dat_regions_samples0 < 0] <- 0; dat_regions_samples0[dat_regions_samples0 > 1] <- 1

dat_regions_samples[,1:24] <- dat_regions_samples0

#saveRDS(dat_regions_samples, 'PDACdata_arcsinh_scaled.rds')


#-------------umaps--------------#

patient_2009_expr <- dat_regions_samples %>%
  filter(sample == '2170')

umap.pdac <- umap(patient_2009_expr[,1:24])

umap_data <- cbind.data.frame(umap.pdac$layout, patient_2009_expr$class)
colnames(umap_data) <- c('X', 'Y', 'class')

p <- ggplot(umap_data) +
  theme_void() +
  theme(legend.position = 'NA') +
  geom_point(aes(X, Y, color = class), size = 0.5) +
  scale_color_manual(values = c('Tumor' = '#fb8648', 'Macrophages' = '#c074eb', 'Bcells' = '#82b4e5',
                                'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0', 'Other Cells' = 'grey', 'CD45 Other' = '#ebe459'))  


ggsave(p, file=paste0("Figures/UMAP_2170.png"), width = 5, height = 5, units = "in", dpi = 300)

#------ end here --------#



# sort by immune cell abundance
PDACdata <- readRDS('./exportedFiles/PDACdata_USETHIS.rds')

Immune_count <- PDACdata %>%
  filter(class != 'Tumor') %>%
  filter(class != 'Other Cells') %>%
  group_by(ROIs, sample) %>%
  tally()



cohort <- read.csv('./sample_key_10.23.csv', encoding = 'UTF-8-BE') %>%
  filter(Tumor.type == 'PDAC') %>%
  filter(Rx.type == 'naive') %>%
  filter(OS_days != 'NA')
colnames(cohort)[1] <- 'sample'
colnames(cohort)[7] <- 'DFS_days'

sampleID <- cohort$sample

# read single-cell dataset
allFiles <- list.files('./Dataset_2021/Classified2021_CSV') 



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
      #theme(legend.position = 'NA') +
      geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +
      #scale_fill_manual(values = c('Tumor' = '#fb8648', 'Macrophages' = '#c074eb', 'Bcells' = '#82b4e5',
      #                             'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#ebe459'))  +
      
      coord_fixed(ratio = 1) 
    p
    ggsave(p, file=paste0("Figures/Voronoi_Tessellation/", sample,'_', ROIs, '.png'), width = 8, height = 6, units = "in", dpi = 300)
    
    
    
  }
}


# -------- sort immune cell counts by patient --------#


# create some empty vectors
density_subregion_sample <- c()
each_density_subregion_sample <- c()

for(sample in sampleID){
  
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


pdf('Survival/inflammation_profile_short.pdf', width=5, height=20) 

q <- Heatmap(as.matrix(each_density_subregion_sample[each_density_subregion_sample$sample %in% long.group, -c(7:8)]), name="rho",
             col = colorRamp2(c(0, 300, 600), c("white", "#377EB8", "#E41A1C")),
             #width = unit(10, "cm"), 
             #height = unit(10, "cm"), 
             row_dend_width = unit(2, 'cm'),
             column_dend_height = unit(2, 'cm'),
             #row_title = NULL,
             #column_title = NULL
             heatmap_legend_param = list(legend_height = unit(4, "cm")),
             column_names_gp = grid::gpar(fontsize = 25),
             row_names_gp = grid::gpar(fontsize = 25),
             #show_heatmap_legend = FALSE,
             #cell_fun = TextFunc(pvalMatrices.vec, numdat = F)
) 

q
dev.off()

# group by sample
colnames(density_subregion_sample)[1] <- 'density'

density_sample_mean <- density_subregion_sample %>%
  group_by(sample) %>%
  summarize(meanDensity = mean(density))


# sort patient by immune cell counts
density_sample_mean <- density_subregion_sample %>%
  group_by(sample) %>%
  summarize(meanDensity = mean(density))

density_sample_mean$sample <- as.factor(density_sample_mean$sample)


#----------------- Chi-square test for immune infiltration level -----------------#


density_sample_mean$sample = with(density_sample_mean, reorder(sample, meanDensity))
density_sample_mean = density_sample_mean[order(density_sample_mean$sample), ]


# thresholding
thresh <- quantile(density_sample_mean$meanDensity, c(0.33, 0.67))

density_sample_mean$label <- ifelse(density_sample_mean$sample %in% long.group, 1, 0)

density_sample_mean$infil <- ifelse(density_sample_mean$meanDensity > thresh[1]  & density_sample_mean$meanDensity <= thresh[2] , 1, 2)
density_sample_mean[density_sample_mean$meanDensity < thresh[1], 'infil'] <- 0

#write.csv(density_sample_mean, 'density_infil_surv.csv')
chisq.test(as.matrix(table(density_sample_mean[, 3:4])))
# plot sorted bar plot
p <-ggplot(data=density_sample_mean, aes(y= reorder(sample, -meanDensity), x=meanDensity)) +
  geom_bar(stat="identity", fill="grey")+
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 24),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = ifelse(density_sample_mean$sample %in% long.group, '#c22828', '#355794')))
p


ggsave(p, file=paste0("Figures/Immune_density_sorted_low_to_high.png"), width = 4, height = 20, units = "in", dpi = 300)



##### stacked bar plot showing immune components

# 4399 patient
PDACdata_stacked <- PDACdata %>%
  group_by(class,sample) %>%
  filter(class != 'Other Cells') %>%
  filter(class != 'Tumor') %>%
  tally()

# reorder sample according to immune density
PDACdata_stacked$sample <- as.factor(PDACdata_stacked$sample)

# use what order to reorder stacked plot?

orders <- density_sample_mean %>% arrange(desc(meanDensity))
PDACdata_stacked$sample <- factor(PDACdata_stacked$sample, levels = orders$sample)




p <- ggplot(PDACdata_stacked, aes(fill=class, x=n, y=as.factor(sample))) + 
  scale_fill_manual(values = c('Tumor' = '#fb8648', 'Macrophages' = '#c074eb', 'Bcells' = '#82b4e5',
                               'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#ebe459'))  +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.ticks = element_blank(),
        legend.position = 'NA')
p
ggsave(p, file=paste0("Figures/Immune_composition_sorted_low_to_high.png"), width = 4, height = 20, units = "in", dpi = 300)


##### stacked bar plot showing immune components

# get ROI stacked bar plot from one patient
PDACdata_sub_stacked <- PDACdata %>%
  filter(sample == '2243') %>%
  group_by(ROIs, class) %>%
  filter(class != 'Other Cells') %>%
  filter(class != 'Tumor') %>%
  tally()

p <- ggplot(PDACdata_sub_stacked, aes(fill=class, x=n, y=as.factor(ROIs))) + 
  scale_fill_manual(values = c('Tumor' = '#fb8648', 'Macrophages' = '#c074eb', 'Bcells' = '#82b4e5',
                               'CD4 Tcells' = '#619e75', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#ebe459'))  +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 24),
        axis.ticks = element_blank(),
        legend.position = 'NA') +
  scale_y_discrete(position = "right")
p
ggsave(p, file=paste0("Figures/2243_composition_sorted_low_to_high.png"), width = 4, height = 8, units = "in", dpi = 300)



#---------------------------------------------#
# Long- versus short term survival -----------#
#---------------------------------------------#

cohort <- read.csv('./exportedFiles/sample_key_10.23.csv', encoding = 'UTF-8-BE') %>%
  filter(Tumor.type == 'PDAC') %>%
  filter(Rx.type == 'naive') %>%
  filter(OS_days != 'NA') 
colnames(cohort)[1] <- 'sample'
colnames(cohort)[7] <- 'DFS_days'

# split to long versus short

thresh <- quantile(cohort$OS_days, c(0.5))
cohort$label <-  1* ifelse(cohort$OS_days > thresh, 1, 0)

#write.csv(cohort, 'cohort_short_long.csv')
km_fit <- survfit(Surv(DFS_days, death) ~ label, data=cohort)

ggsurvplot(
  fit = km_fit, 
  xlab = "Days", 
  ylab = "Overall survival probability",
  risk.table = TRUE,
  pval = TRUE)

median(cohort[cohort$label == 1, 'OS_days'])

#
long.group <- cohort %>%
  filter(label == 1) %>%
  select(sample) %>%
  as.matrix()


short.group <- cohort %>%
  filter(label == 0) %>%
  select(sample) %>%
  as.matrix()





#---------------- pairwise Gcross function ------------------#
Mutual.dist <- data.frame(matrix(nrow = 0, ncol = 0))
for(Sample in unique(PDACdata$sample)){
  
  #ID <- 2089
  
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  for(ROI in unique(sampleDat$ROIs)){
    
    # ROI data
    #ROI <- 'ROI03'
    ROIdat <- sampleDat %>%
      filter(ROIs == ROI)
    
    # B cells coords
    posB <- ROIdat %>%
      filter(class == 'Bcells')
    
    # CD4 cells coords
    posCD4 <- ROIdat %>%
      filter(class == 'CD4 Tcells')
    
    # CD8 cells coords
    posCD8 <- ROIdat %>%
      filter(class == 'CD8 Tcells')
    
    # Macrophage cells coords
    posMpg <- ROIdat %>%
      filter(class == 'Macrophages')
    
    # Tumor cell coords
    posTumor <- ROIdat %>%
      filter(class == 'Tumor')
    
    
    # define region to compute Kcross function
    # Kcross functions
    Rect <- data.frame(cbind(x = c(min(ROIdat$X), max(ROIdat$X), max(ROIdat$X), min(ROIdat$X)), y = c(min(ROIdat$Y), min(ROIdat$Y), max(ROIdat$Y), max(ROIdat$Y))))
    
    # CD4 + CD8
    list1 <- bivarAnalysis.Kcross('posCD4', 'posCD8', Region = list(Rect))
    
    # CD4 + B
    list2 <- bivarAnalysis.Kcross('posCD4', 'posB', Region = list(Rect))
    
    # CD4 + Macrophage
    list3 <- bivarAnalysis.Kcross('posCD4', 'posMpg', Region = list(Rect))
    
    # CD8 + B
    list4 <- bivarAnalysis.Kcross('posCD8', 'posB', Region = list(Rect))
    
    # CD8 + Macrophage
    list5 <- bivarAnalysis.Kcross('posCD8', 'posMpg', Region = list(Rect))
    
    # B + Macrophage
    list6 <- bivarAnalysis.Kcross('posB', 'posMpg', Region = list(Rect))
    
    # Tumor + CD4
    list7 <- bivarAnalysis.Kcross('posCD4', 'posTumor', Region = list(Rect))
    
    # Tumor + CD8
    list8 <- bivarAnalysis.Kcross('posCD8', 'posTumor', Region = list(Rect))
    
    # Tumor + Macrophage
    list9 <- bivarAnalysis.Kcross('posTumor', 'posMpg', Region = list(Rect))
    
    # Tumor + B
    list10 <- bivarAnalysis.Kcross('posB', 'posTumor', Region = list(Rect))
    
    Mutual.dist <- rbind.data.frame(Mutual.dist, cbind.data.frame(Sample, ROI, list1[[1]], list1[[2]],  list2[[1]], list2[[2]], list3[[1]], list3[[2]],
                                                                  list4[[1]], list4[[2]], list5[[1]], list5[[2]], list6[[1]], list6[[2]],
                                                                  list7[[1]], list7[[2]], list8[[1]], list8[[2]],  list9[[1]], list9[[2]], list10[[1]], list10[[2]]))
    
  }
  print(Sample)
  
}




q <- Heatmap(as.matrix(Mutual.dist[Mutual.dist$Sample %in% short.group, -c(1:2)]), name="rho",
             col = colorRamp2(c(0, 60), c("white", "#E41A1C")),
             #width = unit(10, "cm"), 
             #height = unit(10, "cm"), 
             row_dend_width = unit(2, 'cm'),
             column_dend_height = unit(2, 'cm'),
             #row_title = NULL,
             #column_title = NULL
             heatmap_legend_param = list(legend_height = unit(4, "cm")),
             column_names_gp = grid::gpar(fontsize = 25),
             row_names_gp = grid::gpar(fontsize = 25),
             #show_heatmap_legend = FALSE,
             #cell_fun = TextFunc(pvalMatrices.vec, numdat = F)
) 
q





#------------------------------------------------#
#--------------- Shannon Entropy ----------------#
prox_thresh <- 78
ShannonH_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(Sample in unique(PDACdata$sample)){
  
  #ID <- 2089
  
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  for(ROI in unique(sampleDat$ROIs)){
    
    
    #core_id <- 3
    # For each core, create a folder to store all network files
    
    #dir.create(file.path('D:/DP/Projects/HCC/data/Networks/', core_id), showWarnings = FALSE)
    
    ROIDat <- sampleDat[sampleDat$ROIs == ROI, ]
    
    # --------------------------------------------#
    #      compute spatial Shannon's Entropy      #
    # --------------------------------------------#
    
    tensor <- data.frame(table(ROIDat$class)) %>%
      filter(!(Var1 %in% c('Noise', 'Other Cells', 'Tumor'))) #%>%
      #filter(Freq > 100)
    
    ctype_tensor <- as.character(unique(tensor$Var1))
    
    # all cell entropy

    if(nrow(tensor) == 5){
      ShannonH.a <- ShannonE(ctype_tensor, ROIDat)
      ShannonH_all <- rbind.data.frame(ShannonH_all, cbind.data.frame(ShannonH.a, ROI, Sample))
    }
  }
}

ShannonH_all_sample <- ShannonH_all %>%
  group_by(Sample) %>%
  summarise(meanShannon = mean(ShannonH.a)) %>%
  as.data.frame()

ShannonH_all_sample$label <- ifelse(ShannonH_all_sample$Sample %in% long.group, 1, 0)

wilcox.test(ShannonH_all_sample[ShannonH_all_sample$Sample %in% long.group, 2],  ShannonH_all_sample[ShannonH_all_sample$Sample %in% short.group, 2])
write.csv(ShannonH_all_sample, 'ShannonH_all.csv')

#---------------- pairwise Gcross function ------------------#
list1_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list1_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list2_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list2_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list3_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list3_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list4_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list4_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list5_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list5_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list6_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list6_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list7_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list7_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list8_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list8_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list9_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list9_2.all <- data.frame(matrix(nrow = 0, ncol = 0))
list10_1.all <- data.frame(matrix(nrow = 0, ncol = 0))
list10_2.all <- data.frame(matrix(nrow = 0, ncol = 0))


PDACdata <- readRDS('./exportedFiles/PDACdata_USETHIS.rds')

for(Sample in unique(PDACdata$sample)){
  
  #ID <- 2089
  #Sample <- 2089
  sampleDat <- PDACdata %>%
    filter(sample == Sample)
  
  for(ROI in unique(sampleDat$ROIs)){
    
    # ROI data
    #ROI <- 'ROI03'
    ROIdat <- sampleDat %>%
      filter(ROIs == ROI)
    
    # B cells coords
    posB <- ROIdat %>%
      filter(class == 'Bcells')
    
    # CD4 cells coords
    posCD4 <- ROIdat %>%
      filter(class == 'CD4 Tcells')
    
    # CD8 cells coords
    posCD8 <- ROIdat %>%
      filter(class == 'CD8 Tcells')
    
    # Macrophage cells coords
    posMpg <- ROIdat %>%
      filter(class == 'Macrophages')
    
    # Tumor cell coords
    posTumor <- ROIdat %>%
      filter(class == 'Tumor')
    
    
    # define region to compute Kcross function
    # Kcross functions
    Rect <- data.frame(cbind(x = c(min(ROIdat$X), max(ROIdat$X), max(ROIdat$X), min(ROIdat$X)), y = c(min(ROIdat$Y), min(ROIdat$Y), max(ROIdat$Y), max(ROIdat$Y))))
    
    # CD4 + CD8
    if(nrow(posCD4) * nrow(posCD8) > 0){
      list1 <- gcross('posCD4', 'posCD8', Region = list(Rect))
      list1_1.all <- rbind.data.frame(list1_1.all, cbind.data.frame(list1[1], ROI, Sample))
      list1_2.all <- rbind.data.frame(list1_2.all, cbind.data.frame(list1[2], ROI, Sample))
    }
    
    # CD4 + B
    if(nrow(posCD4) * nrow(posB) > 0){
      
      list2 <- gcross('posCD4', 'posB', Region = list(Rect))
      list2_1.all <- rbind.data.frame(list2_1.all, cbind.data.frame(list2[1], ROI, Sample))
      list2_2.all <- rbind.data.frame(list2_2.all, cbind.data.frame(list2[2], ROI, Sample))
      
    }
    
    # CD4 + Macrophage
    if(nrow(posCD4) * nrow(posMpg) > 0){
      list3 <- gcross('posCD4', 'posMpg', Region = list(Rect))
      list3_1.all <- rbind.data.frame(list3_1.all, cbind.data.frame(list3[1], ROI, Sample))
      list3_2.all <- rbind.data.frame(list3_2.all, cbind.data.frame(list3[2], ROI, Sample)) 
    }
    
    # CD8 + B
    if(nrow(posCD8) * nrow(posB) > 0){
      
      list4 <- gcross('posCD8', 'posB', Region = list(Rect))
      list4_1.all <- rbind.data.frame(list4_1.all, cbind.data.frame(list4[1], ROI, Sample))
      list4_2.all <- rbind.data.frame(list4_2.all, cbind.data.frame(list4[2], ROI, Sample))   
    }
    
    # CD8 + Macrophage
    if(nrow(posCD8) * nrow(posMpg) > 0){
      
      list5 <- gcross('posCD8', 'posMpg', Region = list(Rect))
      list5_1.all <- rbind.data.frame(list5_1.all, cbind.data.frame(list5[1], ROI, Sample))
      list5_2.all <- rbind.data.frame(list5_2.all, cbind.data.frame(list5[2], ROI, Sample))   
    }
    
    # B + Macrophage
    if(nrow(posB) * nrow(posMpg) > 0){
      
      list6 <- gcross('posB', 'posMpg', Region = list(Rect))
      list6_1.all <- rbind.data.frame(list6_1.all, cbind.data.frame(list6[1], ROI, Sample))
      list6_2.all <- rbind.data.frame(list6_2.all, cbind.data.frame(list6[2], ROI, Sample)) 
    }
    
    # Tumor + CD4
    if(nrow(posCD4) * nrow(posTumor) > 0){
      
      list7 <- gcross('posCD4', 'posTumor', Region = list(Rect))
      list7_1.all <- rbind.data.frame(list7_1.all, cbind.data.frame(list7[1], ROI, Sample))
      list7_2.all <- rbind.data.frame(list7_2.all, cbind.data.frame(list7[2], ROI, Sample))   
    }
    
    # Tumor + CD8
    if(nrow(posCD8) * nrow(posTumor) > 0){
      
      list8 <- gcross('posCD8', 'posTumor', Region = list(Rect))
      list8_1.all <- rbind.data.frame(list8_1.all, cbind.data.frame(list8[1], ROI, Sample))
      list8_2.all <- rbind.data.frame(list8_2.all, cbind.data.frame(list8[2], ROI, Sample))    
    }
    
    # Tumor + Macrophage
    if(nrow(posTumor) * nrow(posMpg) > 0){
      list9 <- gcross('posTumor', 'posMpg', Region = list(Rect))
      list9_1.all <- rbind.data.frame(list9_1.all, cbind.data.frame(list9[1], ROI, Sample))
      list9_2.all <- rbind.data.frame(list9_2.all, cbind.data.frame(list9[2], ROI, Sample))
    }
    
    # Tumor + B
    if(nrow(posTumor) * nrow(posB) > 0){
      
      list10 <- gcross('posB', 'posTumor', Region = list(Rect))
      list10_1.all <- rbind.data.frame(list10_1.all, cbind.data.frame(list10[1], ROI, Sample))
      list10_2.all <- rbind.data.frame(list10_2.all, cbind.data.frame(list10[2], ROI, Sample))    }
    
    
  }
  print(Sample)
  
}


#----------------------- adjust for p values ---------------------#
test.pval <- c()

for(id in seq_len(10)){
  
  # list 
  #id <- 4
  
  # select cell pair
  
  # select 
  list.cur_1 <- get(eval(paste0('list', id, '_1.all')))
  list.cur_2 <- get(eval(paste0('list', id, '_2.all')))
  
  list.long <- list.cur_1 %>%
    filter(Sample %in% long.group) #%>%
    #filter(r <= 80 & r >= 60)
  #filter(r <= 20)
  
  list.short <- list.cur_1 %>%
    filter(Sample %in% short.group) #%>%
    #filter(r <= 80 & r >= 60)
  
  #filter(r <= 20)
  
  #
  test1 <- wilcox.test(list.short$rs, list.long$rs)$p.value
  
  test.pval <- c(test.pval, test1)
  
  
  list.long <- list.cur_2 %>%
    filter(Sample %in% long.group) %>%
    filter(r <= 100 & r >= 80)
  
  #filter(r <= 20)
  
  list.short <- list.cur_2 %>%
    filter(Sample %in% short.group) %>%
    filter(r <= 100 & r >= 80)
  
  
  #filter(r <= 20)
  test2 <- wilcox.test(list.short$rs, list.long$rs)$p.value
  
  test.pval <- c(test.pval, test2)
}

test <- list.long %>%
  group_by(ROI, Sample) %>%
  mutate(max = max(rs))
p.adjust(test.pval)


#----------------------Plot Spatial G-function -------------------------#
p <-ggplot() +
  theme_test() +
  geom_path(data = list.long[list.long$Sample == 4101 & list.long$ROI == 'ROI02',], aes(r/2, rs), color = '#325698', size = 1) +
  geom_path(data = list.short[list.short$Sample == 2009 & list.short$ROI == 'ROI04',], aes(r/2, rs), color = '#c22828', size = 1) +
  xlab(expression('radius,'~mu~'m')) +
  ylab('Gcross function value') +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))
p
ggsave(p, file=paste0("Figures/G-Function.jpg"), width = 5, height = 5, units = "in", dpi = 300)




#--------------------Plot voronoi tessellation for validation ---------------#
dat <- PDACdata %>%
  filter(sample == '2117') %>%
  filter(ROIs == 'ROI06')


# construct bounding box
Rect <- data.frame(cbind(x = c(min(dat$X), min(dat$X), max(dat$X), max(dat$X)), x = c(min(dat$Y), max(dat$Y), max(dat$Y), min(dat$Y))))

p <- ggplot(dat,aes(X, Y)) +
  theme_void() +
  #theme(legend.position = 'NA') +
  geom_voronoi(aes(fill = class), color = '#878585', outline = Rect, size = .125) +

  # to validate g-function values
  scale_fill_manual(values = c('Tumor' = '#f8f8f8', 'Macrophages' = '#f8f8f8', 'Bcells' = '#82b4e5',
                               'CD4 Tcells' = '#f8f8f8', 'CD8 Tcells' = '#4875b0', 'Other Cells' = '#f8f8f8', 'CD45 Other' = '#f8f8f8'))  +
  
  coord_fixed(ratio = 1) 
p
ggsave(p, file=paste0("Figures/Gfunc_2117_ROI06.png"), width = 8, height = 6, units = "in", dpi = 300)





#----- Save dataset to GraphPad for comparisons at distance intervals ---------#

# prepare data
# select 
list.cur_1 <- get(eval(paste0('list', 4, '_1.all')))
list.cur_2 <- get(eval(paste0('list', 4, '_2.all')))

list.long <- list.cur_1 %>%
  filter(Sample %in% long.group) %>%
  #filter(r <= 60 & r >= 40) %>%
  cbind.data.frame('long') %>%
  select(2, 5) %>%
  `colnames<-`(c("rs", "class"))
write.csv(list.long, './exportedFiles/CD8_Bcells_long.csv')


list.short <- list.cur_1 %>%
  filter(Sample %in% short.group) %>%
  #filter(r <= 60 & r >= 40) %>%
  cbind.data.frame('short') %>%
  select(2, 5) %>%
  `colnames<-`(c("rs", "class"))
write.csv(list.short, './exportedFiles/CD8_Bcells_short.csv')

test.list.short <- list.cur_1 %>%
  filter(Sample %in% short.group) %>%
  group_by(Sample, ROI) %>%
  summarise(max = max(rs))
  


list.comb <- rbind.data.frame(list.short, list.long)

t.test(list.short$rs, list.long$rs)
# add comparisons
mycomparison <- compare_means(rs ~ class,  data = list.comb)
my_comparisons <- list( c("short", "long"))
p <- ggplot(list.comb, aes(class, rs, color = class)) +
  geom_boxplot() +
  theme_classic() +
  scale_color_manual(values = c('short' = '#c22828', 'long' = '#325698')) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22)) +
  geom_signif(comparisons = list(c("short", "long")), 
              map_signif_level=TRUE, size = 1, textsize = 10, color = 'black', 
              y_position = 1, annotations = '****'
              ) +
  xlab('') +
  ylab('G(r) value') +
  #ylab('') +
  ylim(0, 1.1)
p
ggsave(p, file=paste0("Figures/All_comparison_CD4_MPG.jpg"), width = 3, height = 5, units = "in", dpi = 300)


