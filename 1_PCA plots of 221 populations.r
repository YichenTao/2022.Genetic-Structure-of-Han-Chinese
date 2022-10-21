library(dplyr)
library(ggplot2)

# Figure ------------------------------------------------------------------
# input data 
data <- read.csv("TABLE.S3_28_Y-HG_Frequency_of_221_Populations.csv",
                 header = T,stringsAsFactors = F)
info <- read.csv("./TABLE.S1_Classification_of_provinces.changed.csv")
data_het <- data[,-c(1:2)]
data_het_scale <- data_het %>% apply(2,scale) %>% as.data.frame()

# PCA
pc.cr <- prcomp(data_het_scale,retx = TRUE)
pc <- round(summary(pc.cr)$importance[2,],2) 
pca.result <- as.data.frame(pc.cr$x[,1:2])
pca.result <- cbind(pca.result,data[,1:2])
pca.result$PC2 <- -pca.result$PC2

# filter
pca.mydata <- pca.result %>% dplyr::filter(group == "Han")
pca.pubdata <- pca.result %>% dplyr::filter(group != "Han")
myshape <- c(9,15,17:19,0,1,2,3,7,8,14,11,4,5,6,10,12)
factor_df <- data.frame(group = unique(pca.pubdata$group), levels = myshape[1:length(unique(pca.pubdata$group))])
pca.pubdata <- left_join(pca.pubdata,factor_df)
pca.pubdata <- pca.pubdata %>% arrange(group)
pca.mydata <- left_join(pca.mydata,info)
pca.mydata$group <- pca.mydata$Classification.3
pca.mydata <- pca.mydata %>% select(PC1:group)
pca.mydata$group <- factor(pca.mydata$group,
                           levels = unique(pca.mydata$group)[c(2,4,1,3,0)],
                           labels = unique(pca.mydata$group)[c(2,4,1,3,0)])
# Figure 1: visualization
pn <- ggplot() + 
  geom_point(data = pca.pubdata, aes(PC1,PC2,shape=group),
             size =3.2,stroke=1.7,color ="grey") + 
  xlab(paste0("PC1  ",pc[1])) +
  ylab(paste0("PC2  ",pc[2])) +
  theme_bw() +
  scale_shape_manual(values = unique(pca.pubdata$levels)) +
  scale_color_nejm() +
  geom_point(data = pca.mydata, 
             aes(PC1,PC2,color=group),
             size =2.4,stroke=1.7,alpha=.7) +
  guides(color = guide_legend(title = "Geographic Region")) +
  guides(shape = guide_legend(title = "Populations"))  +
  theme(legend.title = element_text(size =16),
        legend.text = element_text(size =14),
        title = element_text(size =16),
        axis.text = element_text(size =14))
pn



