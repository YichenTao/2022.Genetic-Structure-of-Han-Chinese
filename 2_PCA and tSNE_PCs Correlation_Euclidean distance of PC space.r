library(dplyr)
library(plotly)
library(jsonlite)
library(Rtsne)
library(ggplot2)

# 2.3.1/figure 2 3D-PDA -------------------------------------------------------------------
# set colors
mycolor <- c( "#D62728FF", "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#9467BDFF")
# input mtDNA HG dataset
mt.data <- read_csv("./TABLE.S4_854_mt-HG_Frequency_of_33_provinces.csv")
mt.data <- mt.data %>% as.data.frame()
row.names(mt.data) <- mt.data$Populations
mt.data <- mt.data[,-c(1:3)]
# Provinces with a sample size of less than 30 were excluded
provience_rm <- c("Macau","Hongkong","Qinghai")
df_tmp <- mt.data[!(row.names(mt.data) %in% provience_rm),]
df <- df_tmp %>% apply(1,scale) %>% t 
# PCA
pc.cr <- prcomp(df) 
pc <- round(summary(pc.cr)$importance[2,],2) 
pca.result <- as.data.frame(pc.cr$x[,1:3])
pca.result$Populations <- row.names(pca.result)
info <- read.csv("./TABLE.S1_Classification_of_provinces.changed.csv")
pca.result <- left_join(pca.result,info)
pca.result$Classification.Y1.G5 <- factor(pca.result$Classification.Y1.G5,
  levels = c("Eastern Coast","Fujian and Taiwan",
  "Northern China","Lingnan","Upper and Middle Yangtze River"),
   labels = c("Eastern Coast","Fujian and Taiwan",
  "Northern China","Lingnan","Upper and Middle Yangtze River")) 
# visualization
transform(pca.result, am = Classification.Y1.G5) %>%
  plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, 
          color = ~am, colors = mycolor[1:5]) %>%
  add_markers() %>%
  add_text(text = ~Populations) %>% 
  layout(scene = list(
    xaxis = list(title = paste0("PC1 ",pc[1])),
    yaxis = list(title = paste0("PC2 ",pc[2])),
    zaxis = list(title = paste0("PC3 ",pc[3]))
  ))

# Figure 2 tSNE --------------------------------------------------------------------
# Load the corresponding R package into memory
tsne <- Rtsne(df, check_duplicates = FALSE, pca = F,
              perplexity=4, theta=0.0, dims=3,max_iter =2000)
tsne_retult <- as.data.frame(tsne$Y)
tsne_retult$Populations <- YOUR_FRE_DF %>% row.names() 
tsne_retult <- left_join(tsne_retult,info)
transform(tsne_retult, am = Classification.Y1.G5) %>%
    plot_ly(x = ~V1, y = ~V2, z = ~V3, 
            color = ~am, colors = mycolor) %>%
    add_markers() %>%
    add_text(text = ~Populations) %>% 
    layout(scene = list(
        xaxis = list(title = paste0("Dimension 1")),
        yaxis = list(title = paste0("Dimension 2")),
        zaxis = list(title = paste0("Dimension 3"))
    ))

# 2.3.2 Euclidean distance in PCA space------------------------------------------------------------------
pca.result0 <- pca.result %>% select(PC1:Populations)
row.names(pca.result0) <- pca.result0$Populations
pca.result0$Populations <- NULL
dist_mt <- pca.result0 %>% dist() %>% as.matrix()

# Figure 3 ----------------------------------------------------------------
# input geographical coordinates
longlat <- read_csv("./TABLE.S8_Latitude_and_Longitude_of_each_province.csv")
ggdata <- merge(longlat,pca.result)
ggdata <- left_join(ggdata,info)
# visualization function
GGPOINT <- function(X,Y,CLASS){  
  GGDATA <- ggdata 
  P <- cor.test(GGDATA[,X],GGDATA[,Y])$p.value %>% round(4)
  EST <- cor.test(GGDATA[,X],GGDATA[,Y])$estimate %>% round(2)
  ggdata %>% ggplot(aes(x=get(X),y=get(Y))) +
    geom_point(aes(color=get(CLASS)),size=3) +
    geom_smooth(method = 'lm',color='#c93756', linetype='dashed', size=1,level = 0) +
    theme_bw() +
    ggrepel::geom_text_repel(max.overlaps = 40,size =3,aes(label = Populations)) +
    geom_text(aes(x=min(get(X))+(max(get(X))-min(get(X)))/10,y=max(get(Y))-2),
              color = "grey40",size =4,label = paste0("p-value: ",P,"\nPCCs: ",EST)) +
    scale_color_nejm() +
    theme(legend.title = element_text(size =20),
          legend.text = element_text(size =18))
}

# PC1 and longitude
p1 <- GGPOINT("PC1","Longitude","Classification.3") + 
  ggtitle("PC1 with Longitude") + 
  labs(x="PC1",y="Longtitude") +
  guides(color = guide_legend(title = "Geographic Region"))
# PC1 and latitude
p2 <- GGPOINT("PC1","Latitude","Classification.3") + 
  ggtitle("PC1 with Latitude") + 
  labs(x="PC1",y="Latitude") +
  guides(color = guide_legend(title = "Geographic Region"))
# PC2 and longitude
p3 <- GGPOINT("PC2","Longitude","Classification.3") + 
  ggtitle("PC2 with Longitude") + 
  labs(x="PC2",y="Longtitude") +
  guides(color = guide_legend(title = "Geographic Region"))
# PC2 and latitude
p4 <- GGPOINT("PC2","Latitude","Classification.3") + 
  ggtitle("PC2 with Latitude") + 
  labs(x="PC2",y="Latitude") +
  guides(color = guide_legend(title = "Geographic Region"))






