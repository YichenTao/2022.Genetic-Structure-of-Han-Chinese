library(dplyr)
library(ggplot2)
library(vegan)

# Mantel test -------------------------------------------------------------
# read in files
Fst_Fre <- read.csv("./TABLE.S4_854_mt-HG_Frequency_of_33_provinces.csv")
Vincenty_df <- read.csv("./TABLE.S6_Vincenty distance_by_province.csv")
info <- read.csv("./TABLE.S1_Classification_of_provinces.changed.csv")

# Mantel test (Fst vs Vincenty)-----------------------------------------------------------
# Provinces with a sample size of less than 30 were excluded
REMOVE_HK.MC <- function(DF) {
  DF <- DF %>% arrange(X)
  DF[,c("Hongkong","Macau","Qinghai")] <- NULL
  DF <- DF %>% dplyr::filter(!X %in% c("Hongkong","Macau","Qinghai"))
  row.names(DF) <- DF$X
  DF$X <- NULL
  DF <- DF %>% as.dist()  
  return(DF)
}
Vincenty_df <- REMOVE_HK.MC(Vincenty_df)
Fst_Fre <- REMOVE_HK.MC(Fst_Fre)
# Mantel text!
mantel(Vincenty_df,Fst_Fre,permutations=10000)

# intra Northern/Southern China
N_S <- "Northern China"
selected_pro <- info %>% dplyr::filter(Classification.2 %in% N_S) %>% .$Populations
# Provinces with a sample size of less than 30 were excluded
REMOVE_HK.MC <- function(DF) {
  DF <- DF %>% arrange(X)
  DF <- DF %>% select(! c("Hongkong","Macau","Qinghai")) 
  DF <- DF %>% dplyr::filter(!X %in% c("Hongkong","Macau","Qinghai"))
  DF <- DF %>% select(X,selected_pro[!(selected_pro %in% intersect(selected_pro,c("Hongkong","Macau","Qinghai")))]) 
  DF <- DF %>% dplyr::filter(X %in% selected_pro[!(selected_pro %in% intersect(selected_pro,c("Hongkong","Macau","Qinghai")))])
  row.names(DF) <- DF$X
  DF$X <- NULL
  DF <- DF %>% as.dist()  
  return(DF)
}
Vincenty_df <- REMOVE_HK.MC(Vincenty_df)
Fst_Fre <- REMOVE_HK.MC(Fst_Fre)
nrow(Vincenty_df)
# Mantel text!
mantel(Vincenty_df,Fst_Fre,permutations=10000)


# Mantel test (PC vs Vincenty) ------------------------------------------
pca.result <- read.csv("./TABLE.S9_Euclidean distance_in_PCA.MT.csv") 
pca.new <- pca.result %>% select(Populations:PC3)
pca.new <- pca.new %>% arrange(Populations) 
row.names(pca.new) <- pca.new$Populations
pca.new$Populations <- NULL
pca.dist <- dist(pca.new)
mantel(Vincenty_df,pca.dist,permutations=10000)
plot(Vincenty_df,pca.dist)

# intra Northern/Southern China
pca.result <- read.csv("./TABLE.S9_Euclidean distance_in_PCA.MT.csv") 
# choose region
N_S <- "Southern China"
N_S <- "Northern China"
selected_pro <- info %>% dplyr::filter(Classification.2 %in% N_S) %>% .$Populations
# 
REMOVE_HK.MC <- function(DF) {
  DF <- DF %>% arrange(X)
  DF <- DF %>% select(! c("Hongkong","Macau","Qinghai")) 
  DF <- DF %>% dplyr::filter(!X %in% c("Hongkong","Macau","Qinghai"))
  DF <- DF %>% select(X,selected_pro[!(selected_pro %in% intersect(selected_pro,c("Hongkong","Macau","Qinghai")))]) 
  DF <- DF %>% dplyr::filter(X %in% selected_pro[!(selected_pro %in% intersect(selected_pro,c("Hongkong","Macau","Qinghai")))])
  row.names(DF) <- DF$X
  DF$X <- NULL
  DF <- DF %>% as.dist()  
  return(DF)
}
Vincenty_df <- REMOVE_HK.MC(Vincenty_df)
pca.new <- pca.result %>% select(Populations:PC3)
pca.new <- pca.new %>% arrange(Populations) 
pca.new <- pca.new %>% dplyr::filter(Populations %in% selected_pro)
row.names(pca.new) <- pca.new$Populations
pca.new$Populations <- NULL
pca.dist <- dist(pca.new)
mantel(Vincenty_df,pca.dist,permutations=10000)
plot(Vincenty_df,pca.dist)


# Vincent 4 group ------------------------------------------------------------------
Vincenty_df_4group <- read.csv("./TABLE.S7_Vincenty_distance_of_different_groups_4.xlsx")
Vincenty_df_4group$X <- NULL
Vincenty_df_4group <- Vincenty_df_4group %>% as.dist()

# PCA 4 group
pca.result <- read.csv("./TABLE.S9_Euclidean distance_in_PCA.MT.csv") 
pca.result <- pca.result %>% select(Populations:PC3,Classification.mtDNA.G4)
mygroup <- pca.result$Classification.mtDNA.G4 %>% unique() %>% sort()
pc.mean <- data.frame(mygroup,PC1=0,PC2=0,PC3=0)
# Calculate PC distances for several groups
for(i in 1:length(mygroup)){
  select_prov <- pca.result %>% 
    dplyr::filter(Classification.mtDNA.G4 == mygroup[i]) %>% .$Populations
  pc_matrix <- pca.result %>% dplyr::filter(Populations %in% select_prov) %>% 
    select(PC1:PC3)
  pc.mean[i,2:4] <- pc_matrix %>% colSums()
}
row.names(pc.mean) <- pc.mean$mygroup
pc.mean$mygroup <- NULL
pca.dist.4group <- dist(pc.mean)
mantel(pca.dist.4group,Vincenty_df_4group,permutations=10000)


# Fst 4
Fst_Fre <- read.csv("./Pairwise_Fst_MT.csv")
Fst_Fre <- Fst_Fre %>% dplyr::filter(!X %in% c("Hongkong","Macau","Qinghai")) %>% 
  select(!c("Hongkong","Macau","Qinghai"))
info <- info %>% select(Populations,Classification.Y)
mygroup <- info$Classification.mtDNA.G4 %>% unique() %>% sort()
names(Fst_Fre)[1] <- "Populations"
Fst_Fre <- left_join(Fst_Fre,info)
fst_df <- matrix(0,length(mygroup), length(mygroup)) %>% as.data.frame()
row.names(fst_df) <- colnames(fst_df) <- mygroup
for(i in 1:length(mygroup)){
  for(j in 1:length(mygroup)){
    pro_row <- Fst_Fre$Populations[Fst_Fre$Classification.mtDNA.G4 == mygroup[i]]
    pro_col <- Fst_Fre$Populations[Fst_Fre$Classification.mtDNA.G4 == mygroup[j]]
    fst_df[i,j] <- Fst_Fre %>% dplyr::filter(Populations %in% pro_row) %>%  
      select(pro_col) %>% as.matrix() %>% mean()
  }
}
fst_df_dist <- fst_df %>% as.dist()
mantel(Vincenty_df_4group,fst_df_dist,permutations=10000)













































































