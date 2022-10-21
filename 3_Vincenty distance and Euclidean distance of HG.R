library(dplyr)
library(ggplot2)

# read in location of each province --------------------------------------------------------------
data <- read_csv("./TABLE.S8_Latitude_and_Longitude_of_each_province.csv")
data <- data %>% arrange(Populations)
Vincenty_df <- matrix(0,nrow = nrow(data),ncol = nrow(data)) %>% as.data.frame()
row.names(Vincenty_df) <- data$Populations
colnames(Vincenty_df) <- data$Populations

# Vincenty distance, km --------------------------------------------------------------
library(geosphere)
for(i in 1:nrow(data)){
  for(j in 1:nrow(data)){
    Vincenty_df[i,j] <- distVincentyEllipsoid(data[i,2:3],data[j,2:3])/1000
  }
}

# Vincenty_dist_5Groups --------------------------------------------------------
info <- read.csv("./TABLE.S1_Classification_of_provinces.changed.csv")
info_class <- info %>% select(Populations,Classification.Y1.G5)
# Provinces with sample size less than 30 were removed
info_class <- info_class %>% dplyr::filter(!Populations %in% c("Hongkong","Macau","Qinghai"))
mygroup <- info_class$Classification.Y %>% unique() %>% sort()
geo_mean_df <- data.frame(mygroup,long=0,lat=0)
# Calculate the spatial coordinates of 5 groups
for(i in 1:length(mygroup)){
  select_prov <- info_class %>% 
    dplyr::filter(Classification.Y1.G5 == mygroup[i]) %>% .$Populations
  geo_matrix <- data %>% dplyr::filter(Populations %in% select_prov) %>% 
      select(longitude,latitude)
  geo_mean_df[i,2:3] <- geomean(geo_matrix)
}

# calculating
Vincenty_df <- matrix(0,nrow = nrow(geo_mean_df),ncol = nrow(geo_mean_df)) %>% as.data.frame()
row.names(Vincenty_df) <- geo_mean_df$mygroup
colnames(Vincenty_df) <- geo_mean_df$mygroup

for(i in 1:nrow(geo_mean_df)){
  for(j in 1:nrow(geo_mean_df)){
    Vincenty_df[i,j] <- distVincentyEllipsoid(geo_mean_df[i,2:3],geo_mean_df[j,2:3])/1000
  }
}


# Vincenty_dist_4Groups --------------------------------------------------------
info <- read.csv("./TABLE.S1_Classification_of_provinces.changed.csv")
info_class <- info %>% select(Populations,Classification.mtDNA.G4)
# Provinces with sample size less than 30 were removed
info_class <- info_class %>% dplyr::filter(!Populations %in% c("Hongkong","Macau","Qinghai"))
mygroup <- info_class$Classification.mtDNA.G4 %>% unique() %>% sort()
geo_mean_df <- data.frame(mygroup,long=0,lat=0)
# 计算几个组别空间坐标
for(i in 1:length(mygroup)){
  select_prov <- info_class %>% 
    dplyr::filter(Classification.mtDNA.G4 == mygroup[i]) %>% .$Populations
  geo_matrix <- data %>% dplyr::filter(Populations %in% select_prov) %>% 
    select(Longitude,Latitude)
  geo_mean_df[i,2:3] <- geomean(geo_matrix)
}

# Calculate the spatial coordinates of 4 groups
Vincenty_df <- matrix(0,nrow = nrow(geo_mean_df),ncol = nrow(geo_mean_df)) %>% as.data.frame()
row.names(Vincenty_df) <- geo_mean_df$mygroup
colnames(Vincenty_df) <- geo_mean_df$mygroup

for(i in 1:nrow(geo_mean_df)){
  for(j in 1:nrow(geo_mean_df)){
    Vincenty_df[i,j] <- distVincentyEllipsoid(geo_mean_df[i,2:3],geo_mean_df[j,2:3])/1000
  }
}







