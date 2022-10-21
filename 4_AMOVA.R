library(dplyr)
library(ggplot2)

# AMOVA -------------------------------------------------------------------
# input files
Fre <- read.csv("./TABLE.S4_854_mt-HG_Frequency_of_33_provinces.csv")
info_new <- read.csv("./TABLE.S1_Classification_of_provinces.changed.csv")
info_classification <- read.csv("./Data/MBE-18-mtDNA.各省人数.33Prov.csv")

## samples
sample <- info_classification %>% dcast(Populations~Haplogroup) 
row.names(sample) <- sample[,1]
sample <- sample[,-1] 
sample <- sample %>% t %>% as.data.frame() 
sample <- sample %>% select(!c("Hongkong","Macau","Qinghai")) 

## distances
Fre <- Fre %>% dplyr::filter(!Populations %in% c("Hongkong","Macau","Qinghai"))
distances <- t(Fre[,-c(1:3)])
distances <- distances %>% dist()

## structures
info_classification <- left_join(info_classification,info_new)
structures <- info_classification %>% 
  select(Populations,Classification.Xu.G6:Classification.10) %>% 
  unique() %>% arrange(Populations) 
## how many genetic structures there are
structures <- structures %>% dplyr::filter(!Populations %in% c("Hongkong","Macau","Qinghai"))
str_list <- vector("list",length = length(colnames(structures))-1)
for(i in 2:length(colnames(structures))){
  str_list[[i-1]] <- structures %>% select(colnames(structures)[i]) %>% unique()
}

# AMOVA output dataframe
AMOVA_DF <- matrix(0,ncol = 6,length(str_list)) %>% as.data.frame()
colnames(AMOVA_DF) <- c("Group name","Groupings","Number of Groups",
                        "Variations Between Group","Variations Between samples Within Group",
                        "Variations Within samples")
AMOVA_DF_pvalue <- AMOVA_DF
for(i in 1:length(str_list)){
  G_name <- str_list[[i]] %>% colnames()
  G_str <- structures %>% select(G_name) %>% unlist() %>% 
    as.factor() %>% data.frame(Group = .)
  G_type <- G_str %>% unlist() %>% unique() %>% as.character()
  # AMOVing
  amova_G <- amova(samples = sample, 
                   distances = sqrt(distances), 
                   structures = G_str)
  # hypothesis testing
  randtest_G <- randtest(amova_G, nrepet = 100)  
  # output dataframe
  AMOVA_DF_pvalue[i,1] <- AMOVA_DF[i,1] <- G_name
  AMOVA_DF_pvalue[i,2] <- AMOVA_DF[i,2] <- G_type %>% paste0(collapse = ",")
  AMOVA_DF_pvalue[i,3] <- AMOVA_DF[i,3] <- length(G_type)
  AMOVA_DF[i,4:6] <- amova_G$componentsofcovariance[1:3,2]
  AMOVA_DF_pvalue[i,4:6] <- randtest_G$pvalue %>% rev
}

# write.csv(AMOVA_DF,"./Data/AMOVA_value.csv")
# write.csv(AMOVA_DF_pvalue,"./Data/AMOVA_P-value.csv")

# one group amova --------------------------------------------------------------
amova_G <- amova(samples = sample, 
                 distances = sqrt(distances))
randtest_G <- randtest(amova_G, nrepet = 100)
randtest_G$pvalue


