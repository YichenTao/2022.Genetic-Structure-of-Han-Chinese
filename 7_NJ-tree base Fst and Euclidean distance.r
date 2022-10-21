library(rvcheck)
library(ggtree)
library(ape)

# NJ tree of based on Euclidean distance  ---------------------------------
# read in Frequency
Fst_Fre <- read.csv("./TABLE.S4_854_mt-HG_Frequency_of_33_provinces.csv")
Fst_Fre$Num <- NULL
Fst_Fre$Total <- NULL
Fst_Fre <- Fst_Fre %>% dplyr::filter(! Populations %in% c("Hongkong","Macau","Qinghai"))
row.names(Fst_Fre) <- Fst_Fre$Populations
Fst_Fre$Populations <- NULL
# genetation
Fst_dist <- Fst_Fre %>% dist()
tree <-nj(Fst_dist)
tree$edge.length <- tree$edge.length %>% abs
ggtree(tree) + 
  geom_tiplab(size=2)


# NJ tree of based on Fst  -----------------------------------------------
# read in Fst
Fst_Fre <- read.csv("./TABLE.S5_Pairwise_Fst_MT.csv")
Fst_Fre <- Fst_Fre %>% select(!c("Hongkong","Macau","Qinghai")) %>% 
    dplyr::filter(!X %in% c("Hongkong","Macau","Qinghai"))
row.names(Fst_Fre) <- Fst_Fre$X
Fst_Fre$X <- NULL
# genetation
Fst_dist <- Fst_Fre %>% as.dist()
tree <-nj(Fst_dist)
tree$edge.length <- tree$edge.length %>% abs
ggtree(tree) + 
    geom_tiplab(size=2)




