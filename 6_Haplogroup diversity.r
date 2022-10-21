library(dplyr)
library(ggplot2)

# read in mt Fre  ---------------------------------------------------------------
data <- read.csv("./TABLE.S4_854_mt-HG_Frequency_of_33_provinces.csv")
data$Total <- NULL 

# Function of HD value -------------------------------------------------------------
COMPUTE_HVALUE <- function(data){
  data <- data[(data$Num >50),]
  df <- data.frame(Populations = data$Populations,
                   N = data$Num,
                   H_value=NA)
  for(i in 1:nrow(df)){
    df_sub <- data[i,-c(1:2)] 
    df$H_value[i] <- (1-sum((df_sub)^2)) *df$N[i]/(df$N[i]-1)
  }
  df %>% arrange(H_value) %>% return()
}

H_value_df <- COMPUTE_HVALUE(data)
H_value_df <- H_value_df %>% arrange(desc(H_value))








