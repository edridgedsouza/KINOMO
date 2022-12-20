#Integrate factors

library("dplyr")
library("plyr")
library("readr")
library("purrr")

data_all <- list.files(path = ".", pattern = "*.csv", full.names = TRUE) %>% lapply(read_csv) %>% reduce(full_join, by = "Genes")
data_all[is.na(data_all)] <- 0

#write.csv(data_all,file="tumor_gene_top10_correlation.csv")
