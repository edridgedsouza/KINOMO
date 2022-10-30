
#Step 1: Extracting factors and metagenes for each sample as predicted by KINOMO
#Run the following script for each sample separately. 
#Please rename the factor names using the notation: SampleName_Rank_Factor. 
#As an example: SampleName = ABC/ABC_1, Rank = R9, Factor = F1/F2/F3/F4/F5/F6/F7/F8/F9. Thus, it can be: ABC_R9_F1

#Step 1.1: Read the KINOMO rank file and metagene file
ABC_rank_9 <- readRDS(file="ABC_nmf_rank_9.rds")
ABC_top100 <- read.csv(file="ABC_KINOMO_nmf_rank_9_top100_W.csv")
dim(ABC_top100)

ABC_top100
rownames(ABC_top100) <- ABC_top100[,1]
ABC_top100 <- ABC_top100[,-1]
rownames(ABC_top100)

ABC_top100

ABC_genelist_top100 <- c(ABC_top100[,1],ABC_top100[,2],ABC_top100[,3],ABC_top100[,4],ABC_top100[,5],ABC_top100[,6],ABC_top100[,7],ABC_top100[,8],ABC_top100[,9])
ABC_genelist_top100 <- ABC_genelist_top100[!duplicated(ABC_genelist_top100)]
ABC_genelist_top100
saveRDS(ABC_genelist_top100,file="ABC_genelist_top100_rank9_unique.rds")

#Step 1.2: Read the matrix 'W' file
ABC_top100_rank9 <- read.csv(file="ABC_KINOMO_nmf_rank_9_W.csv")
rownames(ABC_top100_rank9) <- ABC_top100_rank9[,1]
ABC_top100_rank9 <- ABC_top100_rank9[,-1]
ABC_top100_rank9_unique <- ABC_top100_rank9[ABC_genelist_top100,]
ABC_top100_rank9_unique

#Step 1.3: Optional step of renaming the factors, if already done before
colnames(ABC_top100_rank9_unique)<-c('ABC_R9_F1','ABC_R9_F2','ABC_R9_F3','ABC_R9_F4','ABC_R9_F5','ABC_R9_F6','ABC_R9_F7','ABC_R9_F8','ABC_R9_F9')
ABC_top100_rank9_unique

saveRDS(ABC_top100_rank9_unique,file="ABC_top100_rank9_unique.rds")

#Step 2: Run Step 1 over all files recursively. 
#Note: Would be difficult if the ranks across different samples are different (which is ideally the case) and if the factors are not properly renamed

#Step 3: Integrate the factors (Assuming Step 2 is complete for all samples). 
best_rank_top100_genes<-c(rownames(ABC_top100_rank2_unique),rownames(ABC_1_top100_rank2_unique))

#Step 3.1: Remove duplicates
best_rank_top100_genes_unique<-best_rank_top100_genes[!duplicated(best_rank_top100_genes)]

out <- bind_rows(as.data.frame(ABC_top100_rank2_unique),as.data.frame(ABC_1_top100_rank2_unique))


#Step 3.2: Minor cleaning (optional)
library(dplyr)
out[is.na(out)] <- 0
best_rank_top100_genes_unique<-out[!duplicated(best_rank_top100_genes_unique),]

#Step 3.3: Save the object
saveRDS(best_rank_top100_genes_unique,file="best_rank_top100_genes_unique.rds")


#library(ComplexHeatmap)
#library(circlize)
#col_fun = colorRamp2(c(0, 20), c("white", "red"))
#col_fun(seq(0, 100))

#Step 3.4: Minor cleaning of the integrated matrix is necessary to check if duplicate genes and/or empty data is present or not
best_rank_top100_genes_unique_updated <- read.csv(file="best_rank_top100_genes_unique_updated.csv")
row.names(best_rank_top100_genes_unique_updated) <- best_rank_top100_genes_unique_updated[,1]
best_rank_top100_genes_unique_updated <-best_rank_top100_genes_unique_updated[,-1]
#rownames(Melanoma_top200_genes_unique)
best_rank_top100_genes_unique_updated[1:5,1:5]

#Step 3.5: Run Co-correlation
#Note: The distance could be calculate by different methods: Spearman, Pearson, Jaccard, Manhattan, etc
dim(mydata.cor)
# mydata.cor.1 <- mydata.cor
# mydata.cor.1[mydata.cor.1 < 0.4] <- 0
mydata.cor = cor(best_rank_top100_genes_unique_updated, method = c("spearman"))

library(RColorBrewer)
#palette = colorRampPalette(c("blue", "white", "red")) (100)
par(mar=c(0,0,0,0)+0.1)

palette = colorRampPalette(c("blue", "white", "red")) (100)

#Step 3.6: Save result
pdf("top100_genes_KINOMO_unique_best_rank_corelation_2.pdf",height=30,width=30)
par(mar=c(0,0,0,0)+0.1)

heatmap(x = mydata.cor, col = palette, symm = TRUE, labRow=rownames(mydata.cor), labCol=colnames(mydata.cor),
        margins=c(10,10))
dev.off()
saveRDS(mydata.cor,file="Correlation.rds")
