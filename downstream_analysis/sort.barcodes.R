#sort the barcodes per factor from the H matrix

df1 <- read.csv(file="~/Downloads/filename_nmf_rank_4_H.csv")
row.names(df1) <- df1[,1]
df1<-df1[,-1]
rownames(df1)
df1 <- t(df1)
dim(df1)
df1[1:5,]
j1 <- max.col(df1, "first")
write.csv(j1,file="j1.csv")
dim(j1)
df2 <- data.frame(matrix(ncol = 0, nrow = length(rownames(df1))))
rownames(df2) <- rownames(df1)
df2$Group <- j1
dim(df2)
df2[1:10,]
write.csv(df2,file="df2.csv")
table(df2$Group)
