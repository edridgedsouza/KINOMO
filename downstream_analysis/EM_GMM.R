#############################
#Script = EM_GMM.R
#Author = Somnath Tagore
#Last Update = 03.20.2019
##############################

#Packages

library(tidyverse)
library(reshape2)
library(matrixStats)
library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(purrr)
library(DropletUtils)
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(celldex)
library(mixtools)
library(gridExtra)
# load gene expression data/raw counts
gene_exp<-readRDS(file="gene_exp.rds")

# normalize the data

### use the seurat workflow
gene_exp <- NormalizeData(gene_exp)
gene_exp <- FindVariableFeatures(gene_exp)
gene_exp <- ScaleData(gene_exp)
gene_exp <- RunPCA(gene_exp)

# or use the tpm normalization
# for(i in 1:ncol(gene_exp)){
#   gene_exp[,i] <- 1E6*gene_exp[,i]/sum(gene_exp[,i])
# }

#Mixture model implementation
#Select one gene at a time
#gene_exp<-matrix(rnorm(36),nrow=1)
gene_exp<-gene_exp[gene_symbol,]

mygrobsem<-list()
gg.mixEM <- function(EM) {
#  require(ggplot2)
  x       <- with(EM,seq(min(x),max(x),len=1000))
  pars    <- with(EM,data.frame(comp=colnames(posterior), mu, sigma,lambda))
  em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
  em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))
  saveRDS(em.df,file="em.df.rds")
  ggplot(data.frame(x=EM$x),aes(x,y=..density..)) + 
    #type = rep(c('tumor', 'normal'), c(503,313)) +
    #geom_histogram(fill=NA,color="black")+
    geom_polygon(data=em.df,aes(x,y,fill=comp),color="grey50", alpha=0.5)+
    #geom_polygon(data=em.df,aes(x,y,fill=type),color="grey50", alpha=0.5)+
    scale_fill_discrete("Component\nMeans",labels=format(em.df$mu,digits=3))+
   # scale_fill_discrete("Component\nSD",labels=format(em.df$sd,digits=3))+
    theme_bw()
#  print(em.df$sd)
  #print(em.df)
   
}

set.seed(1)    # for reproducible example

# k	= Number of components. Initial value ignored unless mu and sigma are both NULL.
#gene_exp_gmm <- normalmixEM(gene_exp, k =2, lambda=NULL)
gene_exp_gmm <- normalmixEM(gene_exp, lambda=NULL, mu=NULL, sigma=NULL)
mygrobsem_gene_exp_gmm<-gg.mixEM(gene_exp_gmm)     

pdf(file = "gmm.pdf", width = 5, height = 5, family = "Times", pointsize = 10)
mygrobsem_gene_exp_gmm
dev.off()

em.df<- readRDS(file="em.df.rds")
# em.df %>%
#   group_by(comp) %>%
#   summarise_at(vars(y), list(name = mean))
aggregate(em.df, list(em.df$comp), mean)

#Select the component with hightest Mean to be the average expression of that gene

