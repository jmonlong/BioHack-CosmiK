library(ggplot2)

dat = read.table("tcga_RSEM_tranform.csv", header=TRUE, sep="\t", as.is=TRUE)
dim(dat)
dat[1:5,1:5]
summary(dat)

## PCA
pca.o = prcomp(as.matrix(dat))
plot(pca.o)
plot(pca.o$x)

## Rescale
pca.mnr = prcomp(dat.mn, scale.=TRUE)
plot(pca.mnr$x)
library(scatterplot3d)
scatterplot3d(pca.mnr$x[,1:3])

## Correlation
cor.mat = cor(t(dat.mn))
