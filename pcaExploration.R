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
pca.mnr = prcomp(dat, scale.=TRUE)
plot(pca.mnr$x)
library(scatterplot3d)
scatterplot3d(pca.mnr$x[,1:3])
library(rgl)
plot3d(pca.mnr$x)

## Remove 0s
med.g = apply(dat,1,median)
hist(med.g, breaks=200)
max.g = apply(dat,1,max)
hist(max.g, breaks=200)
dat.0 = dat[med.g>0.8,]
med.0g = apply(dat.0,1,median)
hist(med.0g, breaks=200)
hist(as.numeric(unlist(dat.0)), breaks=200)

## NA for 0
pca.0 = prcomp(dat.0, scale.=TRUE)
plot(pca.mnr$x)


## Correlation
cor.mat = cor(t(dat.mn))
