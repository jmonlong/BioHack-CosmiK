library(ggplot2)

dat = read.table("tcga.tar.gz", as.is=TRUE, skip=1)
dim(dat)
dat[1:5,1:5]
genes = dat[,1]
dat = dat[,-1]
summary(dat)
dat = dat[which(!is.na(dat[,1])),]

## PCA
pca.o = prcomp(as.matrix(dat))
plot(pca.o)
plot(pca.o$x)

## Median normalization
med.d = apply(dat,2,median)
qplot(med.d)
dat.mn = apply(dat, 2, function(cc)cc/median(cc)*mean(med.d))
pca.mn = prcomp(dat.mn)
plot(pca.o$x)

## Rescale
pca.mnr = prcomp(dat.mn, scale.=TRUE)
plot(pca.mnr$x)
library(scatterplot3d)
scatterplot3d(pca.mnr$x[,1:3])

## Correlation
cor.mat = cor(t(dat.mn))
