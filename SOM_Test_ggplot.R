setwd("~/Desktop/SelfOorganizingMaps")

require(kohonen)

dat = read.table("tcga_RSEM_tranform.csv", header=TRUE, sep="\t", as.is=TRUE)

med.g = apply(dat,1,median)
hist(med.g, breaks=200)
max.g = apply(dat,1,max)
hist(max.g, breaks=200)
dat.0 = dat[max.g>0,]

#training set
dat_subset<-dat.0
dat_train<-as.matrix(scale(dat_subset))

som_grid <- somgrid(xdim = 20, ydim=20, topo="hexagonal")

som_model <- som(dat_train, grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
rm(som_grid, dat_train)

plot(som_model, type="changes")
plot(som_model, type="count")
plot(som_model, type="dist.neighbours")

str(som_model)

## Cluster
som_cluster = cutree(hclust(dist(som_model$codes)),6)

## Effect of expression filters
med.r = apply(dat.0, 1, median)
sd.r = apply(dat.0, 1, sd)
hist(sd.r)
hist(med.r)
dat_train = as.matrix(dat.0)[sd.r>.0 & med.r>0.6,1:30]
som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")
som_model <- som(dat_train, grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
plot(som_model, type="count")



## ggplot2
library(ggplot2)
library(dplyr)

gene.stats = data.frame(gene=rownames(dat_subset), med.exp=rowMeans(dat_subset), sd.exp=apply(dat_subset,1 , sd))

pdf = data.frame(gene.stats,
  somId = som_model$unit.classif,
  dist=som_model$distances)
unit.pdf = data.frame(somId=1:nrow(som_model$grid$pts), som_model$grid$pts,
  clust=cutree(hclust(dist(som_model$codes)),6))
pdf = merge(pdf, unit.pdf)
pdf.s = pdf %>% group_by(x, y, somId, clust) %>% summarize(nb.genes=n(), dist.med=median(dist), med.exp=median(med.exp), sd.exp=mean(sd.exp))

save(pdf, pdf.s, file="somDF-example.RData")

ggplot(pdf.s, aes(x=x,y=y)) + geom_point(size=7) + geom_point(aes(colour=nb.genes),size=6) + theme_bw() + scale_colour_gradient(low="red",high="white")
ggplot(pdf.s, aes(x=x,y=y)) + geom_point(size=7) + geom_point(aes(colour=med.exp),size=6) + theme_bw() + scale_colour_gradient(low="red",high="white")
ggplot(pdf.s, aes(x=x,y=y)) + geom_point(size=7) + geom_point(aes(colour=sd.exp),size=6) + theme_bw() + scale_colour_gradient(low="red",high="white")
ggplot(pdf.s, aes(x=x,y=y)) + geom_point(size=6) + geom_point(aes(colour=dist.med),size=5) + theme_bw() + scale_colour_gradient(low="red",high="white")
ggplot(pdf.s, aes(x=x,y=y)) + geom_point(size=6) + geom_point(aes(colour=factor(clust)),size=5) + theme_bw()

ggplot(pdf.s, aes(x=x,y=y)) + geom_point(aes(size=nb.genes)) + theme_bw()

somId.single = arrange(pdf.s, desc(nb.genes))$somId[1]
pdf$single.unit = pdf$somId == somId.single

ggplot(pdf, aes(x=med.exp, colour=single.unit)) + geom_density()

## Let's say we selected one somId
pdf.si = subset(pdf, somId = 33)
