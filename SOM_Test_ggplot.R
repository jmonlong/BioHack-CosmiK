setwd("~/Desktop/BioHack-CosmiK")

library(kohonen)
library(class)
library(MASS)

require(kohonen)

dat = read.table("tcga_RSEM_tranform.csv", header=TRUE, sep="\t", as.is=TRUE)

#Removing the genes that contains at least a 0
anyZero = apply(dat,1,function(x){min(x)!=0})
dat.0 = dat[anyZero,]

#training set
dat_train<-as.matrix(scale(dat.0))

som_grid <- somgrid(xdim = 40, ydim=40, topo="hexagonal")

som_model <- som(dat_train, grid=som_grid, rlen=60, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
rm(som_grid, dat_train)

plot(som_model, type="changes")
plot(som_model, type="count")
plot(som_model, type="dist.neighbours", palette.name=grey.colors)

str(som_model)

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

save(pdf, pdf.s, file="CosmiK.SOM")

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

## Median expression
ggplot(pdf.si, aes(x=med.exp, y=sd.exp)) + geom_point(alpha=.3) + theme_bw() + xlab("median gene expression") + ylab("standard deviation")

ggplot(pdf, aes(x=med.exp, y=sd.exp, group=factor(somId))) + theme_bw() + xlab("median gene expression") + ylab("standard deviation") + stat_smooth(se=FALSE)

## Clustering


##
library(animint)

som <- ggplot() + geom_point(data=pdf.s, aes(x=x, y=y, colour=factor(clust), clickSelects=somId),size=6) + theme_bw() 

nb.bins = 100
pdf.bar.exp = pdf %>% group_by(somId) %>% mutate(me.bin=cut(med.exp,nb.bins)) %>% group_by(somId, me.bin) %>% summarize(count=n()) %>% group_by(me.bin) %>% mutate(me.x=mean(as.numeric(unlist(strsplit(gsub("\\(","",gsub("]","",as.character(me.bin))),",")))))

bar.exp <- ggplot() + geom_bar(data=pdf.bar.exp, aes(x=me.x, y=count, showSelected=somId), stat="identity", position="identity") + ylim(0,200)

nb.bins = 100
pdf.bar.sd = pdf %>% group_by(somId) %>% mutate(me.bin=cut(sd.exp,nb.bins)) %>% group_by(somId, me.bin) %>% summarize(count=n()) %>% group_by(me.bin) %>% mutate(me.x=mean(as.numeric(unlist(strsplit(gsub("\\(","",gsub("]","",as.character(me.bin))),",")))))

bar.sd <- ggplot() + geom_bar(data=pdf.bar.sd, aes(x=me.x, y=count, showSelected=somId), stat="identity", position="identity") + ylim(0,200) + xlab("standard deviation")


pdf$PC1 = rnorm(nrow(pdf))
pdf$PC2 = rnorm(nrow(pdf))

pca.p = ggplot() + geom_point(data=pdf, aes(x=PC1, y=PC2,showSelected=somId, clickSelects=gene))

gene.name = ggplot() + geom_text(data=pdf, aes(x=0,y=0,label=gene, showSelected=gene))

gg2animint(list(som=som, exp=bar.exp, sd=bar.sd, pca=pca.p, gene=gene.name), "som-bar")

ggplot() + 
  make_text(UStornadoCounts, -100, 50, "year", "Tornadoes in %d") +
  geom_polygon(aes(x=long, y=lat, group=group, clickSelects=state),
               data=USpolygons, fill="black", colour="grey") +
  geom_segment(aes(x=startLong, y=startLat, xend=endLong, yend=endLat,
                   showSelected=year),
               colour="#55B1F7", data=UStornadoes) + 
  theme(axis.line=element_blank(), axis.text=element_blank(), 
        axis.ticks=element_blank(), axis.title=element_blank())
ts <- ggplot() + 
  make_text(UStornadoes, 1980, 200, "state") +
  geom_bar(aes(year, count, clickSelects=year, showSelected=state),
           data=UStornadoCounts, stat="identity", position="identity") + 
  ylab("Number of Tornadoes") + 
  xlab("Year")

tornado.ts.bar <- list(map = map, ts = ts, width=list(map = 970, ts = 500),  height=list(500)) 
gg2animint(tornado.ts.bar, "animint-test-ts-bar")
