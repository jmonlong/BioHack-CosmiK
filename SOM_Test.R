setwd("~/Desktop/SelfOorganizingMaps")

#Load libraries
require(kohonen)
library(dummies)
library(ggplot2)
library(sp)
library(maptools)
library(reshape2)
library(rgeos)

#Source functions
source('coolBlueHotRed.R')
source('plotHeatMap.R')

#nice plots colors
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')

#Read dataset
dat = read.table("tcga_RSEM_tranform 2.csv", header=TRUE, sep="\t", as.is=TRUE)

#Clean dataset
med.g = apply(dat,1,median)
hist(med.g, breaks=200)
max.g = apply(dat,1,max)
hist(max.g, breaks=200)
dat.0 = dat[max.g>0,]

#pick a training set (in this case all of the data)
dat_subset<-dat.0
dat_train<-as.matrix(scale(dat_subset))

#Set up grid size
som_grid <- somgrid(xdim = 20, ydim=20, topo="hexagonal")

#Train SOM
som_model <- som(dat_train, grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
rm(som_grid, dat_train)

#Data visualisation

plot(som_model, type = "counts", main="Node Counts", palette.name=coolBlueHotRed)
plot(som_model, type = "quality", main="Node Quality/Distance", palette.name=coolBlueHotRed)
plot(som_model, type = "counts", main="Node Counts", palette.name=coolBlueHotRed)
plot(som_model, type="dist.neighbours", main = "SOM neighbour distances", palette.name=grey.colors)


