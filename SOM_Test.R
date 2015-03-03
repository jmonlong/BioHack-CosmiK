setwd("~/Desktop/SelfOorganizingMaps")

require(kohonen)

dat = read.table("tcga_RSEM_tranform 2.csv", header=TRUE, sep="\t", as.is=TRUE)

med.g = apply(dat,1,median)
hist(med.g, breaks=200)
max.g = apply(dat,1,max)
hist(max.g, breaks=200)
dat.0 = dat[max.g>0,]

#training set
dat_subset<-dat.0
dat_train<-as.matrix(scale(dat_subset))

#dat_train<-as.matrix(scale(dat.0))
#dat_train<-as.matrix(as.numeric(unlist(dat.0)))
som_grid <- somgrid(xdim = 20, ydim=20, topo="hexagonal")

som_model <- som(dat_train, grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
rm(som_grid, dat_train)

#nice plots colors
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')

plot(som_model, type="changes")
plot(som_model, type="count")
plot(som_model, type="dist.neighbours")

plot(som_model, type = "counts", main="Node Counts", palette.name=coolBlueHotRed)
plot(som_model, type = "quality", main="Node Quality/Distance", palette.name=coolBlueHotRed)
plot(som_model, type = "counts", main="Node Counts", palette.name=coolBlueHotRed)

