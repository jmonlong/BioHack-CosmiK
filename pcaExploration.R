library(ggplot2)
library(scatterplot3d)
#library(rgl)

#Loading data
data = read.table("/home/phc/Desktop/BioHack-CosmiK/tcga_RSEM_tranform.csv", header=TRUE, sep="\t", as.is=TRUE)

#Removing the genes that contains at least a 0
anyZero = apply(data,1,function(x){min(x)!=0})
data.0 = data[anyZero,]

#Puts the data.0 file into matrix format
data.0 = as.matrix(data.0)

#PCA
pca.mnr = prcomp(data.0, scale.=TRUE)

pca.x = pca.mnr$x[,1:3] #Save PCA values
pca.x$new.col = rep(0,nrow(data.0)); #Creates a new collom to add features value

#Plotting PCA
plot(pca.mnr$x)
scatterplot3d(pca.mnr$x[,1:3])

#Interactive 3d plot for better perspective on color
#plot3d(pca.mnr$x)

