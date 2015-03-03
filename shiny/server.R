library(shiny)
library(animint)
library(dplyr)
library(ggplot2)
library(kohonen)

load("dat.RData")

kohUmatrix <- function(x) 
{
    if (x$method != "som" & x$method != "supersom") 
        stop("Neighbour distance plot only implemented for (super)som")
    nhbrdist <- unit.distances(x$grid, x$toroidal)
    nhbrdist[nhbrdist > 1.05] <- NA
    if (x$method == "som") {
        for (i in 2:nrow(nhbrdist)) {
            for (j in 1:(i - 1)) {
                if (!is.na(nhbrdist[i, j])) 
                  nhbrdist[i, j] <- nhbrdist[j, i] <- dist(x$codes[c(i, 
                    j), ])
            }
        }
    }
    else {
        if (x$method == "supersom") {
            nhbrdist[!is.na(nhbrdist)] <- 0
            for (k in 1:length(x$data)) {
                for (i in 2:nrow(nhbrdist)) {
                  for (j in 1:(i - 1)) {
                    if (!is.na(nhbrdist[i, j])) 
                      nhbrdist[i, j] <- nhbrdist[i, j] + x$weights[k] * 
                        dist(x$codes[[k]][c(i, j), ])
                  }
                }
                nhbrdist[j, i] <- nhbrdist[i, j]
            }
        }
    }
    colSums(nhbrdist, na.rm = TRUE)
}


shinyServer(function(input, output) {
  som.comp <- shiny::reactive({
    som_grid <- somgrid(xdim = input$nb.units, ydim=input$nb.units, topo="hexagonal")
    som_model <- som(dat, grid=som_grid, rlen=60, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
    som_cluster = cutree(hclust(dist(som_model$codes)),6)
    ## Summary stats
    gene.stats = data.frame(gene=rownames(dat), med.exp=rowMeans(dat), sd.exp=apply(dat,1 , sd), stringsAsFactors=FALSE)
    pdf = data.frame(gene.stats,
      somId = som_model$unit.classif,
      dist=som_model$distances)
    unit.pdf = data.frame(somId=1:nrow(som_model$grid$pts), som_model$grid$pts,
      clust=cutree(hclust(dist(som_model$codes)),6))
    ##
    unit.pdf$dist.neighbours = kohUmatrix(som_model)
    ##
    pdf = merge(pdf, unit.pdf)
    pdf.s = pdf %>% group_by(x, y, somId, clust, dist.neighbours) %>% summarize(nb.genes=n(), dist.med=median(dist), med.exp=median(med.exp), sd.exp=mean(sd.exp))
    return(list(pdf=pdf, pdf.s=pdf.s))
  })

  output$compile = shiny::renderText({
    som.o = som.comp()
    pdf = som.o$pdf
    pdf.s = som.o$pdf.s
    
    som <- ggplot() + geom_point(data=pdf.s, aes(x=x, y=y, colour=dist.neighbours, clickSelects=somId, clickSelects2=clust),size=6) + theme_bw()  + theme(line = element_blank(),rect=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) + xlab("") + ylab("") + scale_colour_gradient(low="grey90",high="grey30")
    med.sd = ggplot() + geom_point(data=pdf,alpha=.3, aes(x=med.exp, y=sd.exp, showSelected=somId))  + xlab("median expression") + ylab("standard deviation")
    ##PCA
    pca.mnr = prcomp(as.matrix(dat), scale.=TRUE)
    pca.x = pca.mnr$x[,1:3] #Save PCA values
pca.x<-data.frame(pca.x)
    pca.x$gene<-row.names(pca.x)
    pdf$gene<-as.character(pdf$gene)
    pdf<-merge(pdf,pca.x,by="gene",all.x=T)
    pca.p = ggplot() + geom_point(data=pdf,colour="red", aes(x=PC1, y=PC2,showSelected=somId, clickSelects=gene))
    gene.name = ggplot() + geom_text(data=pdf, aes(x=0,y=0,label=gene, showSelected=gene)) + theme(line = element_blank(),rect=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) + xlab("") + ylab("")
    
    print(animint2dir(list(som=som, msd=med.sd, pca=pca.p, gene=gene.name), ".", open.browser=FALSE))
    return("done")
  })
  output$SOM = shiny::renderPlot({
    som.o = som.comp()
    pdf = som.o$pdf
    pdf.s = som.o$pdf.s
    pdf.s$col = as.data.frame(pdf.s)[,input$col]
    return(ggplot(pdf.s, aes(x=x, y=y, colour=col)) + geom_point(size=8, colour="black")  + geom_point(size=7) + theme_bw()  + theme(line = element_blank(),rect=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="top") + xlab("") + ylab("") + scale_colour_gradient(low="blue",high="red"))
  })
})
