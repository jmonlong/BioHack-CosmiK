library(shiny)
library(animint)
library(dplyr)
library(ggplot2)
library(kohonen)
library(reshape)

load("dat.RData")

gene.stats = data.frame(gene=rownames(dat), med.exp=rowMeans(dat), sd.exp=apply(dat,1 , sd), stringsAsFactors=FALSE)
    

kohUmatrix <- function(x){
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

  get.dat <- reactive({
    if(input$med.norm){
      med.global = median(dat, na.rm=TRUE)
      dat = t(apply(dat, 1, function(ee)ee/median(ee, na.rm=TRUE)*med.global))
    }
    return(dat)
  })
  
  som.model.comp <- reactive({
    som_grid <- somgrid(xdim = input$nb.units, ydim=input$nb.units, topo="hexagonal")
    dat = get.dat()
    som_model <- som(dat[, sample.int(ncol(dat), input$nb.samp)], grid=som_grid, rlen=input$nb.iter, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
    return(som_model)
  })
  
  som.comp <- reactive({
    som_model = som.model.comp()
    ## Summary stats
    pdf = data.frame(gene.stats,
      somId = som_model$unit.classif,
      dist=som_model$distances)
    unit.pdf = data.frame(somId=1:nrow(som_model$grid$pts), som_model$grid$pts,
      clust=cutree(hclust(dist(som_model$codes)),input$nb.clust))
    ##
    unit.pdf$dist.neighbours = kohUmatrix(som_model)
    ##
    pdf = merge(pdf, unit.pdf)
    pdf.s = pdf %>% group_by(x, y, somId, clust, dist.neighbours) %>% summarize(nb.genes=n(), dist.med=median(dist), med.exp=median(med.exp), sd.exp=mean(sd.exp))
    iter.df = data.frame(iter=1:length(som_model$changes), changes=som_model$changes)
    return(list(pdf=pdf, pdf.s=pdf.s, iter.df=iter.df))
  })

  output$anim = renderAnimint({
    som.o = som.comp()
    pdf = as.data.frame(som.o$pdf)
    pdf.s = as.data.frame(som.o$pdf.s)
    pdf.s$col = pdf.s[,input$col]
    if(input$col=="clust"){
      scCol = scale_colour_hue()
      pdf.s$col = factor(pdf.s$col)
    } else {
      scCol = scale_colour_gradient(name="", low="yellow",high="red")
    }
    
    som <- ggplot() + geom_point(data=pdf.s, aes(x=x, y=y), size=7, colour="black") + geom_point(data=pdf.s, aes(x=x, y=y, colour=col, clickSelects=somId),size=6)   + theme_bw()  + theme(line = element_blank(),rect=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) + xlab("") + ylab("") + scCol
    
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
    
    list(som=som, msd=med.sd, pca=pca.p, gene=gene.name)
  })

  output$som = renderPlot({
    som.o = som.comp()
    pdf.s = as.data.frame(som.o$pdf.s)
    pdf.s$col = pdf.s[,input$col]
    if(input$col=="clust"){
      scCol = scale_colour_hue(name="")
      pdf.s$col = factor(pdf.s$col)
    } else {
      scCol = scale_colour_gradient(name="", low="yellow",high="red")
    }
    
    ggplot(pdf.s, aes(x=x, y=y, colour=col)) +
      geom_point(size=8, colour="black")  + geom_point(size=7) +
        theme_bw()  + theme(line = element_blank(),rect=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="top") +
          xlab("") + ylab("") + scCol

  })

  output$som.conv = renderPlot({
    som.o = som.comp()
    
    ggplot(som.o$iter.df, aes(x=iter, y=changes)) +
      geom_point(size=2) + geom_line() + 
        theme_bw()  +
          xlab("Iteration") + ylab("Changes")
  })

  output$som.codes = renderPlot({
    som_model = som.model.comp()
    codes.df = melt(som_model$codes)
    colnames(codes.df) = c("somId","sample","expression")
    grid.df = data.frame(somId = 1:nrow(som_model$grid$pts), som_model$grid$pts)
    codes.df = merge(codes.df, grid.df)
    codes.df = codes.df %>% mutate(exp.n=expression/max(expression), samp.n=as.numeric(factor(sample))/length(unique(sample))) %>% mutate(x.n=x+samp.n, y.n=y+exp.n)
    
    ggplot(codes.df) + geom_rect(data=grid.df, aes(xmin=x, xmax=x+1, ymin=y, ymax=y+1), alpha=0, colour="black") + geom_point(aes(x=x.n, y=y.n)) + theme_bw() + theme(line = element_blank(),rect=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="top") + xlab("") + ylab("expression")
  })
  
})
