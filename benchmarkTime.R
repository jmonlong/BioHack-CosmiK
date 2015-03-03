library(kohonen)
library(ggplot2)

dat = read.table("tcga_RSEM_tranform.csv", header=TRUE, sep="\t", as.is=TRUE)
max.g = apply(dat,1,max)
dat.0 = dat[max.g>0,]

## Time benchmark
library(parallel)
nbG = 5e3
grid.size=10

nbs.int = seq(10,100,10)
st.nbs = mclapply(nbs.int, function(nbs){
  st.o = system.time({som_grid <- somgrid(xdim = grid.size, ydim=grid.size, topo="hexagonal")
                      som_model <- som(as.matrix(scale(dat.0[1:nbG,1:nbs])), grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
                    })
  st.o[3]
}, mc.cores=4)
nbs.df = data.frame(nbs=nbs.int, time=unlist(st.nbs))
ggplot(nbs.df, aes(x=nbs, y=time)) + geom_point()

nbG = 5e3
grid.size.int=seq(5,20,2)
nbs = 50
st.grid = mclapply(grid.size.int, function(grid.size){
  st.o = system.time({som_grid <- somgrid(xdim = grid.size, ydim=grid.size, topo="hexagonal")
                      som_model <- som(as.matrix(scale(dat.0[1:nbG,1:nbs])), grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
                    })
  st.o[3]
}, mc.cores=4)
gs.df = data.frame(grid.size=grid.size.int, time=unlist(st.grid))
ggplot(gs.df, aes(x=grid.size, y=time)) + geom_point()

nbG.int = c(1:10)*1e3
grid.size = 10
nbs = 50
st.nbg = mclapply(nbG.int, function(nbG){
  st.o = system.time({som_grid <- somgrid(xdim = grid.size, ydim=grid.size, topo="hexagonal")
                      som_model <- som(as.matrix(scale(dat.0[1:nbG,1:nbs])), grid=som_grid, rlen=100, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
                    })
  st.o[3]
}, mc.cores=4)
nbg.df = data.frame(nbg=nbG.int, time=unlist(st.nbg))
ggplot(nbg.df, aes(x=nbg, y=time)) + geom_point()

nbG = 5e3
grid.size = 10
nbs = 50
nb.iter = seq(10,100,10)
st.iter = mclapply(nb.iter, function(rlen){
  st.o = system.time({som_grid <- somgrid(xdim = grid.size, ydim=grid.size, topo="hexagonal")
                      som_model <- som(as.matrix(scale(dat.0[1:nbG,1:nbs])), grid=som_grid, rlen=rlen, alpha=c(0.05,0.01), keep.data = TRUE, n.hood="circular" )
                    })
  st.o[3]
}, mc.cores=4)
iter.df = data.frame(nb.iter=nb.iter, time=unlist(st.iter))
ggplot(iter.df, aes(x=nb.iter, y=time)) + geom_point()
