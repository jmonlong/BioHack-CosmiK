#Remove all genes with zeros in the dataset

anyZ = apply(dat,1,function(x){min(x)!=0})

new_data<-dat[anyZ,]