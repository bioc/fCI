
deg.up.down.info<-function(wt.index.in.list, df.index.in.list, npci,
                           use.normalization=FALSE, target.ratio=0.5){
  
  wt.index.more=do.call(expand.grid, wt.index.in.list) 	
  df.index.more=do.call(expand.grid, df.index.in.list) 
  wt.index.more=(lapply(1:dim(wt.index.more)[1], function(x) 
    as.numeric(wt.index.more[x,])))
  df.index.more=(lapply(1:dim(df.index.more)[1], function(x) 
    as.numeric(df.index.more[x,])))
  
  pairwise.up.down=list()
  pairwise.index=list()
  pairwise.wt.up.down.fold=list()
  pairwise.df.up.down.fold=list()
  k=1
  
  if(use.normalization==TRUE){npci=normalization(npci)} ## just do once
  
  for(i in 1:length(wt.index.more)){
    for(j in 1:length(df.index.more)){
      npci@wt.index=wt.index.more[[i]]
      npci@df.index=df.index.more[[j]]
      
      npci=populate(npci)
      npci=compute(npci)
      cat("Control-Control Used : ", npci@wt.index, " & Control-Case Used : ", 
          npci@df.index , " ; Fold_Cutoff=", npci@result[2], "; Num_Of_DEGs=", 
          npci@result[1], "; Divergence=",npci@result[3], "\n")
      npci=summarize(npci)
      
      
      diff.gene.ids=npci@diff.gene.ids[[1]] 
      pairwise.up.down[k][[1]]=rep(-1, length(diff.gene.ids))
      pairwise.up.down[k][[1]][which(npci@expr.by.fold[1,diff.gene.ids]>1)]=1
      pairwise.index[k][[1]]= diff.gene.ids
      pairwise.wt.up.down.fold[k][[1]]=as.double(npci@null.data.start) 
      pairwise.df.up.down.fold[k][[1]]=as.double(npci@diff.data.start)
      k=k+1
    }
  }
  
  common.ids=Reduce(intersect, pairwise.index)  
  union.ids=Reduce(union, pairwise.index)  
  up.down.sum=rep(0, length(common.ids))
  
  return(list(common.ids, up.down.sum, pairwise.wt.up.down.fold, 
              pairwise.df.up.down.fold, union.ids, pairwise.index, pairwise.up.down))
}


pairwise.change.occupancy<-function(common.ids, pairwise.index, 
	pairwise.up.down, target.ratio){
  
  up.down.sum=c()
  for(gid in common.ids){
    deg.sample.ids=c()
    for(j in 1:length(pairwise.index)){
      fid=which(pairwise.index[[j]]==gid)
      if(length(fid)>0){
        deg.sample.ids=c(deg.sample.ids, j)
      }
    }
    up.down.label=0
    for(sid in deg.sample.ids){
      diff.gene.ids=pairwise.index[[sid]]
      up.down.label=
		up.down.label+pairwise.up.down[[sid]][which(diff.gene.ids==gid)]
    }
    up.down.sum=c(up.down.sum, up.down.label)
  }
  
  condition=which(abs(up.down.sum)>=round(target.ratio*length(pairwise.index)))
  common.ids=common.ids[condition]
  up.down.sum=up.down.sum[condition]
  
  return (list(common.ids, up.down.sum))
}


deg.pairwise.fold.change<-function(pairwise.wt.up.down.fold, 
	pairwise.df.up.down.fold, d=1, min.fold=1.2){
  gene.num=length(pairwise.wt.up.down.fold[[1]])/d
  good.ids=c()
  for(i in 1:gene.num){
    sample.ids=unlist(lapply(1:d, function(n) (i+(n-1)*gene.num)))
    wt.fold=unlist(lapply(1:length(pairwise.wt.up.down.fold), function(n) 
		pairwise.wt.up.down.fold[[n]][sample.ids]))
    df.fold=unlist(lapply(1:length(pairwise.df.up.down.fold), function(n) 
		pairwise.df.up.down.fold[[n]][sample.ids]))
    if(all(df.fold>min.fold)==TRUE | all(df.fold<(1/min.fold))==TRUE){
      good.ids=c(good.ids, i)
    }
  }
  return(good.ids)
}


