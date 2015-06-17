
get.npci.data<-function(sample.data.normalized, wt.index, df.index){
  
  wt.comb=as.numeric(wt.index)
  df.comb=as.numeric(df.index)
  num.var=length(wt.comb)/2
  null.data.start=matrix(0, dim(sample.data.normalized)[1], num.var)
  diff.data.start=matrix(0, dim(sample.data.normalized)[1], num.var)
  for(k in 1:num.var){
    null.data.start[,k]=sample.data.normalized[,wt.comb[2*k]]/
		sample.data.normalized[,wt.comb[2*k-1]]
    diff.data.start[,k]=sample.data.normalized[,df.comb[2*k]]/
		sample.data.normalized[,df.comb[2*k-1]]
  }
  return(cbind(null.data.start, diff.data.start))
}


get.rank.combinations<-function(rank.index.to.be.removed, symmetric.fold){
  
  list.comb=list()
  if(symmetric.fold==FALSE){
    l.rank.index.to.be.removed=rank.index.to.be.removed[[1]]
    r.rank.index.to.be.removed=rank.index.to.be.removed[[2]]
    
    d=length(l.rank.index.to.be.removed)
    k=1
    for(i in 1:d){
      list.comb[[k]]=l.rank.index.to.be.removed[[i]] ; k=k+1
      list.comb[[k]]=r.rank.index.to.be.removed[[i]] ; k=k+1		
    }
  }else{
    list.comb=rank.index.to.be.removed
  }
  combinations=do.call(expand.grid, list.comb)
  return(combinations)
}


get.npci.distance.matrix<-function(npci.data, null.data.start, diff.data.start, 
	choice=2, rank.index.to.be.removed, expr.by.fold, ctr.indexes, trt.indexes, 
	use.intersect=FALSE, symmetric.fold=TRUE, fold.cutoff.list){	
  
  l=dim(null.data.start)[1]
  d=if(symmetric.fold==TRUE){length(rank.index.to.be.removed)}
	else{length(rank.index.to.be.removed[[1]])}
  combinations=get.rank.combinations(rank.index.to.be.removed, symmetric.fold)
  combinations.num=dim(combinations)[1]	
  distance.matrix=c()
  
  for(i in 1:combinations.num){
    index.to.be.removed=if(symmetric.fold==TRUE){
      lapply(1:d, FUN=function(j){max.rank=combinations[i,j];
        this.expr.by.fold=expr.by.fold[j,]
        this.expr.by.fold[which(this.expr.by.fold<1)]=
		1/(this.expr.by.fold[which(this.expr.by.fold<1)])
        order(this.expr.by.fold, decreasing=TRUE)[0:max.rank]})
    }else{
      lapply(1:d, FUN=function(j){
        l.max.rank=combinations[i,2*j-1]
        r.max.rank=combinations[i,2*j  ] 
        union(order(expr.by.fold[j,],	decreasing=TRUE )[0:r.max.rank],
              order(expr.by.fold[j,],	decreasing=FALSE)[0:l.max.rank])})
    }
    index.to.be.removed=if(use.intersect==TRUE){
		intersect.of.lists(index.to.be.removed)}
		else{unique(unlist(index.to.be.removed))}
    if(use.intersect==TRUE){num.removed=0; if(length(index.to.be.removed)>0)
		{num.removed=length(index.to.be.removed)};
		combinations[i,]=rep(num.removed, d)}
    case.control.remain.index=setdiff(1:l, index.to.be.removed)	
    null.data=null.data.start
    diff.data=cbind(diff.data.start[case.control.remain.index,])
    min.kl.dist=divergence.multivariate.distributions(null.data, 
		diff.data, choice)  
 
    distance.matrix=c(distance.matrix, min.kl.dist)
  }
  return(distance.matrix)
}	



npci.index.to.be.removed<-function(expr.by.fold, d, symmetric.fold, 
	max.rank, l.max.rank, r.max.rank){
  index.to.be.removed=if(symmetric.fold==TRUE){
    lapply(1:d, FUN=function(j){
      this.expr.by.fold=expr.by.fold[j,]
      this.expr.by.fold[which(this.expr.by.fold<1)]=
		1/(this.expr.by.fold[which(this.expr.by.fold<1)])
      order(this.expr.by.fold, decreasing=TRUE)[0:max.rank]})
  }else{
    lapply(1:d, FUN=function(j){
      union(order(expr.by.fold[j,],	decreasing=TRUE )[0:r.max.rank],
            order(expr.by.fold[j,],	decreasing=FALSE)[0:l.max.rank])})
  }
  return(index.to.be.removed[[1]])
}



npci.index.reconsidered<-function(npci.data, expr.by.fold, null.data.start, 
	diff.data.start, gene.indexes, ctr.indexes, 
	trt.indexes, left.fold, right.fold){
  
  left.size=round(length(which(null.data.start<(1/left.fold)))*
	(dim(null.data.start)[1]-length(gene.indexes))/(dim(null.data.start)[1]))
  right.size=round(length(which(null.data.start>right.fold))*
	(dim(null.data.start)[1]-length(gene.indexes))/(dim(null.data.start)[1]))
  pvalues=npci.gene.by.pvalues(npci.data, gene.indexes, 
	ctr.indexes, trt.indexes)
  
  cat(left.size, " = ", right.size, "\n")
  
  chosen.index=which(expr.by.fold[, gene.indexes]<1)
  big.pvalues.order=order(pvalues[chosen.index], 
	decreasing=TRUE)[1:(if(left.size>length(chosen.index))
		{length(chosen.index)}else{left.size})]
  down.regulation=gene.indexes[chosen.index][big.pvalues.order]
  
  chosen.index=which(expr.by.fold[, gene.indexes]>1)
  big.pvalues.order=order(pvalues[chosen.index], 
	decreasing=TRUE)[1:(if(right.size>length(chosen.index))
	{length(chosen.index)}else{right.size})]
  up.regulation=gene.indexes[chosen.index][big.pvalues.order]
  
  indexes.chosen=c()
  for(index in c(down.regulation, up.regulation)){
    if(log(diff.data.start[index,],2)>(mean(log(diff.data.start[-gene.indexes,],
		2)))-4*sd(log(diff.data.start[-gene.indexes,],2)) & 
         log(diff.data.start[index,],2)<(mean(log(diff.data.start[-gene.indexes,
		 ],2)))+4*sd(log(diff.data.start[-gene.indexes,],2))){
      indexes.chosen=c(indexes.chosen, index)
    }
  }
  return(indexes.chosen)
}

