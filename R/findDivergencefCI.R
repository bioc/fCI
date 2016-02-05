


npci.gene.by.pvalues<-function(npci.data, gene.indexes,
	ctr.indexes, trt.indexes){
  pvalues=c()
  if(dim(npci.data)[2]>= length(unique(ctr.indexes, trt.indexes))){
    for(i in gene.indexes){
      pvalue=two.sample.log.ratio(as.numeric(npci.data[i, ctr.indexes]), 
		as.numeric(npci.data[i, trt.indexes]))
      pvalues=c(pvalues, pvalue)
    }
  }
  return(pvalues)
}


divergence.multivariate.distributions<-function(null.data, diff.data, choice=2){
  
  null.data=log(null.data,2)
  diff.data=log(diff.data,2)
  null.sample.size=dim(null.data)[1]
  diff.sample.size=dim(diff.data)[1]
  
  min.kl.dist=0
  kl.dist=0
  if(choice==1){	   
	kl.dist=crossentropy(null.data, diff.data, k=2, algorithm="kd_tree")  
     min.kl.dist=abs(kl.dist[1]) + sqrt(sum(kl.dist^2))
  }	
  else{	
    d=dim(null.data)[2] 
    mu1=apply(null.data, 2, mean)
    sigma1=var(null.data) * (null.sample.size-1) / (null.sample.size) 
    mu2=apply(diff.data, 2, mean)
    sigma2=var(diff.data) * (diff.sample.size-1) / (diff.sample.size)  
    
    if(choice==2){	# Kullback-Leibler distance
      kl.dist=0.5*(log(det(sigma2)/det(sigma1)) - d + 
		tr(solve(sigma2)%*%sigma1) +
		t(mu2-mu1)%*%(solve(sigma2))%*%(mu2-mu1))	*
        0.5*(log(det(sigma1)/det(sigma2)) - d +tr(solve(sigma1)%*%sigma2) + 
		t(mu1-mu2)%*%(solve(sigma1))%*%(mu1-mu2))
    }
    if(choice==3){	# Hellinger distance
      kl.dist=1-(((det(sigma1)^0.25)*(det(sigma2)^0.25))/
		(det(0.5*(sigma1+sigma2))^0.5))*(exp(-0.125*(t(mu1-mu2))%*%
			(solve(0.5*(sigma1+sigma2)))%*%(mu1-mu2)))
    }
    min.kl.dist=kl.dist	
  }
  return(min.kl.dist)
}	



fCI.call.by.index<-function(wt.indexes, df.indexes, data.file, 
                            use.normalization=FALSE, npci=NULL, short.report=TRUE){
  
  need.initialize=FALSE
  if(is.null(npci)){
    npci=new("NPCI")
    initialize("NPCI")
  }
  
  if(is.matrix(data.file) | is.data.frame(data.file)){
    npci@sample.data.normalized=data.file
  }else{
    npci@sample.data.file=data.file  
    sample.data.comb=read.csv(file=data.file, sep="\t", 
                         stringsAsFactors=FALSE) 
    none.zero.index=1:dim(sample.data.comb)[1]		
    none.zero.index=which(lapply(1:dim(sample.data.comb)[1], FUN=function(i){ 
      values=sample.data.comb[i,]; 
      any(is.na(values))==FALSE & min(values)>0 & max(values)<Inf })==TRUE)  
    sample.data.comb=sample.data.comb[none.zero.index,]
    npci@sample.data.normalized=sample.data.comb
    
  }  
  
  #print(npci@diff.gene.ids)
  
  npci@ctr.indexes=wt.indexes
  npci@trt.indexes=df.indexes
  
  wt.index.more=combinations(length(wt.indexes),2, 
                             v=wt.indexes, repeats.allowed=FALSE)
  wt.index.more=(lapply(1:dim(wt.index.more)[1], function(x) 
    as.numeric(wt.index.more[x,])))# list(a=wt.index,b=df.index))
  df.index.more=do.call(expand.grid, 
                        list(a=wt.indexes, b=df.indexes)) 
  df.index.more=(lapply(1:dim(df.index.more)[1], function(x) 
    as.numeric(df.index.more[x,])))
  
  pairwise.up.down=list()
  pairwise.index=list()
  pairwise.wt.up.down.fold=list()
  pairwise.df.up.down.fold=list()
  k=1
  
  if(use.normalization==TRUE){npci=normalization(npci)} ## just do once
  
  #print(npci@diff.gene.ids)
  
  for(i in 1:length(wt.index.more)){
    for(j in 1:length(df.index.more)){
      npci@wt.index=wt.index.more[[i]]
      npci@df.index=df.index.more[[j]]
      
      npci=populate(npci)
      npci=compute(npci)
      npci=summarize(npci)
      
      if(npci@result[2]<max(npci@fold.cutoff.list[[1]])){
        
        cat("Control-Control Used : [", npci@wt.index, "] & Control-Case Used : [", 
            npci@df.index , "];  Fold_Cutoff=", npci@result[2], "; Num_Of_DEGs=", 
            npci@result[1], "; Divergence=",npci@result[3], "\n")
        diff.gene.ids=npci@diff.gene.ids[[1]] 
        diff.gene.names=rownames(npci@sample.data.normalized)[diff.gene.ids]
        pairwise.up.down[k][[1]]=rep(-1, length(diff.gene.ids))
        pairwise.up.down[k][[1]][which(npci@expr.by.fold[1,diff.gene.ids]>1)]=1
        pairwise.index[k][[1]]= diff.gene.names
        pairwise.wt.up.down.fold[k][[1]]=as.double(npci@null.data.start) 
        pairwise.df.up.down.fold[k][[1]]=as.double(npci@diff.data.start)
        k=k+1
        
      }else{
        ## fci already incorporated gene ids with cutoff=10 into the results
        npci@diff.gene.ids=list()
      }
    }
  }
  #print(npci@diff.gene.ids)
  
  result=list(pairwise.index, pairwise.wt.up.down.fold, 
              pairwise.df.up.down.fold,  pairwise.up.down, npci)
  if(short.report==TRUE){
    result=report.target.summary(pairwise.index) 
    index.of.DEGs=unlist(lapply(as.vector(result[,1]), function(n){which(rownames(npci@sample.data.normalized)==n)}))
    mean.a=round(apply(cbind(npci@sample.data.normalized[index.of.DEGs, wt.indexes]), 1, mean),3)
    mean.b=round(apply(cbind(npci@sample.data.normalized[index.of.DEGs, df.indexes]), 1, mean),3)
    log.fc=round(log(mean.b/mean.a,2),3)
    result=data.frame(cbind(as.vector(result[,1]), mean.a, mean.b, log.fc, round(as.numeric(result[,2]),3)))
    colnames(result)=c("DEG_Names", "Mean_Control", "Mean_Case", "Log2_FC", "fCI_Prob_Score")
    if(dim(result)[1]>0){rownames(result)=1:dim(result)[1]}
  }
  return(result)
}


