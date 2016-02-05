### =========================================================================
### build the fCI object
### -------------------------------------------------------------------------


setMethod("initialize", "NPCI", function(.Object){
  
  if(length(.Object@sample.data.file)>0 & dim(.Object@sample.data.normalized)[1]
	==0){
    sample.data=read.csv(file=.Object@sample.data.file, sep="\t", 
		stringsAsFactors=FALSE) 
    if(dim(sample.data)[1]<1){
    }
    else{
      .Object@sample.data.normalized=sample.data
    }
    sample.data.comb=.Object@sample.data.normalized
    none.zero.index=1:dim(sample.data.comb)[1]		
    none.zero.index=which(lapply(1:dim(sample.data.comb)[1], FUN=function(i){ 
		values=sample.data.comb[i,]; 
		any(is.na(values))==FALSE & min(values)>0 & max(values)<Inf })==TRUE)  
    sample.data.comb=sample.data.comb[none.zero.index,]
    if(dim(sample.data.comb)[1]<1){
      print("After filtering input data with missing values or zero values, 
		there is no data left! Please provide another valid file!")
    }
    if(dim(.Object@attr.info)[1]>0 & dim(.Object@sample.data.normalized)[1]
		==dim(.Object@attr.info)[1]){
      .Object@attr.info=.Object@attr.info[none.zero.index,]
    }
    .Object@sample.data.normalized=sample.data.comb
    
  }
  if(length(.Object@sample.data.file)==0 
	& dim(.Object@sample.data.normalized)[1]==0){
    #print("You didn't provide either the input file name or set the sample data 
		#from the object!")
  }
  
  .Object@method.option=2
  .Object@percent.genes.to.scan=0.9 
  .Object@num.genes.to.skip.each=5
  .Object@use.ratio=TRUE
  .Object@use.fold.change=TRUE
  .Object@symmetric.fold=TRUE
  .Object@center.by.gaussian.kernel=FALSE
  
  return(.Object)
})				


setGeneric("normalization",
           function(.Object)
             standardGeneric("normalization"))	
setMethod("normalization", "NPCI", 
          function(.Object){
            sample.data=.Object@sample.data.normalized
            sample.data.normalized=trim.size.normalization(sample.data) 
            .Object@sample.data.normalized=sample.data.normalized
            return(.Object)
          })


setGeneric("setfCI",
           function(.Object, wt.index, df.index, fold.cutoff.list, center.distribution)
             standardGeneric("setfCI"))  
setMethod("setfCI", "NPCI", 
          function(.Object, wt.index, df.index, fold.cutoff.list, center.distribution){
            .Object@wt.index=wt.index
            .Object@df.index=df.index
            .Object@fold.cutoff.list=list(fold.cutoff.list)
            .Object@center.by.gaussian.kernel=center.distribution
            return(.Object)
          })				





total.library.size.normalization<-function(sample.data){
  sample.1.sum=sum(sample.data[,1])
  sample.data.normalized=sample.data
  for(i in 2:dim(sample.data)[2]){
    sample.data.normalized[,i]=
		as.numeric(sample.data[,i])*sample.1.sum/(sum(sample.data[,i]))
  }
  return(sample.data.normalized)
}

trim.size.normalization<-function(sample.data){
  trim.cutoff=quantile(sample.data[,1], c(0.05,0.95))
  cutoff.index=which(sample.data[,1]>trim.cutoff[[1]] 
	& sample.data[,1]<trim.cutoff[[2]])
  sample.1.sum=sum(sample.data[cutoff.index,1])
  sample.data.normalized=sample.data
  for(i in 2:dim(sample.data)[2]){
    sample.data.normalized[,i]=
	 as.numeric(sample.data[,i])*sample.1.sum/(sum(sample.data[cutoff.index,i]))
  }
  return(sample.data.normalized)
}

deseq.median.ratio.normalization<-function(npci.data){
  npci.scale=matrix(0, dim(npci.data)[1], dim(npci.data)[2])
  for(i in 1:dim(npci.data)[1]){
    for(j in 1:dim(npci.data)[2]){
      npci.scale[i,j]=npci.data[i,j]/mean(as.numeric(npci.data[i,]))
    }
  }
  median.scales=unlist(lapply(1:dim(npci.data)[2], 
	FUN=function(i){median(npci.scale[,i])}))
  npci.rescaled.data=sweep(npci.data,MARGIN=2,median.scales,`/`)
}


setGeneric("populate",
           function(.Object)
             standardGeneric("populate"))	
setMethod("populate", "NPCI", 
          function(.Object){
            
            use.ratio=.Object@use.ratio
            percent.genes.to.scan=.Object@percent.genes.to.scan
            num.genes.to.skip.each=.Object@num.genes.to.skip.each
            use.fold.change=.Object@use.fold.change
            
            wt.comb=as.numeric(.Object@wt.index)
            df.comb=as.numeric(.Object@df.index)
            num.var=if(use.ratio==TRUE){length(wt.comb)/2}else{length(wt.comb)}
            sample.data.normalized=.Object@sample.data.normalized
            null.data.start=matrix(0, dim(sample.data.normalized)[1], num.var)
            diff.data.start=matrix(0, dim(sample.data.normalized)[1], num.var)
            
            if(use.ratio==TRUE){
              for(k in 1:num.var){
                null.data.start[,k]=sample.data.normalized[,wt.comb[2*k]] /
					sample.data.normalized[,wt.comb[2*k-1]]
                diff.data.start[,k]=sample.data.normalized[,df.comb[2*k]] /
					sample.data.normalized[,df.comb[2*k-1]]
              }
            }else{
              null.data.start=as.matrix(sample.data.normalized[, wt.comb])
              diff.data.start=as.matrix(sample.data.normalized[, df.comb])									
            }
            
             
            if(.hasSlot(.Object, "center.by.gaussian.kernel")==TRUE &
				length(slot(.Object, "center.by.gaussian.kernel"))==0){
              .Object@center.by.gaussian.kernel=FALSE
            }
            if(.hasSlot(.Object, "center.by.gaussian.kernel")==TRUE  &
				slot(.Object, "center.by.gaussian.kernel")==TRUE){
              null.peak.val=unlist(lapply(1:dim(null.data.start)[2], function(x) 
				{tdensity=density(log(null.data.start[,x],2), bw=0.5); 
					tdensity$x[which(tdensity$y==max(tdensity$y))] })) 
              diff.peak.val=unlist(lapply(1:dim(diff.data.start)[2], function(x) 
				{tdensity=density(log(diff.data.start[,x],2), bw=0.5); 
					tdensity$x[which(tdensity$y==max(tdensity$y))] }))
              for(i in 1:dim(null.data.start)[2]){null.data.start[,i]=
				null.data.start[,i]*2^(-null.peak.val[i])}
              for(i in 1:dim(diff.data.start)[2]){diff.data.start[,i]=
				diff.data.start[,i]*2^(-diff.peak.val[i])}  
                            
            }
            
            .Object@null.data.start=null.data.start
            .Object@diff.data.start=diff.data.start
            
            
            d=num.var
            l=dim(diff.data.start)[1]
            exclusion.max.num=round(l*percent.genes.to.scan)
            expr.by.fold=matrix(0, d, l)
            
            for(i in 1:d){
              expr.by.fold[i,]=as.numeric(lapply(1:l, FUN=function(j){fold=
				if(use.ratio==FALSE){diff.data.start[j,i]}
					else{diff.data.start[j,i]};fold})) 
            }
            
            expr.by.fold=rbind(unlist(lapply(1:l, FUN=function(j)
				{as.numeric(prod(expr.by.fold[,j]))})))
            d=1
            
            rank.index.to.be.removed=list()  
            l.rank.index.to.be.removed=list()
            r.rank.index.to.be.removed=list()
            both.rank.index.to.be.removed=list()
            fold.cutoff.list=list()
            cutoffs=NULL
            k=1
            for(i in 1:d){
              both.index=c()
              l.index=c()
              r.index=c()
              if(use.fold.change==TRUE){
                cutoffs= get.rna.fold.step() 
                both.index=as.numeric(lapply(cutoffs, FUN=function(f){
					cutoff.index=length(which(expr.by.fold[i,]>=f | 
						expr.by.fold[i,]<=(1/f)))}))
                r.index=as.numeric(lapply(cutoffs, FUN=function(f)
					{cutoff.index=length(which(expr.by.fold[i,]>=f))}))     
                l.index=as.numeric(lapply(cutoffs, FUN=function(f){
					cutoff.index=length(which(expr.by.fold[i,]<=(1/f)))})) 
              }else{
                step.list=rev(seq(from=1, to=exclusion.max.num, 
					by=num.genes.to.skip.each))
                cutoffs=both.index=l.index=r.index=step.list
              }
              l.rank.index.to.be.removed[[i]]=l.index
              r.rank.index.to.be.removed[[i]]=r.index
              both.rank.index.to.be.removed[[i]]=both.index
              fold.cutoff.list[[k]]=cutoffs ; k=k+1
              if(.Object@symmetric.fold==FALSE){fold.cutoff.list[[k]]=cutoffs;
				k=k+1}
            }
            
            if(.Object@symmetric.fold==FALSE){
              rank.index.to.be.removed=list(l.rank.index.to.be.removed, 
				r.rank.index.to.be.removed)
            }else{
              rank.index.to.be.removed=both.rank.index.to.be.removed
            }
                  
            .Object@rank.index.to.be.removed=rank.index.to.be.removed
            .Object@fold.cutoff.list=fold.cutoff.list
            .Object@expr.by.fold=cbind(expr.by.fold)
            return (.Object)
          })
