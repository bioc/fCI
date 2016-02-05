
setGeneric(name="find.fci.targets",
           def=function(.Object, wt.indexes, df.indexes, 
                        data.file, use.normalization=TRUE)
           {
             standardGeneric("find.fci.targets")
           }
)

setMethod(f="find.fci.targets",
          signature="NPCI",
          definition=function(.Object, wt.indexes, df.indexes, 
                              data.file, use.normalization=TRUE)
          {
            
            if(is.list(wt.indexes) & is.list(df.indexes)){
            #if(use.normalization==TRUE){
            #  wt.indexes=list(wt.indexes[1:(length(wt.index)/2)],
            #                  wt.indexes[(length(wt.index)/2+1):length(wt.index)])
            #  wt.indexes=list(df.indexes[1:(length(df.index)/2)],
            #                  df.indexes[(length(df.index)/2+1):length(df.index)])              
              
              if(is.matrix(data.file) | is.data.frame(data.file)){
                .Object@sample.data.normalized=data.file
              }else{

                .Object@sample.data.file=data.file
                .Object=initialize(.Object)
              }
              
              .Object@ctr.indexes=unlist(wt.indexes)
              .Object@trt.indexes=unlist(df.indexes)
              
              wt.index.more=combinations(length(wt.indexes[[1]]),2, 
                                         v=wt.indexes[[1]], repeats.allowed=FALSE)
              wt.index.more.1=(lapply(1:dim(wt.index.more)[1], function(x) 
                as.numeric(wt.index.more[x,])))
              
              wt.index.more=combinations(length(wt.indexes[[2]]),2, 
                                         v=wt.indexes[[2]], repeats.allowed=FALSE)
              wt.index.more.2=(lapply(1:dim(wt.index.more)[1], function(x) 
                as.numeric(wt.index.more[x,])))
              
              
              df.index.more=combinations(length(df.indexes[[1]]),2, 
                                         v=df.indexes[[1]], repeats.allowed=FALSE)
              df.index.more.1=(lapply(1:dim(df.index.more)[1], function(x) 
                as.numeric(df.index.more[x,])))
              
              df.index.more=combinations(length(df.indexes[[2]]),2, 
                                         v=df.indexes[[2]], repeats.allowed=FALSE)
              df.index.more.2=(lapply(1:dim(df.index.more)[1], function(x) 
                as.numeric(df.index.more[x,])))    
              
              
              pairwise.index=list()
              cnt=1
              l=1
              for(i in 1:length(wt.index.more.1)){
                for(j in 1:length(wt.index.more.2)){
                  for(k in 1:length(df.index.more.1)){
                    for(p in 1:length(df.index.more.2)){
                      if(cnt<1000){
                        wt.index.in.list=list(a=wt.index.more.1[[i]][1],
                                              b=wt.index.more.1[[i]][2],
                                              d=wt.index.more.2[[j]][1],
                                              f=wt.index.more.2[[j]][2])
                        df.index.in.list=list(a=wt.index.more.1[[i]], 
                                              d=df.index.more.1[[k]],
                                              b=wt.index.more.2[[j]], 
                                              f=df.index.more.2[[p]])
                        
                        result.lists=deg.up.down.info(wt.index.in.list, 
                                                      df.index.in.list, .Object, TRUE)
                        
                        good.degs=deg.pairwise.fold.change(result.lists[[3]], 
                                                           result.lists[[4]],
                                                           d=2,
                                                           1.1)
                        deg.occupancy=pairwise.change.occupancy(good.degs, result.lists[[6]],
                                                                result.lists[[7]], 0.5)
                        
                        #for(m in 1:length(result.lists[[6]])){
                        #  pairwise.index[[l]]=result.lists[[6]][[m]]
                        #  l=l+1
                        #}
                        pairwise.index[[l]]=rownames(.Object@sample.data.normalized)[deg.occupancy[[1]]]
                        l=l+1
                        cnt=cnt+1
                      }
                    }
                  }
                }
              }
              
              .Object@pairwise.diff.gene.ids=pairwise.index
              
              
            }else{
              
              result=fCI.call.by.index(unlist(wt.indexes), unlist(df.indexes), data.file, 
                                       use.normalization, .Object, short.report=FALSE)
    
              .Object=result[[5]]  ## which is actually the object fci
              .Object@pairwise.diff.gene.ids=result[1][[1]]  ## actually pairwise.index
              .Object@ctr.indexes=unlist(wt.indexes)
              .Object@trt.indexes=unlist(df.indexes)
            }
            return(.Object)
          }
)


setGeneric("show.targets",
           function(.Object)
             standardGeneric("show.targets"))  
setMethod("show.targets", "NPCI", 
          function(.Object){  
            # if(length(.Object@pairwise.diff.gene.ids)==0){
            #  .Object@pairwise.diff.gene.ids=.Object@diff.gene.ids
            # }
            result=NULL
            if(length(.Object@pairwise.diff.gene.ids)==0 & length(.Object@diff.gene.ids)==0){
              print("No differentially expressed genes are found!")
            }else{
              if(length(.Object@pairwise.diff.gene.ids)==0){
                .Object@pairwise.diff.gene.ids=.Object@diff.gene.ids
              }
              
              result=report.target.summary(.Object@pairwise.diff.gene.ids)
              index.of.DEGs=unlist(lapply(as.vector(result[,1]), function(n)
                {which(rownames(.Object@sample.data.normalized)==n)}))
              
              mean.a=round(apply(cbind(.Object@sample.data.normalized[index.of.DEGs, 
                                                                      unlist(.Object@ctr.indexes)]), 1, mean),3)
              mean.b=round(apply(cbind(.Object@sample.data.normalized[index.of.DEGs, 
                                                                      unlist(.Object@trt.indexes)]), 1, mean),3)
              log.fc=round(log(mean.b/mean.a,2),3)
              result=data.frame(cbind(as.vector(result[,1]), mean.a, mean.b, 
                                      log.fc, round(as.numeric(result[,2]),3)))
              colnames(result)=c("DEG_Names", "Mean_Control", "Mean_Case", "Log2_FC", "fCI_Prob_Score")
              if(dim(result)[1]>0){rownames(result)=1:dim(result)[1]}
              
              #   cat("The optimal cutoff fold change is ", as.vector(
              #    .Object@fold.cutoff.list[[1]][which(as.vector(.Object@distance.matrix)==
              #          min(as.vector(.Object@distance.matrix)))]), "\n")
              cat("A total of ", dim(result)[1], 
                  " genes were identified as differentially expressed.\n\n")            
              
            }
            return(result)
          }
)



report.target.summary=function(pairwise.diff.gene.ids){
  final.results=c()
  if(length(pairwise.diff.gene.ids)>0){
    targets=Reduce(c, pairwise.diff.gene.ids)
    targets.table=table(targets)
    ratio=as.numeric(targets.table)/length(pairwise.diff.gene.ids)
    DEGs=rep(0, length(targets))
    DEGs[which(ratio>0.5)]=1
    final.results=cbind(as.data.frame(targets.table), ratio)
    final.results=final.results[which(DEGs==1), c(1,3)]
    colnames(final.results)=c("DEG_Name", "fCI_Probability_Score")
  }
  return(final.results)
}


setGeneric("call.npci",
           function(.Object)
             standardGeneric("call.npci"))	
setMethod("call.npci", "NPCI", 
          function(.Object){					
            wt.comb.list=.Object@wt.comb
            df.comb.list=.Object@df.comb
            wt.combinations=do.call(expand.grid, wt.comb.list)
            df.combinations=do.call(expand.grid, df.comb.list)
            
            pairwise.diff.gene.ids=list()
            k=1
            for(i in 1:dim(wt.combinations)[1]){
              for(j in 1:dim(df.combinations)[1]){
                wt.comb=as.numeric(wt.combinations[i,])
                df.comb=as.numeric(df.combinations[j,])							
                .Object@wt.index=wt.comb
                .Object@df.index=df.comb
                .Object=populate(.Object)
                .Object=compute(.Object)
                .Object=summarize(.Object)
                pairwise.diff.gene.ids[[k]]=.Object@diff.gene.ids
                k=k+1
              }
            }
            .Object@pairwise.diff.gene.ids=pairwise.diff.gene.ids
            return(.Object)
})


setGeneric("compute",
           function(.Object)
             standardGeneric("compute"))	
setMethod("compute", "NPCI", #signature(obj="NPCI"),
          
    function(.Object){
            
        if(is.installed('FNN')==FALSE){
             install.packages('FNN')
        }
        if(is.installed('gtools')==FALSE){
          install.packages('gtools')
        }
        if(is.installed('psych')==FALSE){
            install.packages('psych')
        }
        distance.matrix=get.npci.distance.matrix(.Object@sample.data.normalized,
            .Object@null.data.start, 
            .Object@diff.data.start, 
            .Object@method.option,
            .Object@rank.index.to.be.removed,
            .Object@expr.by.fold,
            .Object@ctr.indexes,
            .Object@trt.indexes,
            FALSE,
            .Object@symmetric.fold,
            .Object@fold.cutoff.list)
            fold.cutoff.list=.Object@fold.cutoff.list
            distance.matrix=matrix(unlist(distance.matrix), 
				nrow = length(fold.cutoff.list[[1]]), byrow = TRUE, 
				dimnames=fold.cutoff.list)  
            .Object@distance.matrix=distance.matrix
            return (.Object)
})
