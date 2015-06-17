
setGeneric(name="find.fci.targets",
           def=function(.Object, wt.indexes, df.indexes, 
			data.file, use.normalization)
           {
             standardGeneric("find.fci.targets")
           }
)

setMethod(f="find.fci.targets",
          signature="NPCI",
          definition=function(.Object, wt.indexes, df.indexes, 
			data.file, use.normalization)
          {
            result=fCI.call.by.index(wt.indexes, df.indexes, data.file, 
				use.normalization, .Object, short.report=FALSE)
            .Object=result[[5]]
            .Object@pairwise.diff.gene.ids=result[1][[1]]
            return(.Object)
          }
)


setGeneric("show.targets",
           function(.Object)
             standardGeneric("show.targets"))  
setMethod("show.targets", "NPCI", 
          function(.Object){  
            final.results=report.target.summary(.Object@pairwise.diff.gene.ids)
            return(final.results)
          }
)

report.target.summary=function(pairwise.diff.gene.ids){
  final.results=c()
  if(length(pairwise.diff.gene.ids)>0){
    targets=Reduce(c, pairwise.diff.gene.ids)
    targets.table=table(targets)
    ratio=as.numeric(targets.table)/length(pairwise.diff.gene.ids)
    final.results=cbind(as.data.frame(targets.table), ratio)
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
        if(is.installed('psych')==FALSE){
            install.packages('psych')
        }
        if(is.installed('gtools')==FALSE){
            install.packages('gtools')
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
