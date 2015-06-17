
setGeneric("figures",
           function(.Object)
             standardGeneric("figures"))					
setMethod("figures", "NPCI",
          function(.Object){
            null.data.start=.Object@null.data.start
            diff.data.start=.Object@diff.data.start
            print(max(diff.data.start))
            distance.matrix=log(as.vector(t(.Object@distance.matrix)),2)
            d=dim(diff.data.start)[2]
            par(mfrow=c(d,2))	
            for(i in 1:d){
              y.up.lim=max(density(log(null.data.start[,i],2), bw=0.5)$y, 
				density(log(diff.data.start[,i],2), bw=0.5)$y)
              x.lf.lim=min(density(log(null.data.start[,i],2), bw=0.5)$x, 
				density(log(diff.data.start[,i],2), bw=0.5)$x) - 0.5
              x.rt.lim=max(density(log(null.data.start[,i],2), bw=0.5)$x, 
				density(log(diff.data.start[,i],2), bw=0.5)$x) + 0.5
              plot(density(log(null.data.start[,i],2), bw=0.5), col="blue", 
				ylim=c(0, y.up.lim), xlim=c(x.lf.lim, x.rt.lim), 
				xlab="Ratio of Repl1/Repl2 in Log2", ylab="Density", 
				main="Density of Control(blue) and Treatment(Red)")
              lines(density(log(diff.data.start[,i],2), bw=0.5), col="red")								
            }
            plot(unlist(distance.matrix), col="red", xlab="Cutoff Index", 
				ylab="Divergence", 
				main="Divergence between Control and Treatment Groups")
            lines(unlist(distance.matrix), col="red")
            if(d>1){
              for(i in 2:d){
                plot(log(null.data.start[,1],2), log(null.data.start[,i],2), 
					col="blue", xlab="Rep11 Expression",ylab="Repl2 Expression", 
					main="Expression Correlation by Replicates")
                plot(log(diff.data.start[,1],2), log(diff.data.start[,i],2), 
					col="red", xlab="Rep11 Expression", ylab="Repl2 Expression", 
					main="Expression Correlation by Replicates")
              }
            }		
            if(length(.Object@fold.cutoff.list)>1){
              x=as.numeric(.Object@fold.cutoff.list[[1]])
              y=as.numeric(.Object@fold.cutoff.list[[2]])
              x=sort(x, decreasing=FALSE)
              y=sort(y, decreasing=FALSE)								
              z=log(.Object@distance.matrix,2)
              nbcol = 100
              color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
              zcol  = cut(z, nbcol)
              persp3d(x, y, z, 
                      theta=50, phi=25, expand=0.75, col=color[zcol], 
					  shade = 0.4, ticktype="detailed", xlab="", ylab="time", 
					  zlab="",axes=TRUE, nticks=5)
            }
        }
)



setGeneric("venndiagram",
           function(.Object)
             standardGeneric("venndiagram"))					
setMethod("venndiagram", "NPCI",
          function(.Object){	
            
            d=dim(.Object@diff.data.start)[2]
            num.pairs=length(.Object@pairwise.diff.gene.ids)
            for(i in 1:d){
              for(j in 1:round(num.pairs/4)){
                diff.gene.ids=list()
                n=1
                for(k in (4*j-3):(min(4*j, num.pairs))){
                  diff.gene.ids[[n]]=.Object@pairwise.diff.gene.ids[[k]][i][[1]]
                  n=n+1
                }	
                npci.venn.diagram(diff.gene.ids,i,(4*j-3))
              }
            }
        }
)


setGeneric("summarize",
           function(.Object)
             standardGeneric("summarize"))					
setMethod("summarize", "NPCI",
          function(.Object){
            distance.matrix=.Object@distance.matrix
            distance.matrix=as.vector(t(distance.matrix))
            min.dist.index=which(distance.matrix==min(distance.matrix))
            min.dist.index=min.dist.index[1]
            
            fold.combinations=do.call(expand.grid, .Object@fold.cutoff.list)	
            combinations=get.rank.combinations(.Object@rank.index.to.be.removed,
				.Object@symmetric.fold)
            
            diff.gene.num=as.numeric(combinations[min.dist.index,])
            diff.gene.ids=list()
            
            d=if(.Object@symmetric.fold==TRUE){length(diff.gene.num)}
				else{(length(diff.gene.num)/2)}
            if(!(is.na(diff.gene.num))){
              for(i in 1:d){
                diff.gene.ids[[i]]=if(.Object@symmetric.fold==TRUE){
                  this.expr.by.fold=.Object@expr.by.fold[i,]
                  this.expr.by.fold[which(this.expr.by.fold<1)]=
					1/(this.expr.by.fold[which(this.expr.by.fold<1)])
                  order(this.expr.by.fold, decreasing=TRUE)[0:diff.gene.num[i]]
                }else{
                  union(order(.Object@diff.data.start[,i], 
					decreasing=FALSE)[0:diff.gene.num[2*i-1]],  
                        order(.Object@diff.data.start[,i], 
						decreasing=TRUE )[0:diff.gene.num[2*i ]])	
                }
              }
            }
            .Object@diff.gene.ids=diff.gene.ids
            fold.change=as.numeric(as.matrix(
				fold.combinations[min.dist.index,]))
            result=c(diff.gene.num, fold.change, round(min(distance.matrix), 8))
            .Object@result=result
            
            if(FALSE & dim(.Object@null.data.start)[2]==1 
				& .Object@symmetric.fold==TRUE 
				& dim(.Object@sample.data.normalized)[2]>=
				length(unique(.Object@ctr.indexes, .Object@trt.indexes))){
              index.to.be.removed=npci.index.to.be.removed(.Object@expr.by.fold,
				1, TRUE, .Object@rank.index.to.be.removed[[1]][which(round(
				as.numeric(.Object@fold.cutoff.list[[1]]),1)==fold.change)],1,1)
				indexes.reconsidered=npci.index.reconsidered(
				.Object@sample.data.normalized, .Object@expr.by.fold, 
				.Object@null.data.start, .Object@diff.data.start, 
				index.to.be.removed, 
				.Object@ctr.indexes, 
				.Object@trt.indexes, 
				fold.change, 
				fold.change)
              if(length(indexes.reconsidered)>0){
                .Object@indexes.reconsidered=indexes.reconsidered
              }
            }
            
            return (.Object)
    }
)



npci.venn.diagram<-function(diff.gene.ids, i=1, k=1){
  
  if(length(diff.gene.ids)==4){
    set1=diff.gene.ids[[1]]
    set2=diff.gene.ids[[2]]
    set3=diff.gene.ids[[3]]
    set4=diff.gene.ids[[4]]
    venn.plot <-draw.quad.venn(length(set1), length(set2), 
		length(set3), length(set4),
        length(intersect(set1, set2)), length(intersect(set1, set3)), 
		length(intersect(set1, set4)), length(intersect(set2, set3)), 
		length(intersect(set2, set4)), length(intersect(set3, set4)),
        length(intersect(intersect(set1, set2), set3)), 
		length(intersect(intersect(set1, set2), set4)), 
		length(intersect(intersect(set1, set3), set4)), 
		length(intersect(intersect(set2, set3), set4)), 
        length(intersect(intersect(intersect(set1, set2), set3), set4)),
        category = c("First", "Second", "Third", "Fourth"),
        fill = c("orange", "red", "green", "blue"),
			lty = "dashed",
			cex = 2,
			cat.cex = 2,
			cat.col = c("orange", "red", "green", "blue")
		)
   }				
  
  if(length(diff.gene.ids)==3){
    set1=diff.gene.ids[[1]]
    set2=diff.gene.ids[[2]]
    set3=diff.gene.ids[[3]]
    venn.plot <-draw.triple.venn(length(set1), length(set2), length(set3), 
        length(intersect(set1, set2)), 
		length(intersect(set2, set3)), 
		length(intersect(set1, set3)), 
        length(intersect(intersect(set1, set2), set3)),
        category = c("First", "Second", "Third"),
        fill = c("orange", "red", "green"),
        lty = "dashed",
        cex = 2,
        cat.cex = 2,
        cat.col = c("orange", "red", "green")
    )
  }				
  if(length(diff.gene.ids)==2){
    set1=diff.gene.ids[[1]]
    set2=diff.gene.ids[[2]]
    venn.plot <-draw.pairwise.venn(length(set1), length(set2),
        length(intersect(set1, set2)),
        category = c("First", "Second"),
        fill = c("orange", "red"),
			lty = "dashed",
			cex = 2,
			cat.cex = 2,
			cat.col = c("orange", "red")
		)
  }
  
  tiff(filename = paste("Quad_Venn_diagram_", i, "_", k, "_", 
	(k+length(diff.gene.ids)-1), ".tiff", sep=""), 
	compression = "lzw");
  grid.draw(venn.plot);
  dev.off();
  
}

