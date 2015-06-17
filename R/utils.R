

two.sample.log.ratio<-function(a,b){
  
  return(1/(abs(log(mean(as.numeric(b))/mean(as.numeric(a)),2))))
}

two.sample.permutation.test<-function(a,b){
  
  # Observed difference
  diff.observed = mean(as.numeric(b)) - mean(as.numeric(a))
  c=c(a,b)
  num.comb=choose(length(a)+length(b), length(a))
  number_of_permutations = if(num.comb>1000){1000}else{num.comb}
  
  diff.random = NULL
  for (i in 1 : number_of_permutations) {
    
    # Sample from the combined dataset
    a.random = sample (c, length(a), TRUE)
    b.random = sample (c, length(b), TRUE)
    
    # Null (permuated) difference
    diff.random[i] = mean(b.random) - mean(a.random)
  }
  
  pvalue = sum(abs(diff.random) >= abs(diff.observed)) / number_of_permutations
  return(pvalue)
}


get.rna.fold.step<-function(){
  fold.lower.end=seq(from=11, to=30, by=1)/10	 # seq(from=1.1, to=3, by=0.1)
  fold.middle=seq(from=32, to=50, by=2)/10
  fold.upper.end=seq(from=6, to=10, by=1)
  fold.cutoffs=(c(fold.lower.end, fold.middle, fold.upper.end))
  return(fold.cutoffs)
}

get.protein.fold.step<-function(){
  fold.lower.end=seq(from=11, to=20, by=1)/10 #fold.lower.end[[2]]==1.2 is FALSE
  fold.middle=seq(from=22, to=30, by=2)/10
  fold.upper.end=seq(from=4, to=10, by=1)
  fold.cutoffs=(c(fold.lower.end, fold.middle, fold.upper.end))
  return (fold.cutoffs)  
}

get.fold.large.step<-function(){
  fold.lower.end=seq(from=11, to=20, by=2)/10
  fold.middle=seq(from=22, to=35, by=5)/10
  fold.upper.end=seq(from=4, to=10, by=1)
  fold.cutoffs=(c(fold.lower.end, fold.middle, fold.upper.end))
  return(fold.cutoffs)
}

intersect.of.lists <- function(vectorlist) {
  return(Reduce(intersect, vectorlist))
}


get.outline.index<-function(values){
  
  x.order=order(values, decreasing=FALSE)
  x=sort(values, decreasing=FALSE)
  x.quantile=quantile(x, c(0.005, 0.1, 0.9, 0.995))
  x.mean=mean(x)
  x.sd=sd(x)
  outline.sd.cutoff=3.37  
  
  removed.index=c()
  for(i in 1:length(x)){
    if((abs(x[i]-mean(x))/sd(x))<(-1*outline.sd.cutoff) | 
		(abs(x[i]-mean(x))/sd(x))>outline.sd.cutoff){
      removed.index=c(removed.index, x.order[i])
    }
  }
  return(removed.index)
}


find.mid.point=function(Y){
  # Y is the density of a vector Y=density(vals,bw=0.5)
  best.point=0
  distance=1
  for(x in Y$x[2:length(Y$x)]){
    xt <- diff(Y$x[Y$x<=x])
    yt <- rollmean(Y$y[Y$x<=x],2)
    x.sum=sum(xt*yt)
    if(abs(x.sum-0.5)<distance){
      best.point=x
      distance=abs(x.sum-0.5)
    }
  }
  return(best.point)
}

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 



