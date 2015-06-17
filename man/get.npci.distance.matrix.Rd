\name{get.npci.distance.matrix}
\alias{get.npci.distance.matrix}
\title{generate the divergence estimation based of fold change cutoff values
}
\description{
  generate the divergence estimation based of fold change cutoff values
}
\usage{
get.npci.distance.matrix(npci.data, null.data.start, diff.data.start, choice = 2, rank.index.to.be.removed, expr.by.fold, ctr.indexes, trt.indexes, use.intersect = FALSE, symmetric.fold = TRUE, fold.cutoff.list)
}
\arguments{
  \item{npci.data}{
}
  \item{null.data.start}{
}
  \item{diff.data.start}{
}
  \item{choice}{
}
  \item{rank.index.to.be.removed}{
}
  \item{expr.by.fold}{
}
  \item{ctr.indexes}{
}
  \item{trt.indexes}{
}
  \item{use.intersect}{
}
  \item{symmetric.fold}{
}
  \item{fold.cutoff.list}{
}
}
\details{TBD
}
\value{

  \item{divergence }{A matrix of computed divergences}

}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{


data.file=data.frame(matrix(sample(3:100, 100*4, replace=TRUE), 100,4))
wt.index=c(1,2)
df.index=c(1,3)
npci=new("NPCI")
npci@wt.index=wt.index
npci@df.index=df.index
npci@sample.data.normalized=data.file
npci=initialize(npci)
npci=normalization(npci)
npci=populate(npci)

null.data.start=npci@null.data.start
diff.data.start=npci@diff.data.start
choice=2
rank.index.to.be.removed=npci@rank.index.to.be.removed
expr.by.fold=npci@expr.by.fold
ctr.indexes=npci@wt.index
trt.indexes=npci@df.index
use.intersect=FALSE
symmetric.fold=TRUE
fold.cutoff.list=npci@fold.cutoff.list

get.npci.distance.matrix(npci.data, null.data.start, diff.data.start, 
  choice = 2, rank.index.to.be.removed, expr.by.fold, ctr.indexes, trt.indexes, 
  use.intersect, symmetric.fold, fold.cutoff.list)

}
\keyword{ ~kwd1 }
