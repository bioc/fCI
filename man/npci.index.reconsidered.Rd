\name{npci.index.reconsidered}
\alias{npci.index.reconsidered}
\title{find targets that have little evidence to be differentially expressed}
\description{
  the function will be depreciated
}
\usage{
npci.index.reconsidered(npci.data, expr.by.fold, null.data.start, diff.data.start, gene.indexes, ctr.indexes, trt.indexes, left.fold, right.fold)
}
\arguments{
  \item{npci.data}{a data frame containing non-zero numeric values (the data
        frame must contain more than one row and one column)
}
  \item{expr.by.fold}{a 1xN matrix of case-control fold changes for every gene
      of the total N genes
}
  \item{null.data.start}{a Nx1 matrix of control-control fold changes
}
  \item{diff.data.start}{a Nx1 matrix of case-control fold changes
}
  \item{gene.indexes}{ the genes used for differential expression analysis.
}
  \item{ctr.indexes}{ the control sample column indexes
}
  \item{trt.indexes}{ the case sample column indexes
}
  \item{left.fold}{ the minimum fold changes for downregulation
}
  \item{right.fold}{ the minimum fold changes for upregulation
}
}
\details{TBD
}
\value{
\item{values }{genes wrongly considered as differentially expressed}
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
npci=compute(npci)
npci=summarize(npci)

npci.data=npci@sample.data.normalized
null.data.start=npci@null.data.start
diff.data.start=npci@diff.data.start
choice=2
rank.index.to.be.removed=npci@rank.index.to.be.removed
expr.by.fold=npci@expr.by.fold

ctr.indexes=1:2
trt.indexes=3:4
use.intersect=FALSE
symmetric.fold=TRUE
fold.cutoff.list=npci@fold.cutoff.list
gene.indexes=npci@diff.gene.ids
left.fold=2
right.fold=2


}
\keyword{ ~kwd1 }
