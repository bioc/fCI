\name{npci.gene.by.pvalues}
\alias{npci.gene.by.pvalues}
\title{find most signficantly change fCI targets}
\description{
  identify the genes that change most significantly using inverse of log ratio
  the smaller the results, the more signficant the changes.
}
\usage{
npci.gene.by.pvalues(npci.data, gene.indexes, ctr.indexes, trt.indexes)
}
\arguments{
  \item{npci.data}{a data frame containing non-zero numeric values (the data
        frame must contain more than one row and one column)
}
  \item{gene.indexes}{ the row ids of genes used for p-value calculation
}
  \item{ctr.indexes}{The wild type sample column indexes in the matrix, 
      i.e. 1,2
}
  \item{trt.indexes}{The experimental sample column indexes in the matrix, 
      i.e. 1,2
}
}
\details{TBD
}
\value{
 \item{pvalues}{a vector of pvalues}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{

  npci.data=data.frame(matrix(sample(3:100, 100*4, replace=TRUE), 100,4))
  gene.indexes=sample(1:97, 25)
  ctr.indexes=c(1,2)
  trt.indexes=c(3,4)
  npci.gene.by.pvalues(npci.data, gene.indexes, ctr.indexes, trt.indexes)
  

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
