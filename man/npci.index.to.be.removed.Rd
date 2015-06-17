\name{npci.index.to.be.removed}
\alias{npci.index.to.be.removed}
\title{gene indexes that will be considered as targets}
\description{
  This function will be depreciated.
}
\usage{
npci.index.to.be.removed(expr.by.fold, d, symmetric.fold, max.rank, 
l.max.rank, r.max.rank)
}
\arguments{
  \item{expr.by.fold}{ a 1xN matrix of fold change between case and control for
    every genes in N genes
}
  \item{d}{ the dimension of the data, if RNA-Seq or LC-MS/MS data, d=1
}
  \item{symmetric.fold}{ a booleam valuable indicating whether to use the same
   fold change cutoff for upregulation and downregulation
}
  \item{max.rank}{ the maximum fold change, i.e  3 fold
}
  \item{l.max.rank}{ the maximum fold change for downregulation, i.e  1.5 fold
}
  \item{r.max.rank}{ the maximum fold change for upregulation, i.e  1.5 fold
}
}
\details{TBD
}
\value{
\item{indexes}{gene (indexes) considered as differentially expressed}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
 print("Function to be discarded!")
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
