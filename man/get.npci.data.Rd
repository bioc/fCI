\name{get.npci.data}
\alias{get.npci.data}
\title{return a fCI object given the gene expression data}
\description{
return a fCI object given the gene expression data}
\usage{
get.npci.data(sample.data.normalized, wt.index, df.index)
}
\arguments{
  \item{sample.data.normalized}{
}
  \item{wt.index}{
}
  \item{df.index}{
}
}
\details{TBD
}
\value{
\item{expression ratio}{a dataframe of fCI gene expression ratios (folds)
  defined by control-control index and control-case index}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
  sample.data.normalized=data.frame(matrix(sample(3:100, 100*4, replace=TRUE),
    100,4))
  wt.index=c(1,2)
  df.index=c(1,3)
  get.npci.data(sample.data.normalized, wt.index, df.index)
}
\keyword{ ~kwd1 }
