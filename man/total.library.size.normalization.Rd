\name{total.library.size.normalization}
\alias{total.library.size.normalization}
\title{normalize the gene expression based on the library size (summation)
  of the first sample replicate}
\description{
normalize the gene expression based on the library size (summation) of the
first sample replicate
}
\usage{
total.library.size.normalization(sample.data)
}
\arguments{ 
  \item{sample.data}{a data frame of gene expression (noen-zero) with columns 
  being the sample and rows being genes
}
}
\details{TBD
}
\value{
\item{dataframe}{a data frame where column values were normalized by total 
library size}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
  sample.data=data.frame(matrix(sample(3:100, 100*4, replace=TRUE), 100,4))
  total.library.size.normalization(sample.data)

}
\keyword{ ~kwd1 }
