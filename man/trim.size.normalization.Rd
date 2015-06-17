\name{trim.size.normalization}
\alias{trim.size.normalization}
\title{normalize gene expression by exluding genes on the top 5 and bottom
5 percentage}
\description{normalize gene expression by exluding genes on the top 5 and 
bottom 5 percentage
}
\usage{
trim.size.normalization(sample.data)
}
\arguments{
  \item{sample.data}{a data frame of gene expression (noen-zero) with columns 
  being the sample and rows being genes
}
}
\details{TBD
}
\value{
\item{dataframe}{a data frame where column values were normalized by all genes
except the top 5 percent and bottom 5 percent genes}
}
 
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
  sample.data=data.frame(matrix(sample(3:100, 100*4, replace=TRUE), 100,4))
  trim.size.normalization(sample.data)
}
\keyword{ ~kwd1 }

