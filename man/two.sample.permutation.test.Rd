\name{two.sample.permutation.test}
\alias{two.sample.permutation.test}
\title{perform permuation test on two vectors}
\description{perform permuation test on two vectors
}
\usage{
two.sample.permutation.test(a, b)
}
\arguments{
  \item{a}{ a vector of numeric values (value must be greater than 0)
}
  \item{b}{ a vector of numeric values (value must be greater than 0)
}
}
\details{TBD
}
\value{
\item{pvalue}{the pvalue of permutation test}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
  
  two.sample.permutation.test(sample(1:100, 20), sample(5:104, 20))

}
\keyword{ ~kwd1 }
