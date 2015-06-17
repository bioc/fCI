\name{get.outline.index}
\alias{get.outline.index}
\title{find the outline genes of a given distribution}
\description{find the outline genes of a given distribution
}
\usage{
get.outline.index(values)
}
\arguments{
  \item{values}{
}
}
\details{TBD
}
\value{

  \item{indexes}{remove the index of values that are outliners based on the 
    t-test with alpha=0.05}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{

  values=rnorm(100)
  get.outline.index(values)

}
\keyword{ ~kwd1 }
