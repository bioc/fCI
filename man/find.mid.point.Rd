\name{find.mid.point}
\alias{find.mid.point}
\title{find the middle value of the density distribution
}
\description{find the middle value of the density distribution
}
\usage{
find.mid.point(Y)
}
\arguments{
  \item{Y}{
}
}
\details{TBD
}
\value{
\item{position }{The value the separates density into two halves}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
Y=density(sample(1:100, 50), bw=0.5)
find.mid.point(Y)
}
\keyword{ ~kwd1 }
