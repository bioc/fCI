\name{divergence.multivariate.distributions}
\alias{divergence.multivariate.distributions}
\title{estimate fCI divergence for given samples of aritrary dimensions
}
\description{
estimate fCI divergence for given samples of aritrary dimensions}
\usage{
divergence.multivariate.distributions(null.data, diff.data, choice = 2)
}
\arguments{
  \item{null.data}{
   the empirical null dataset (a dataframe of none-zero ratio values)
}
  \item{diff.data}{
  the case-control dataset (a dataframe of none-zero ratio values)
  
  }
  \item{choice}{
  choice=1 => cross entropy
  choice=2 => Helligan distance
  choice=3 => KL distance
}
}
\details{TBD
}
\value{
  \item{divergences}{The estimated divergence given control-control and case-
  control expression ratios}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{

null.data=data.frame(matrix(sample(seq(from=0.1,to=10, by=0.01), 100), 100,1))
diff.data=data.frame(matrix(sample(seq(from=0.1,to=10, by=0.01), 100), 100,1))
divergence.multivariate.distributions(null.data, diff.data, choice = 2)

}

\keyword{ ~kwd1 }
