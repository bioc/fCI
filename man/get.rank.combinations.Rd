\name{get.rank.combinations}
\alias{get.rank.combinations}
\title{fold change values}
\description{
  identify the fold change value indexes beyond the fCI estimation 
}
\usage{
  get.rank.combinations(rank.index.to.be.removed, symmetric.fold)
}
\arguments{
  \item{rank.index.to.be.removed}{a list of integers representing the genes
        to be removed because it exceeds the predefined fold change, i.e 
        1.2 fold
}
  \item{symmetric.fold}{ a boolean value indicating the upregulation and
      downregulation are treatedly equally
}
}
\details{TBD
}
\value{
   \item{combinations}{a data frame of gene indexes}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
  rank.index.to.be.removed=list(sample(1:100, 20))
  symmetric.fold=TRUE
  get.rank.combinations(rank.index.to.be.removed, symmetric.fold)

}
\keyword{ ~kwd1 }
