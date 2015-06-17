\name{npci.venn.diagram}
\alias{npci.venn.diagram}
\title{generate venn diagram for multiple fCI analysis}
\description{
  plot the overlap differentially expressed genes by pairwise fCI analysis
}
\usage{
npci.venn.diagram(diff.gene.ids, i = 1, k = 1)
}
\arguments{
  \item{diff.gene.ids}{ gene ids for genes that are differentially expressed
}
  \item{i}{ number of comparisons for fCI analysis, i,e  1 or 2
}
  \item{k}{ number of genes for fCI analysis
}
}
\details{TBD
}
\value{
\item{figure}{the venn diagram plot}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{

  targets.run1=c(2:10)
  targets.run2=c(1:8)
  targets.run3=c(6:12)
  diff.gene.ids=list(targets.run1, targets.run2, targets.run3)
  npci.venn.diagram(diff.gene.ids)
  
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
