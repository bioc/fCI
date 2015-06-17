\name{pairwise.change.occupancy}
\alias{pairwise.change.occupancy}
\title{find the targets whose fold changes occur consistently (upregulated or 
  downregulated) in all fCI analysis}
\description{
  find the targets whose fold changes occur consistently 
  (upregulated or downregulated) in all fCI analysis 
}
\usage{
pairwise.change.occupancy(common.ids, pairwise.index,
  pairwise.up.down, target.ratio)
}
\arguments{
  \item{common.ids}{ the gene ids that are differentially expressed
}
  \item{pairwise.index}{ a list of the genes ids that differentially expressed in 
    each of the fCI analysis
}
  \item{pairwise.up.down}{ a list of up regulatio (+1) or downregulation (-1) 
   for each gene in fCI analysis
}
  \item{target.ratio}{ the expected fold changes 
}
}
\details{TBD
}
\value{
  \item{consistent targets}{Gene (indexes) that are consistently changed in 
    fCI pairwise analysis}
  \item{direction}{Gene (indexes) that are consistently upregulated (if < 0) or
    upregulated (if > 0)}
      
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
  common.ids=6:13
  pairwise.index=list(c(4:13), c(6:15))
  pairwise.up.down=list(c(sample(c(-1,1), 10, replace=TRUE)),
                        c(sample(c(-1,1), 10, replace=TRUE)))
  target.ratio=0.5
  pairwise.change.occupancy(common.ids, pairwise.index,
    pairwise.up.down, target.ratio)

}
\keyword{ ~kwd1 }
