\name{deg.up.down.info}
\alias{deg.up.down.info}
\title{find targets and their detailed expression changes
}
\description{
given expression matrix, find targets and their detailed expression changes
}
\usage{
deg.up.down.info(wt.index.in.list, df.index.in.list, npci, 
use.normalization = FALSE, target.ratio = 0.5)
}
\arguments{
  \item{wt.index.in.list}{  a list of numeric values representing
        the column indexes for control samples
}
  \item{df.index.in.list}{  a list of numeric values representing
        the column indexes for experimental samples
}
  \item{npci}{
        the object npci
}
  \item{use.normalization}{ a boolean value indicating if the normalization
    will be applied or not
}
  \item{target.ratio}{
    a numeric value indicating the expected fold changes, i.e, 1.5
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
  print("TBC")
}
\keyword{ ~kwd1 }
