\name{deseq.median.ratio.normalization}
\alias{deseq.median.ratio.normalization}
\title{ data matrix normalization method }
\description{
normalize expression matrix by first replicate's median gene expression values
}
\usage{
deseq.median.ratio.normalization(npci.data)
}
\arguments{
  \item{npci.data}{a data frame containing non-zero numeric values (the data
        frame must contain more than one row and one column) }
}
\details{TBD
}
\value{
  \item{data.frame}{a new dataframe with each column having the same median 
  value}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{

udata=data.frame(matrix(sample(3:100, 6*4), 6,4))
normalized.udata=deseq.median.ratio.normalization(udata)

}

\keyword{ ~normalization }

