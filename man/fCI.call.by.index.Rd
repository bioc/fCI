\name{fCI.call.by.index}
\alias{fCI.call.by.index}
\title{top level function call to find targets based on expression data and 
  control & case indexes
}
\description{
top level function call to find targets based on expression data and control 
  & case indexes}
\usage{
fCI.call.by.index(wt.indexes, df.indexes, data.file, use.normalization = FALSE, 
  npci=NULL, short.report=TRUE)
}
\arguments{
  \item{wt.indexes}{The wild type sample column indexes in the matrix, i.e. 1,2}
  \item{df.indexes}{The diseases type sample column indexes in the matrix, 
    i.e. 3,4}
  \item{data.file}{The expression matrix}
  \item{use.normalization}{boolean value whether you want the data to be 
  normalized or not}
  \item{npci}{the fCI object}
  \item{short.report}{whether you want to have a report summary}
}
\details{TBD}
\value{
 \item{rtable }{A data frame of the detected targets}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{

wt.indexes=1:2
df.indexes=3:4
data.file=data.frame(matrix(sample(3:100, 100*4, replace=TRUE), 100,4))
use.normalization=FALSE
npci=NULL
short.report=TRUE
fCI.call.by.index(wt.indexes, df.indexes, data.file, use.normalization, 
  npci, short.report)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
