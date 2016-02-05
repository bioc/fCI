\name{find.fci.targets}
\alias{find.fci.targets}
\title{identify differentially expressed genes
}
\description{identify differentially expressed genes
}
\usage{
find.fci.targets(.Object, wt.indexes, df.indexes, data.file, use.normalization)
}
\arguments{
  \item{.Object}{the fCI object
}
  \item{wt.indexes}{The wild type sample column indexes in the matrix, i.e. 1,2}
  \item{df.indexes}{The diseases type sample column indexes in the matrix, 
    i.e. 3,4}
  \item{data.file}{The expression matrix}
  \item{use.normalization}{boolean value whether you want the data to be 
  normalized or not}
}
\details{
TBD
}
\value{
  \item{NA }{No values will be returned}

}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
 fci=new("NPCI")
 fci.data=data.frame(matrix(sample(3:100, 1043*6, replace=TRUE), 1043,6))
 targets=find.fci.targets(fci, c(1,2,3), c(4,5,6), fci.data)
 head(show.targets(targets))
  
}
\keyword{ ~kwd1 }
