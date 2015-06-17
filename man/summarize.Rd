\name{summarize}
\alias{summarize}
\title{
result summerization}
\description{
 summerize the result after fCI computation is done
}
\usage{
summarize(.Object)
}
\arguments{
  \item{.Object}{ the class object, for exaple, fci = new("NPCI")
}
}
\details{TBD
}
\value{
\item{NA}{No values will be returned}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
 data.file=data.frame(matrix(sample(3:100, 100*4, replace=TRUE), 100,4))
wt.index=c(1,2)
df.index=c(1,3)
npci=new("NPCI")
npci@wt.index=wt.index
npci@df.index=df.index
npci@sample.data.normalized=data.file
npci=initialize(npci)
npci=normalization(npci)
npci=populate(npci)
npci=summarize(npci)
}
\keyword{ ~kwd1 }
