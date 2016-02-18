\name{setfCI}
\alias{setfCI}
\title{ the generic function 'setfCI' for s4 class
}
\description{the generic function 'setfCI' for s4 class
}
\usage{
setfCI(.Object, wt.index, df.index, fold.cutoff.list, 
  center.distribution)
}
\arguments{
  \item{.Object}{the fCI object}
  \item{wt.index}{the control sample column ids,
    such as c(1,2)}
  \item{df.index}{the case sample column ids, 
    such as c(1,2)}
  \item{fold.cutoff.list}{the predefined fold change cut-off
    such as list(seq(from=1.1, to=3.0, by=0.1))}
  \item{center.distribution}{a boolean value showing that if
    the users want to center the distribution or not}
  
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
 fci=setfCI(fci, 7:8, 11:12, seq(from=1.1,to=3,by=0.1), TRUE)
}
\keyword{ ~kwd1 }