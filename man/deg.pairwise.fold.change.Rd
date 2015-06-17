\name{deg.pairwise.fold.change}
\alias{deg.pairwise.fold.change}
\title{find targets that have a consistent fold change in the same direction 
  (either up- or down-regulation)}
\description{
  find targets that have a consistent fold change in the same direction
}
\usage{
deg.pairwise.fold.change(pairwise.wt.up.down.fold, pairwise.df.up.down.fold, 
  d = 1, min.fold = 1.2)
}
\arguments{
  \item{pairwise.wt.up.down.fold}{a list of numeric values representing
        the fold changes between control replicates for every gene
}
  \item{pairwise.df.up.down.fold}{a list of numeric values representing
        the fold changes between case and control replicates for every gene
}
  \item{d}{ the dimensionality of the database, if the dataset is from 
           proteogenomics, then d=2
}
  \item{min.fold}{minimum fold change to declare a gene to be dysregulated, by 
      default, min.fold=2
}
}
\details{TBD
}
\value{
\item{expression ratio}{a dataframe of fCI gene expression ratios (folds) with 
   none zero values defined by given control-control index (i.e. 1 & 2) and 
   control-case index (i.e. 3&4)}
}
\references{http://software.steenlab.org/fCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{

wt.fold.changes=list(c(1.2,1.3,1.5,1.6))
df.fold.changes=list(c(1.1,1.3,1.4,1.6))
deg.pairwise.fold.change(wt.fold.changes,df.fold.changes)
}
\keyword{ ~kwd1 }
