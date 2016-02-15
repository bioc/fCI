\name{NPCI-class}
\Rdversion{1.1}
\docType{class}
\alias{NPCI-class}

\title{Class \code{"NPCI"}}
\description{
The main Class that defines the slots values
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NPCI", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{sample.data.file}:}{Object of class \code{"character"} ~~ }
    \item{\code{distance.matrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{sample.data.normalized}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{attr.info}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{null.data.start}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{diff.data.start}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{expr.by.fold}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{fold.cutoff.list}:}{Object of class \code{"list"} ~~ }
    \item{\code{rank.index.to.be.removed}:}{Object of class \code{"list"} ~~ }
    \item{\code{diff.gene.ids}:}{Object of class \code{"list"} ~~ }
    \item{\code{wt.index}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{df.index}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ctr.indexes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{trt.indexes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{method.option}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{use.ratio}:}{Object of class \code{"logical"} ~~ }
    \item{\code{percent.genes.to.scan}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{num.genes.to.skip.each}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{use.fold.change}:}{Object of class \code{"logical"} ~~ }
    \item{\code{wt.comb}:}{Object of class \code{"list"} ~~ }
    \item{\code{df.comb}:}{Object of class \code{"list"} ~~ }
    \item{\code{diff.ids}:}{Object of class \code{"list"} ~~ }
    \item{\code{result}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{indexes.reconsidered}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{center.by.gaussian.kernel}:}{Object of class \code{"logical"} ~~ }
    \item{\code{symmetric.fold}:}{Object of class \code{"logical"} ~~ }
    \item{\code{pairwise.diff.gene.ids}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "NPCI" in the signature.
}
\references{http://software.steenlab.org/NPCI/}
\author{Shaojun Tang}
\note{TBD}

\seealso{TBD}
\examples{
showClass("NPCI")
}
\keyword{classes}
