\name{triplex.3D}
\alias{triplex.3D}

\title{
Triplex visualization, 3D representation 
}

\description{
This function visualizes a \code{\link{TriplexViews}} object as a 3D model.
Its structure can be drawn with automatic optimalizations.
}

\details{
The input \code{\link{TriplexViews}} object is required to provide additional
algorithm options (see \code{\link{triplex.search}}). These are used for
proper computation of triplex alignment.

An example of a graphical output corresponding to a triplex type 3 with DNA
sequence "GGAAAGCAATGCCAGGCAGGG" is shown in the following figure

\if{text}{Figure: Example 3D visualization (see HTML or Latex version of this page)}
\figure{triplex3d.png}
}

\usage{
triplex.3D(triplex, opt = TRUE, A.col = "red", T.col = "brown", 
           G.col = "green", C.col = "blue", bbone.col = "violet", 
           bgr.col = "white", bbone.n = 20)
}

\arguments{
  \item{triplex}{
    \code{\link{TriplexViews}} object including only one triplex.
  }
  \item{opt}{
    TRUE or FALSE: TRUE - structure of triplex will be optimalized, FALSE -
    structure will be drawn without optimalization.
  }
  \item{A.col}{
    Color of Adine base.
  }
  \item{T.col}{
    Color of Thymin base.
  }
  \item{G.col}{
    Color of Guanine base.
  }
  \item{C.col}{
    Color of Cytosine base.
  }
  \item{bgr.col}{
    Color of background.
  }
  \item{bbone.col}{
    Color of backbone.
  }
  \item{bbone.n }{
    Number of sides of backbone bonds.
  }
}

\value{
Instance of \code{\link{DNAStringSet}} object with computed alignment.
}

\author{
Kamil Rajdl, Jiri Hon
}

\seealso{
\code{\link{triplex.diagram}},
\code{\link{triplex.search}},
\code{\link{triplex.alignment}}
}

\examples{
seq <- DNAString("GGAAAGCAATGCCAGGCAGGG")
t <- triplex.search(seq)
triplex.3D(t[1])
}
