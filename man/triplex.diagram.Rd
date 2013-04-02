\name{triplex.diagram}
\alias{triplex.diagram}

\title{
Triplex visualization, diagram representation 
}

\description{
This function visualizes a \code{\link{TriplexViews}} object as a 2D diagram.
Nucleotides are drawn as characters in circles and bonds as lines between
 them (Watson-Crick or Hogsteen).
}

\details{
The input \code{\link{TriplexViews}} object is required to provide additional
algorithm options (see \code{\link{triplex.search}}). These are used for
proper computation of triplex alignment.

An example of a graphical output corresponding to a triplex of type 3 with 
DNA sequence "GGAAAGCAATGCCAGGCAGGG" is shown in the following figure

\if{text}{Figure: Example 2D visualization (see HTML or Latex version of this page)}
\figure{triplex2d.png}
}

\usage{
triplex.diagram(triplex, circles = TRUE, mbonds.lty = 1, 
                mbonds.lwd = 2.5, wcbonds.lty = 1, wcbonds.lwd = 1, 
                hbonds.lty = 2, hbonds.lwd = 1, labels.cex = 1, circles.cex = 1,
                margin = 0.1, bonds.length = 0.07)
}

\arguments{
  \item{triplex}{
    \code{\link{TriplexViews}} object including only one triplex.
  }
  \item{circles}{
    TRUE or FALSE: TRUE - nucleotides are drawn as characters in circles, FALSE -
    nucleotides are drawn just as characters.
  }
  \item{mbonds.lty}{
    Type of main (skelet) bonds lines.
  }
  \item{mbonds.lwd}{
    Width of main (skelet) bonds lines.
  }
  \item{wcbonds.lty}{
    Type of Watson-Crick bonds lines.
  }
  \item{wcbonds.lwd}{
    Width of Watson-Crick bonds lines.
  }
  \item{hbonds.lty}{
    Type of Hoogsteen bonds lines.
  }
  \item{hbonds.lwd}{
    Width of Hoogsteen bonds lines.
  }
  \item{labels.cex}{
    Multiplier of size of labels of nucleotides.
  }
  \item{circles.cex}{
    Multiplier of size of nucleotides.
  }
  \item{margin}{
    Left and right margin of the picture.
  }
  \item{bonds.length}{
    Length of lines representing Watson-Crick and Hoogsteen bonds.
  }
}

\value{
Instance of \code{\link{DNAStringSet}} object with computed alignment.
}

\author{
Kamil Rajdl, Jiri Hon
}

\seealso{
\code{\link{triplex.3D}},
\code{\link{triplex.search}},
\code{\link{triplex.alignment}}
}

\examples{
seq <- DNAString("GGAAAGCAATGCCAGGCAGGG")
t <- triplex.search(seq)
triplex.diagram(t[1])
}
