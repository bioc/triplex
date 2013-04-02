\name{triplex.alignment}
\alias{triplex.alignment}

\title{
Triplex alignment
}

\description{
This function computes best triplex alignment.
}

\details{

Similarly to other DNA multiple sequence alignments the output of the
\code{triplex.alignment} method is stored as \code{\link{DNAStringSet}}
object. This object consists of four sequences: \code{plus} and \code{minus}
sequences representing 5' to 3' and 3' to 5' DNA strand of detected triplex;
one of the \code{anti-plus}, \code{anti-minus}, \code{para-plus} or
\code{para-minus} sequence representing the third triplex strand aligned to
\code{plus} or \code{minus} strand in \code{antiparallel} or \code{parallel}
fashion; and finally \code{loop} sequence representing unpaired loop.  

Please note that all eight triplex types shown in following
figure can be represented using four types of alignments, because each
alignment can correspond to triplex detected either on forward or reverse DNA
strand.

\if{text}{Figure: Triplex types (see HTML or Latex version of this page)}
\if{html}{\figure{types.png}{Figure: Triplex types}}
\if{latex}{\figure{types.pdf}{Figure: Triplex types}}

The input \code{\link{TriplexViews}} object is required to provide additional
algorithm options (see \code{\link{triplex.search}}). These are used for
proper computation of triplex alignment.
}

\usage{
triplex.alignment(triplex)
}

\arguments{
  \item{triplex}{
    \code{\link{TriplexViews}} object including only one triplex.
  }
}

\value{
Instance of \code{\link{DNAStringSet}} object.
}

\author{
Jiri Hon
}

\seealso{
\code{\link{triplex.diagram}},
\code{\link{triplex.3D}},
\code{\link{triplex.search}}
}

\examples{
seq <- DNAString("GGAAAGCAATGCCAGGCAGGG")
t <- triplex.search(seq)
triplex.alignment(t[1])
}
