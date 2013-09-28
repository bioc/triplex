\name{triplex-package}
\docType{package}

\alias{triplex-package}
\alias{triplex}

\title{Triplex search and visualization package}

\description{
This package provides functions for the identification and visualization
of potential intramolecular triplex (H-DNA) patterns in DNA sequences.
}

\details{
This package is essentially an R interface to the underlying C
implementation of a dynamic-programming search startegy of the same name
(Lexa et al., 2011). The main functionality of the original program was to detect
the positions of subsequences in a much larger sequence capable of folding
into an intramolecular triplex (H-DNA) made of as many cannonical nucleotide
triplets as possible (see \code{\link{triplex.search}}). In creating its
incarnation in R, we extended this
basic functionality, to include the calculation of exact base-pairing in the
triple helices, which allowed us to extend the functionality of the package
towards visualization showing the exact base-pairing in 1D, 2D or 3D (see
\code{\link{triplex.diagram}} and \code{\link{triplex.3D}}).
}

\author{
Matej Lexa, Tomas Martinek, Kamil Rajdl, Jiri Hon

Maintainer: Jiri Hon <jiri.hon@google.com>
}

\references{
Lexa, M., Martinek, T., Burgetova, I., Kopecek, D., Brazdova, M.: \emph{A
dynamic programming algorithm for identification of triplex-forming
sequences}, In: Bioinformatics, Vol. 27, No. 18, 2011, Oxford, GB, p.
2510-2517, ISSN 1367-4803
}

\seealso{
\code{\link{DNAString}},
\code{\link{triplex.search}},
\code{\link{triplex.diagram}},
\code{\link{triplex.3D}}
}

\examples{
seq <- DNAString("GGAAAGCAATGCCAGGCAGGG")
t <- triplex.search(seq, min_score=10, p_value=1)
triplex.diagram(t[1])
\dontrun{
triplex.3D(t[1])
}
triplex.score.table()
triplex.group.table()
}

\keyword{package}
