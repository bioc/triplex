\name{triplex.group.table}
\alias{triplex.group.table}

\title{Get default isogroup tables}

\description{
The \code{triplex.group.table} function returns default isogroup tables
for parallel and antiparallel triplex types.
}

\usage{
triplex.group.table()
}

\details{

This function is used by \code{\link{triplex.search}} function to get
default isogroup tables. These tables correspond exactly with Table 1 published
in (Lexa, 2011) and they represent the isomorphic groups of triplets.

As a common triplet structure is H.WC:WC, to customize isogroup tables
just use H as a row index, WC as a column index, then set desired group number and
pass such modified tables through \code{group_table} option of \code{\link{triplex.search}}
interface.

}

\value{
List of two matrixes, one for parallel triplex types and the other one
for antiparallel.
}

\note{

If you modify the isogroup tables,
you should consider changing also default P-value constants (\code{lambda},
\code{mu} and \code{rn}), because these are valid just for the default isogroup tables.

}

\references{
Lexa, M., Martinek, T., Burgetova, I., Kopecek, D., Brazdova, M.: \emph{A
dynamic programming algorithm for identification of triplex-forming
sequences}, In: Bioinformatics, Vol. 27, No. 18, 2011, Oxford, GB, p.
2510-2517, ISSN 1367-4803
}

\author{
Tomas Martinek, Jiri Hon
}

\seealso{
\code{\link{triplex.score.table}}
\code{\link{triplex.search}},
\code{\link{TriplexViews}},
\code{\link{triplex.diagram}},
\code{\link{triplex.3D}},
\code{\link{triplex.alignment}}
}

\examples{
triplex.group.table()
}
