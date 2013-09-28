\name{triplex.score.table}
\alias{triplex.score.table}

\title{Get default scoring tables}

\description{
The \code{triplex.score.table} function returns default scoring tables
for parallel and antiparallel triplex types.
}

\usage{
triplex.score.table()
}

\details{

This function is used by \code{\link{triplex.search}} function to get
default scoring tables. These tables correspond exactly with Table 1 published
in (Lexa, 2011) and they represent the strength of Hoogsteen bonds between
thirdstrand nucleotide (row index) and a nucleotide from duplex (column index).

As a common triplet structure is H.WC:WC, to customize scoring tables
just use H as a row index, WC as a column index, then set desired score value and
pass such modified tables through \code{score_table} option of \code{\link{triplex.search}}
interface.

Please keep in mind that for a mismatch (no bond at all) special value -9 is used.
If you want to change mismatch penalization, please use \code{mis_pen} option of
\code{\link{triplex.search}} function.

}

\value{
List of two matrixes, one for parallel triplex types and the other one
for antiparallel.
}

\note{

If you modify the scoring tables,
you should consider changing also default P-value constants (\code{lambda},
\code{mu} and \code{rn}), because these are valid just for the default scoring tables.

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
\code{\link{triplex.group.table}}
\code{\link{triplex.search}},
\code{\link{TriplexViews}},
\code{\link{triplex.diagram}},
\code{\link{triplex.3D}},
\code{\link{triplex.alignment}}
}

\examples{
triplex.score.table()
}
