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


\note{

If you change the scoring tables,
you should consider recalculation of P-value constants (\code{lambda},
\code{mu} and \code{rn}).

}

\details{

The \code{triplex.search} function identifies potential intramolecular
triplex-forming sequences in DNA sequence represented as
a \code{\link{DNAString}} object.

}

\value{
List of two matrixes, one for parallel triplex types and the other one
for antiparallel.
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
\code{\link{triplex.search}},
\code{\link{TriplexViews}},
\code{\link{triplex.diagram}},
\code{\link{triplex.3D}},
\code{\link{triplex.alignment}}
}

\examples{
# GAA triplet repeats involved in Friedreichs's ataxia
seq <- DNAString("GAAGAAGAAGAAGAAGAAGAAGAAGAAGAA")

# Search specific triplex types (see details section)
triplex.search(seq, type=c(2,3))

# Search all triplex types
t <- triplex.search(seq)

# Sort triplexes by score
t[order(score(t), decreasing=TRUE)]
}

\keyword{interface}
