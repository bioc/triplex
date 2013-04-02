\name{lwidth}

\alias{lwidth}

\title{Loop width accessor}

\description{
  Gets the loop widths of an object.
}

\usage{
lwidth(x, ...)
}

\arguments{
  \item{x}{An object to get the loop width values of.}
  \item{\dots}{Additional arguments.}
}

\examples{
t <- triplex.search(DNAString("TATTTATTTTTTCATCTTCTTTTTTTATTTTT"))
lwidth(t);
}

\keyword{methods}
