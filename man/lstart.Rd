\name{lstart}

\alias{lstart}

\title{Loop start accessor}

\description{
  Gets the loop starts of an object.
}

\usage{
lstart(x, ...)
}

\arguments{
  \item{x}{An object to get the loop start values of.}
  \item{\dots}{Additional arguments.}
}

\examples{
t <- triplex.search(DNAString("TATTTATTTTTTCATCTTCTTTTTTTATTTTT"))
lstart(t);
}

\keyword{methods}
