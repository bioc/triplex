\name{ins}

\alias{ins}

\title{Insertions accessor}

\description{
  Gets the insertion vector of an object.
}

\usage{
ins(x, ...)
}

\arguments{
  \item{x}{An object to get the insertion values of.}
  \item{\dots}{Additional arguments.}
}

\examples{
t <- triplex.search(DNAString("TATTTATTTTTTCATCTTCTTTTTTTATTTTT"), max_len=11)
ins(t);
}

\keyword{methods}
