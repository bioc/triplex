\name{lend}

\alias{lend}

\title{Loop end accessor}

\description{
  Gets the loop ends of an object.
}

\usage{
lend(x, ...)
}

\arguments{
  \item{x}{An object to get the loop end values of.}
  \item{\dots}{Additional arguments.}
}

\examples{
t <- triplex.search(DNAString("TATTTATTTTTTCATCTTCTTTTTTTATTTTT"))
lend(t);
}

\keyword{methods}
