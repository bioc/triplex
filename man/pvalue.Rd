\name{pvalue}

\alias{pvalue}

\title{P-value accessor}

\description{
  Gets the P-values of an object.
}

\usage{
pvalue(x, ...)
}

\arguments{
  \item{x}{An object to get the P-values of.}
  \item{\dots}{Additional arguments.}
}

\examples{
t <- triplex.search(DNAString("TATTTATTTTTTCATCTTCTTTTTTTATTTTT"))
pvalue(t);
}

\keyword{methods}
