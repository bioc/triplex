\name{TriplexViews-class}
\docType{class}

% Class:
\alias{class:TriplexViews}
\alias{TriplexViews-class}
\alias{TriplexViews}
\alias{TriplexViews-constructor}

% Methods:
\alias{show,TriplexViews-method}
\alias{score,TriplexViews-method}
\alias{pvalue,TriplexViews-method}
\alias{ins,TriplexViews-method}
\alias{type,TriplexViews-method}
\alias{lstart,TriplexViews-method}
\alias{lwidth,TriplexViews-method}
\alias{lend,TriplexViews-method}
\alias{strand,TriplexViews-method}
\alias{as.character,TriplexViews-method}
\alias{toString,TriplexViews-method}
\alias{[,TriplexViews-method}

\title{The TriplexViews class}

\description{
The TriplexViews class is a container for storing a set of triplexes
identified in the same DNA sequence (an \code{\link{DNAString}} object).
Each triplex is defined by its start/end locations, score, P-value,
insertion number, type, loop start, loop end and strand identification.

TriplexViews object contains also a parameter vector plus score and group tables
that stores custom algorithm options that were used for triplex search.
This is necessary for proper triplex visualization by \code{\link{triplex.diagram}}
and \code{\link{triplex.3D}} functions.
}

\details{
A TriplexViews object is in fact a particular case of an
\link[Biostrings:XStringViews-class]{XStringViews} 
object (the TriplexViews class contains the
\link[Biostrings:XStringViews-class]{XStringViews} class) so it 
can be manipulated in a similar manner. See
\code{?\link[Biostrings:XStringViews-class]{XStringViews}} for 
detailed information.

If you are interested in algorithm options that are stored in
the TriplexViews object, see parameters of \code{\link{triplex.search}}
function. These options are required by visualization functions for
proper computation of triplex alignment.
}

\section{Constructor}{
\describe{
There is no public constructor for TriplexViews object as it
stores search algorithm options. TriplexViews object would not be useful
without algorithm options attached. For more information see description
or details.
}
}

\section{Accessor-like methods}{
All the accessor-like methods defined for
\link[Biostrings:XStringViews-class]{XStringViews} objects
work on TriplexViews objects. In addition, the following accessors are
defined for TriplexViews objects:

	\describe{
	\item{}{
	\code{score(x)}:
	A vector of non-negative integers containing the scores of triplexes.
	}
	\item{}{
	\code{pvalue(x)}:
	A vector of non-negative doubles containing the P-values of triplexes.
	}
	\item{}{
	\code{ins(x)}:
	A vector of non-negative integers containing the number of insertions/deletions
	in triplexes.
	}
	\item{}{
	\code{type(x)}:
	A vector of non-negative integers containing the triplex type.
	}
	\item{}{
	\code{lstart(x)}:
	A vector of non-negative integers containing the triplex loop starts.
	}
	\item{}{
	\code{lwidth(x)}:
	A vector of non-negative integers containing the triplex loop widths.
	}
	\item{}{
	\code{lend(x)}:
	A vector of non-negative integers containing the triplex loop ends.
	}
	\item{}{
	\code{strand(x)}:
	A vector of '+' or '-' signs to identify on which strand the triplex was found.
	}
	\item{}{
	\code{toString(x)}:
	Converts \code{TriplexViews} object into vector of strings.
	}
	\item{}{
	\code{[x]}:
	Index accessor.
	}
	}
}

\note{
The only standard way to create a TriplexViews object
is to use \code{\link{triplex.search}} function.
}

\author{Jiri Hon}

\seealso{
\code{\link{triplex.search}},
\code{\link{triplex.diagram}},
\code{\link{triplex.3D}},
\code{\link{XStringViews}},
\code{\link{triplex.alignment}}
}

\examples{  
seq <- DNAString("GAAGAAGAAGAAGAAGAAGAAGAAGAAGAA")
t <- triplex.search(seq, min_score=10, p_value=1)
start(t)
end(t)
score(t)
pvalue(t)
ins(t)
type(t)

# Search triplex with maximal score
t[score(t) == max(score(t))]

# Sort triplexes by score
t[order(score(t), decreasing=TRUE)]
}

\keyword{methods}
\keyword{classes}
