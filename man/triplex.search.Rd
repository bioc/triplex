\name{triplex.search}
\alias{triplex.search}

\title{Search intramolecular triplex-forming sequences in DNA}

\description{
The \code{triplex.search} function identifies potential intramolecular
triplex-forming sequences in DNA.
}

\usage{
triplex.search(
  dna, 
  type        = 0:7,
  min_score   = 15,
  p_value     = 0.05,
  min_len     = 6,
  max_len     = 25,
  min_loop    = 3,
  max_loop    = 10,
  seq_type    = 'eukaryotic',
  score_table = 'default',
  group_table = 'default',
  lambda_par  = 'default',
  lambda_apar = 'default',
  mu_par      = 'default',
  mu_apar     = 'default',
  rn_par      = 'default',
  rn_apar     = 'default',
  dtwist_pen  = 'default',
  ins_pen     = 'default',
  iso_pen     = 'default',
  iso_bonus   = 'default',
  mis_pen     = 'default')
}

\arguments{
  \item{dna}{
    A \code{\link{DNAString}} object.
  }
  \item{type}{
    Vector of triplex types (0..7) to be searched for.
  }
  \item{min_score}{
    Minimal score treshold.
  }
  \item{p_value}{
    Acceptable P-value.
  }
  \item{min_len}{
    Minimal triplex length.
  }
  \item{max_len}{
    Maximal triplex length. 
  }
  \item{min_loop}{
    Minimal triplex loop length. Can not be lower than one.
  }
  \item{max_loop}{
    Maximal triplex loop length.
  }
  \item{seq_type}{
    Type of input sequence. Possible options: prokaryotic, eukaryotic.
  }
  \item{score_table}{
    Scoring table for parallel and antiparallel triplex types. Default
    is the same as \code{\link{triplex.score.table}} output.
    Before changing this option, please read \code{\link{triplex.score.table}}
    help carefully.
  }
  \item{group_table}{
    Isomorphic group table for parallel and antiparallel triplex types. Default
    is the same as \code{\link{triplex.group.table}} output.
    Before changing this option, please read \code{\link{triplex.group.table}}
    help carefully.
  }
  \item{lambda_par}{
    Lambda for parallel triplex types 0,1,2,3. Default for prokaryotic
    sequence is 0.8892, for eukaryotic 0.8433.
  }
  \item{lambda_apar}{
    Lambda for antiparallel triplex types 4,5,6,7. Default for prokaryotic
    sequence is 0.8092, for eukaryotic 0.6910.
  }
  \item{mu_par}{
    Mu for parallel triplex types 0,1,2,3. Default for prokaryotic
    sequence is 7.4805, for eukaryotic 7.5835.
  }
  \item{mu_apar}{
    Mu for antiparallel triplex types 4,5,6,7. Default for prokaryotic
    sequence is 7.6569, for eukaryotic 7.9611.
  }
  \item{rn_par}{
    Hit ratio (reported hits to sequence length) for parallel triplex
    types 0,1,2,3.
    
    Default for prokaryotic sequence is 0.0406, for eukaryotic 0.0304.
  }
  \item{rn_apar}{
    Hit ratio (reported hits to sequence length) for antiparallel triplex
    types 4,5,6,7.
    
    Default for prokaryotic sequence is 0.0273, for eukaryotic 0.0405.
  }
  \item{dtwist_pen}{
    Dtwist penalization, default is 7.
  }
  \item{ins_pen}{
    Insertion penalization, default is 9.
  }
  \item{iso_pen}{
    Isomorphic group change penalization, default is 5.
  }
  \item{iso_bonus}{
    Isomorphic group stay bonus, default is 0.
  }
  \item{mis_pen}{
    Mismatch penalization, default is 7.
  }
}


\note{

If you modify the penalization options (\code{dtwist_pen}, \code{ins_pen},
\code{iso_pen}, \code{iso_bonus}, \code{mis_pen}), scoring tables (\code{score_table})
or isogroup tables (\code{group_table}),
you should consider changing also default P-value constants (\code{lambda},
\code{mu} and \code{rn}) to get relevant P-values.

}

\details{

The \code{triplex.search} function identifies potential intramolecular
triplex-forming sequences in DNA sequence represented as
a \code{\link{DNAString}} object.

Based on triplex position (forward or reverse strand) and its third strand
orientation, up to 8 types of triplexes are distinguished by the function (see
the following figure). By default, the function detects all 8 types, however 
this behavior can be changed by setting the \code{type} parameter to any value
or a subset of values in the range 0 to 7.

\if{text}{Figure 1: Triplex types (see HTML or Latex version of this page)}
\if{html}{\figure{types.png}{Figure 1: Triplex types}}
\if{latex}{\figure{types.pdf}{Figure 1: Triplex types}}

Detected triplexes are returned as instances of the 
\code{\link{TriplexViews}} class,
which represents the basic container for storing a set of views on the same
input sequence similarly to the \code{\link{XStringViews}} object (in fact
\code{\link{TriplexViews}} only extends the \code{\link{XStringViews}} class
with a number of displayed columns). Each triplex view is defined by start
and end locations, width, score, P-value, number of insertions, type, strand, loop
start and loop end.  Please note, that the strand orientation depends on
triplex type only. The \code{triplex.search} function assumes that the input
DNA sequence represents the forward strand.

Basic requirements for the shape or length of detected triplexes can be
defined using four parameters: \code{min_len}, \code{max_len},
\code{min_loop} and \code{max_loop}. While \code{min_len} and \code{max_len}
specify the length of the triplex stem composed of individual triplets,
\code{min_loop} and \code{max_loop} parameters define the range of lengths
for the unpaired loop at the top of the triplex. A graphical representation of
these parameters is shown in the following figure.  Please note, these 
parameters also impact the overall computation time. For longer triplexes, 
larger space has to be explored and thus more computation time is consumed. 

\if{text}{Figure 2: Triplex scheme (see HTML or Latex version of this page)}
\if{html}{\figure{triplex.png}{Figure 2: Triplex scheme}}
\if{latex}{\figure{triplex.pdf}{Figure 2: Triplex scheme}}

The quality of each triplex is represented by its score value. A higher score
value represents a higher-quality triplex. This quality is decreased by several
types of imperfections at the level of triplets, such as character mismatch,
insertion, deletion, isomorphic group change etc. Penalization constants for
these imperfections can be setup using the following parameters: \code{mis_pen},
\code{ins_pen}, \code{iso_pen}, \code{iso_bonus} and
\code{dtwist_pen}. Detailed information about the scoring function and
penalization parameters can be found in (Lexa et al., 2011). It is highly
recommended to see (Lexa et al., 2011) prior to changing any penalization 
parameters.

The \code{triplex.search} function can output a large list containing tens of
thousands of potential triplexes. The size of these results can be reduced
using two filtration mechanisms: (1) by specifying the minimal acceptable
score value using the \code{min_score} parameter or (2) by specifying the
maximum acceptable P-value of results using the \code{p_value} parameter. The
P-value represents the probability of occurrence of detected triplexes in
random sequence. By default, only triplexes with P-value equal or less than
0.05 are reported. Calculation of P-value depends on two extreme value
distribution parameters \code{lambda} and \code{mi}.  By default, these
parameters are set up for searching in human genome sequences. It is highly
recommended to see (Lexa et al., 2011) prior to changing either of the
\code{lambda} and \code{mi} parameters.

}

\value{
Instance of \code{\link{TriplexViews}} object based on
\code{\link{XStringViews}} class.
}

\references{
Lexa, M., Martinek, T., Burgetova, I., Kopecek, D., Brazdova, M.: \emph{A
dynamic programming algorithm for identification of triplex-forming
sequences}, In: Bioinformatics, Vol. 27, No. 18, 2011, Oxford, GB, p.
2510-2517, ISSN 1367-4803
}

\author{
Matej Lexa, Tomas Martinek, Jiri Hon
}

\seealso{
\code{\link{TriplexViews}},
\code{\link{triplex.score.table}}
\code{\link{triplex.group.table}}
\code{\link{triplex.diagram}},
\code{\link{triplex.3D}},
\code{\link{triplex.alignment}}
}

\examples{
# GAA triplet repeats involved in Friedreichs's ataxia
seq <- DNAString("GAAGAAGAAGAAGAAGAAGAAGAAGAAGAA")

# Search specific triplex types (see details section)
triplex.search(seq, type=c(2,3), min_score=10, p_value=1)

# Search all triplex types
t <- triplex.search(seq, min_score=10, p_value=1)

# Sort triplexes by score
t[order(score(t), decreasing=TRUE)]
}

\keyword{interface}
