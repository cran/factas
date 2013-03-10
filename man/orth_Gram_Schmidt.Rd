\name{orth_Gram_Schmidt}
\alias{orth_Gram_Schmidt}

\title{ Gram-Schmidt orthogonalization}
\description{this function orthogonalizes a family of vectors in the Gram-Schmidt sense : the metric is a matrix  }
\usage{
orth_Gram_Schmidt(M, Y)
}

\arguments{
  \item{M}{ matrix : the metric}
  \item{Y}{ list : family of vectors to orthogonalize}
}



\value{ returns a list : orthogonal family of vectors built with the initial family
of vectors (2nd argument : Y) thanks to the orthogonalization
process of Gram-Schmidt }

\author{
Romain Bar, Université de Lorraine, IECL, INRIA: BIGS group
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 

\code{\link{orth_Gram_Schmidt_metrique_diag}}

\code{\link{orth_norm_Gram_Schmidt}}

\code{\link{orth_norm_Gram_Schmidt_metrique_diag}}

}

\examples{
orth_Gram_Schmidt(diag(c(3,2)),list(c(1,2),c(2,3)))

stopifnot( orth_Gram_Schmidt(diag( c(1,2)),list( c(1,2),c(7,5) ) ) [[2]]== c(4,-1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Gram-Schmidt orthogonalization }

