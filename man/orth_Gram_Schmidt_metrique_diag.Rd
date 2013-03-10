\name{orth_Gram_Schmidt_metrique_diag}
\alias{orth_Gram_Schmidt_metrique_diag}

\title{ Gram-Schmidt orthogonalization when the metric is diagonal}
\description{this function orthogonalizes a family of vectors in the Gram-Schmidt sense : the metric is a diagonal matrix stored as a vector }
\usage{
orth_Gram_Schmidt_metrique_diag(M, Y)
}

\arguments{
  \item{M}{ vector : diagonal of the metric}
  \item{Y}{ list : family of vectors to orthogonalize}
}


\value{ returns a list : orthogonal family of vectors built with the initial family
of vectors (2nd argument : Y) thanks to the orthogonalization
process of Gram-Schmidt. }


\author{
BAR Romain, Université de Lorraine, IECL, Université de Lorraine, INRIA : BIGS group.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{orth_Gram_Schmidt}},
\code{\link{orth_norm_Gram_Schmidt}},
\code{\link{orth_norm_Gram_Schmidt_metrique_diag}}
}

\examples{

orth_Gram_Schmidt_metrique_diag(c(3,2),list(c(1,2),c(2,3)))

stopifnot( orth_Gram_Schmidt_metrique_diag( c(1,2),list( c(1,2),c(7,5) ) ) [[2]]== c(4,-1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Gram-Schmidt orthogonalization }

