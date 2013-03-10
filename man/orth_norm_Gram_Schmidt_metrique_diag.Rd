\name{orth_norm_Gram_Schmidt_metrique_diag}
\alias{orth_norm_Gram_Schmidt_metrique_diag}

\title{ Gram-Schmidt orthonormalization when the metric is diagonal}
\description{this function orthonormalizes a family of vectors in the Gram-Schmidt sense :  the metric is a diagonal matrix stored as a vector  }
\usage{
orth_norm_Gram_Schmidt_metrique_diag(M, Y)
}

\arguments{
  \item{M}{  vector : diagonal of the metric}
  \item{Y}{ list : family of vectors to orthonormalize}
}


\value{ returns a list : orthonormal family of vectors built with the initial family
of vectors (2nd argument: Y) thanks to the orthonormalization
process of Gram-Schmidt. }

\author{
BAR Romain, Université de Lorraine, IECL, Université de Lorraine, INRIA : BIGS group.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ 
\code{\link{orth_Gram_Schmidt_metrique_diag}},
\code{\link{orth_Gram_Schmidt}},
\code{\link{orth_norm_Gram_Schmidt}}
}

\examples{

orth_norm_Gram_Schmidt_metrique_diag(c(3,2),list(c(1,2),c(2,3)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Gram-Schmidt orthogonalization }

