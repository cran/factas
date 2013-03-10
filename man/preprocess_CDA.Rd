\name{preprocess_CDA}
\alias{preprocess_CDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
preprocessing for Canonical Discriminant Analysis
}
\description{
performs an encoding of the original data set before launching CDA
}
\usage{
preprocess_CDA(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{matrix : data set before encoding
}
}

\value{
matrix : data set after encoding
}

\author{
Romain Bar, Universite de Lorraine, IECL, INRIA : BIGS group
}
\seealso{ 
\code{\link{preprocess_MCA}},
\code{\link{preprocess_CA}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ canonical discriminant analysis }
\keyword{ encoding }
