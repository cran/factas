\name{PCA_iter}
\alias{PCA_iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs a normed (iterative) PCA on a given data table.
}
\description{
This function permits user to perform fast Principal Components Analysis on (high dimensional) data.  
Using stochastic processes, several estimation may be calculated depending on the user goals:
principal factors, principal axes, eigenvalues and correlation coefficients between the original
variables and the principal factors can be computed.
If \code{stream=TRUE}, \code{PCA_iter} performs an initialisation of the function \code{\link{PCA}} : considering \code{\link[=PCA]{data_init}} observations, the function
\code{PCA_iter} calculates a first estimation of the desired elements and returns them to function \code{PCA}.
If \code{stream=FALSE}, \code{PCA_iter} performs a Principal Components Analysis on the (entire) given matrix and returns the results to function \code{PCA}.
Most of the arguments of \code{PCA} are inherited from the \code{PCA_iter} function.
}
\usage{
PCA_iter(data,stream=TRUE,nb_fact,principal_factors=TRUE,
principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix :
if \code{stream}, \code{PCA_iter} computes a first estimation of the elements of interest on the matrix composed of the aggregation of the
\code{\link[=PCA]{data_init}} first observations; if \code{stream=FALSE}, \code{data} is the entire given data table.
}
  \item{stream}{boolean : stream=TRUE if you consider data streams; stream=FALSE otherwise.
}
 \item{nb_fact}{integer : number of elements to be calculated. For instance, if you need
 the three first principal factors, nb_fact=3.
}
  \item{principal_factors}{ boolean : principal_factors=TRUE if you need the value of some first
principal factors; principal_factors=FALSE otherwise.
}
  \item{principal_axes}{ boolean : principal_axes=TRUE if you need the value of some first
principal axes; principal_axes=FALSE otherwise.
}
  \item{eigenvalues}{boolean : eigenvalues=TRUE if you need the value of some first
eigenvalues; eigenvalues=FALSE otherwise.
}
  \item{corr}{boolean : corr=TRUE if you need the value of correlation coefficient between the
original variables and some first principal factors; corr=FALSE otherwise.
}
  \item{exec_time}{real : execution time. In this function, this parameter is only relevant if stream=FALSE ( if stream=TRUE, see \code{\link{PCA}} ).
}
}

\value{ the value returned to the PCA function is a list: 
 \item{X}{list composed of the \code{\link[=PCA]{nb_fact}} first principal factors}
 \item{L1 }{list composed of the \code{\link[=PCA]{nb_fact}} first eigenvalues}
 \item{A }{list composed of the \code{\link[=PCA]{nb_fact}} first principal axes}
 \item{M1 }{vector : classical metric of the analysis stored as a vector }
 \item{Corr }{list composed of the \code{\link[=PCA]{nb_fact}} first correlation coefficient}
}

\author{BAR Romain, Université de Lorraine, IECL, INRIA: BIGS group.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{PCA}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{principal components analysis}
