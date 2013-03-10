\name{GCCA_iter}
\alias{GCCA_iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs a (iterative) GCCA on a given data table.
}
\description{

This function permits user to perform fast Generalized Canonical Correlation Analysis on (high dimensional) data.  
Using stochastic processes, several estimation may be calculated depending on the user goals:
principal general factors, principal general axes, eigenvalues and correlation coefficients between the original
variables and the principal general factors can be computed.
If \code{stream=TRUE}, \code{GCCA_iter} performs an initialisation of the function \code{\link{GCCA}} : considering \code{\link[=GCCA]{data_init}} observations, the function
\code{GCCA_iter} calculates a first estimation of the desired elements and returns them to function \code{GCCA}.
If \code{stream=FALSE}, \code{GCCA_iter} performs a Generalized Canonical Correlation Analysis on the (entire) given matrix and returns the results to function \code{GCCA}.
Most of the arguments of \code{GCCA} are inherited from the \code{GCCA_iter} function.


Reminding that GCCA is an analysis achieved on q groups of variables (q > or = 2), we precise that the \code{GCCA_iter} function computes the general elements
(principal factors, principal axes, ...) . If the user needs the canonical elements associated with the $k^{th}$ group (k=1,\dots,q), it's possible to deduce them from the
 general elements: the principal canonical factors associated with the $k^{th}$ group are 
 colinears to the vectors composed of the relevant components of the principal general factors.
}
\usage{
GCCA_iter(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,
principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix :
if \code{stream=TRUE}, \code{GCCA_iter} computes a first estimation of the elements of interest on the matrix composed of the aggregation of the
\code{\link[=GCCA]{data_init}} first observations; if \code{stream=FALSE}, \code{data} is the entire given data table.
}
  \item{groups}{vector :
each component of \code{groups} gives the length of the associated group of variables.
}
  \item{stream}{boolean : stream=TRUE if you consider data streams; stream=FALSE otherwise.
}
 \item{nb_fact}{integer : number of elements to be calculated. For instance, if you need
 the three first principal general factors, nb_fact=3
}
  \item{principal_factors}{ boolean : principal_factors=TRUE if you need the value of some first
principal general factors; principal_factors=FALSE otherwise.
}
  \item{principal_axes}{ boolean : principal_axes=TRUE if you need the value of some first
principal general axes; principal_axes=FALSE otherwise.
}
  \item{eigenvalues}{boolean : eigenvalues=TRUE if you need the value of some first
eigenvalues; eigenvalues=FALSE otherwise.
}
  \item{corr}{boolean : corr=TRUE if you need the value of correlation coefficient between the
original variables and some first principal general factors; corr=FALSE otherwise.
}
  \item{exec_time}{real : execution time. In this function, this parameter is only relevant if stream=FALSE ( if stream=TRUE, see \code{\link{GCCA}} ).
}

}

\value{the value returned to the GCCA function is a list:

 \item{A}{list composed of the \code{\link[=GCCA]{nb_fact}} first principal canonical axes}
 \item{C}{covariance matrix of the entire dataset}
 \item{M}{list composed of the inverses of the covariance matrices of the $k^{th}$ group (k=1,\dots,d)}
 \item{M1}{matrix : metric associated to the GCCA : blockwise diagonal matrix composed of the elements of M}
 \item{L1}{list composed of the \code{\link[=GCCA]{nb_fact}} first eigenvalues}
 \item{X}{list composed of the \code{\link[=GCCA]{nb_fact}} first principal canonical factors}
 \item{Corr}{list composed of the \code{\link[=GCCA]{nb_fact}} first correlation coefficient}
}
\references{
Stochastic approximation of the factors of a generalized canonical correlation analysis, Monnez J-M., Statistics and Probability Letters 78, 14 (2008) 2210-2216
}
\author{
Bar Romain, Université de Lorraine, IECL, INRIA : BIGS group.
}
\seealso{ 
\code{\link{GCCA_iter}},
\code{\link{CCA}},
\code{\link{CCA_iter}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{generalized canonical correlation analysis}

