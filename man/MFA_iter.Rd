\name{MFA_iter}
\alias{MFA_iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs a (iterative) MFA on a given data table.
}
\description{
This function permits user to perform fast Multiple Factorial Analysis on (high dimensional) data.  
Using stochastic processes, several estimation may be calculated depending on the user goals:
principal factors, principal axes, eigenvalues and correlation coefficients between the original
variables and the principal factors can be computed.
If \code{stream=TRUE}, \code{MFA_iter} performs an initialisation of the function \code{\link{MFA}} : considering \code{\link[=MFA]{data_init}} observations, the function
\code{MFA_iter} calculates a first estimation of the desired elements and returns them to function \code{MFA}.
If \code{stream=FALSE}, \code{MFA_iter} performs a Multiple Factorial Analysis on the (entire) given matrix and returns the results to function \code{MFA}.
Most of the arguments of \code{MFA} are inherited from the \code{MFA_iter} function.
}
\usage{
MFA_iter(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,
principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix :
if \code{stream=TRUE}, \code{MFA_iter} computes a first estimation of the elements of interest on the matrix composed of the aggregation of the
\code{\link[=MFA]{data_init}} first observations; if \code{stream=FALSE}, \code{data} is the entire given data table.
}
  \item{groups}{vector :
each component of \code{groups} gives the dimension/length of the associated group.
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
  \item{exec_time}{real : execution time. In this function, this parameter is only relevant if stream=FALSE ( if stream=TRUE, see \code{\link{MFA}} ).
}

}

\value{ the value returned to the MFA function is a list: 
 \item{X}{list composed of the \code{\link[=MFA]{nb_fact}} first principal factors}
 \item{Xk}{list composed of the first principal factors of the PCA on the $k^{th}$ group, useful to calculate the metric}
 \item{Ck}{list composed of the covariance matrices associated with the $k^{th}$ group, useful to calculate the metric }
 \item{Mk }{list composed of the classical metric of the PCA on the $k^{th}$ group (stored as a vector), usefull to calculate the MFA metric }
 \item{L1 }{list composed of the \code{\link[=MFA]{nb_fact}} first eigenvalues}
 \item{A }{list composed of the \code{\link[=MFA]{nb_fact}} first principal axes}
 \item{Corr }{list composed of the \code{\link[=MFA]{nb_fact}} first correlation coefficient}
}
\references{
Multiple Factor Analysis : Main Features and Application to Sensory Data, Jérôme Pagès, 2004.
Approximation stochastique en analyse factorielle multiple, Jean-Marie Monnez, Publications de l'Institut de Statistique de l'Université de Paris L, 3 (2006) 27-45
}
\author{
BAR Romain, IECL, Université de Lorraine, INRIA : BIGS group.
}
\seealso{ 
\code{\link{MFA}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multiple factor analysis}
