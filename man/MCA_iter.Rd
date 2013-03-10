\name{MCA_iter}
\alias{MCA_iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs a (iterative) MCA on a given data table.
}
\description{
This function permits user to perform fast Multiple Correspondence Analysis on (high dimensional) data.  
Using stochastic processes, several estimation may be calculated depending on the user goals:
principal canonical factors, principal canonical axes, eigenvalues and correlation coefficients between the original
variables and the principal canonical factors can be computed.
If \code{stream=TRUE}, \code{MCA_iter} performs an initialisation of the function \code{\link{MCA}} : considering \code{\link[=MCA]{data_init}} observations, the function
\code{MCA_iter} calculates a first estimation of the desired elements and returns them to function \code{MCA}.
If \code{stream=FALSE}, \code{MCA_iter} performs a Multiple Correspondence Analysis on the (entire) given matrix and returns the results to function \code{MCA}.
Most of the arguments of \code{MCA} are inherited from the \code{MCA_iter} function.

}
\usage{
MCA_iter(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,
principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix :
if \code{stream=TRUE}, \code{MCA_iter} computes a first estimation of the elements of interest on the matrix composed of the aggregation of the
\code{\link[=MCA]{data_init}} first observations; if \code{stream=FALSE}, \code{data} is the entire given data table.
}
  \item{groups}{vector :
each component of \code{groups} gives the number (>1) of categories for each nominative variable.
}
  \item{stream}{boolean : stream=TRUE if you consider data streams; stream=FALSE otherwise.
}
 \item{nb_fact}{integer : number of elements to be calculated. For instance, if you need
 the three first principal canonical factors, nb_fact=3.
}
  \item{principal_factors}{ boolean : principal_factors=TRUE if you need the value of some first
principal canonical factors; principal_factors=FALSE otherwise.
}
  \item{principal_axes}{ boolean : principal_axes=TRUE if you need the value of some first
principal canonical axes; principal_axes=FALSE otherwise.
}
  \item{eigenvalues}{boolean : eigenvalues=TRUE if you need the value of some first
eigenvalues; eigenvalues=FALSE otherwise.
}
  \item{corr}{boolean : corr=TRUE if you need the value of correlation coefficient between the
original variables and some first principal canonical factors; corr=FALSE otherwise.
}
  \item{exec_time}{real : execution time. In this function, this parameter is only relevant if stream=FALSE ( if stream=TRUE, see \code{\link{CA}} ).
}

}

\value{the value returned to the MCA function is a list:
 \item{rbar}{ mean vector of the data set }
 \item{A}{list composed of the \code{\link[=MCA]{nb_fact}} first principal canonical axes}
 \item{C}{covariance matrix of the entire dataset}
 \item{M}{list composed of the inverses of the covariance matrices of the $k^{th}$ group (k=1,\dots,d)}
 \item{M1}{metric associated to the MCA: blockwise diagonal matrix composed of the elements of M}
 \item{L1}{list composed of the \code{\link[=MCA]{nb_fact}} first eigenvalues}
 \item{X}{list composed of the \code{\link[=MCA]{nb_fact}} first principal canonical factors}
 \item{Corr}{list composed of the \code{\link[=MCA]{nb_fact}} first correlation coefficient}
}
\references{
http://en.wikipedia.org/wiki/Multiple_correspondence_analysis
}
\author{
Bar Romain, Université de Lorraine, IECL, INRIA : BIGS group.
}
\seealso{ 
\code{\link{MCA}},
\code{\link{CA}},
\code{\link{CA_iter}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multiple correspondence analysis}
