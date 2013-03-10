\name{CDA_iter}
\alias{CDA_iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs a (iterative) CDA on a given data table.
}
\description{
This function permits user to perform fast Canonical Discriminant Analysis on (high dimensional) data.  
Using stochastic processes, several estimation may be calculated depending on the user goals:
principal canonical factors, principal canonical axes, eigenvalues and correlation coefficients between the original
variables and the principal canonical factors can be computed.
If \code{stream=TRUE}, \code{CDA_iter} performs an initialisation of the function \code{\link{CDA}} : considering \code{\link[=CDA]{data_init}} observations, the function
\code{CDA_iter} calculates a first estimation of the desired elements and returns them to function \code{CDA}.
If \code{stream=FALSE}, \code{CDA_iter} performs a Canonical Discriminant Analysis on the (entire) given matrix and returns the results to function \code{CDA}.
Most of the arguments of \code{CDA} are inherited from the \code{CDA_iter} function.

Reminding that CDA is an analysis achieved on two groups of variables only (one quantitative and the other nominative or qualitative), we precise that the \code{CDA_iter} function computes the canonical elements
(principal factors, principal axes, ...) associated with the first group of (quantitative) variables only.

To compute the function, the group of quantitative variables must be the first group in the dataset and the values taken by the nominative variables must be exclusive categories.
Indeed, a pre-process is done on the nominative variables to obtain a codage composed of 1 and 0, essential to achieve the analysis.


}
\usage{
CDA_iter(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,
principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix :
if \code{stream=TRUE}, \code{CDA_iter} computes a first estimation of the elements of interest on the matrix composed of the aggregation of the
\code{\link[=CDA]{data_init}} first observations; if \code{stream=FALSE}, \code{data} is the entire given data table.
}
  \item{groups}{vector ( of dimension 2 ) : the first component (>1) corresponds to the dimension/length of the group of quantitative variables, the
second component (>1) gives the number of categories for the nominative variable.
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
  \item{exec_time}{real : execution time. In this function, this parameter is only relevant if stream=FALSE ( if stream=TRUE, see \code{\link{CDA}} )
}

}

\value{the value returned to the CDA function is a list:


 \item{rbar}{mean vector of the first group}
 \item{sbar}{mean vector of the second group}
 \item{A}{list composed of the \code{\link[=CDA]{nb_fact}} first principal canonical axes}
 \item{N_calc}{covariance matrix of the first group}
 \item{M_calc}{inverse of the covariance matrix of the first group}
 \item{J_calc, K_calc, A2_calc, B_calc, Cov_SR_calc, Cov_S_calc}{ useful matrices to define an estimation processus}
 \item{D_calc}{vector usefull to define an estimation processus}
 \item{L1}{list composed of the \code{\link[=CDA]{nb_fact}} first eigenvalues}
 \item{X}{list composed of the \code{\link[=CDA]{nb_fact}} first principal canonical factors}
 \item{Corr}{list composed of the \code{\link[=CDA]{nb_fact}} first correlation coefficient}
}
\references{ ACP projetée de données séquentielles, Monnez J-M., 42èmes Journées de Statistique (2010) 
}
\author{
Bar Romain, Université de Lorraine, IECL, INRIA : BIGS group
}
\seealso{ 
\code{\link{CDA}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{canonical discriminant analysis}
