\name{CA_iter}
\alias{CA_iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs a (iterative) Correpondence Analysis on a given data table. 
}
\description{
This function permits user to perform fast  Correspondence Analysis on (high dimensional) data.  
Using stochastic processes, several estimation may be calculated depending on the user goals:
principal canonical factors, principal canonical axes, eigenvalues and correlation coefficients between the original
variables and the principal canonical factors can be computed.
If \code{stream=TRUE}, \code{CA_iter} performs an initialisation of the function \code{\link{CA}} : considering \code{\link[=CA]{data_init}} observations, the function
\code{CA_iter} calculates a first estimation of the desired elements and returns them to the \code{CA} function.
If \code{stream=FALSE}, \code{CA_iter} performs a Correspondence Analysis on the (entire) given matrix and returns the results to the \code{CA} function.
Most of the arguments of \code{CA} are inherited from the \code{CA_iter} function.

Reminding that CA is an analysis achieved on two groups of nominative variables only, we precise that the \code{CA_iter} function computes the canonical elements
(principal factors, principal axes, ...) associated with the first group of nominative variables only. If the user needs the canonical elements associated with
the second group only, he may permute the two groups of variables before launching the analysis on that new dataset. If the user needs both canonical elements
associated with the two groups, an option is to launch the \code{\link{MCA_iter}} function (generalization of the CA_iter function
to more than two groups) on the dataset: indeed, \code{MCA_iter} returns the general elements (principal factors, principal axes, ...) that are strongly linked to the canonical elements: for instance, if the two groups are respectively
composed of p and q variables, the principal canonical factors associated with the first (resp. second) group returned by the CA_iter function are of length p (resp. q) whereas 
the principal general factors returned by the MCA_iter function are of length p+q. The principal canonical factors associated with the first (resp. second) group are 
colinears to the vectors composed of the p first (resp. q last) components of the principal general factors.
}
\usage{
CA_iter(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,
principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix :
if \code{stream=TRUE}, \code{CA_iter} computes a first estimation of the elements of interest on the matrix composed of the aggregation of the
\code{\link[=CA]{data_init}} first observations. 
if \code{stream=FALSE}, \code{data} is the entire given data table.
}
  \item{groups}{vector ( of dimension 2) : each component corresponds to the number (>1) of categories for each nominative variable.
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
  \item{exec_time}{real : execution time. In this function, this parameter is only relevant if stream=FALSE ( if stream=TRUE, see \code{\link{CA}} )
}

}


\value{the value returned to the CA function is a list:

 \item{N}{inverse of the covariance matrix of the first group}
 \item{rbar}{mean vector of the first group}
 \item{sbar}{mean vector of the second group}
 \item{A}{list composed of the \code{\link[=CCA]{nb_fact}} first principal canonical axes}
 \item{F, F2, G, B, E_SR}{ useful matrices to define an estimation process }
 \item{E_SS}{(diagonal) covariance matrix of the second group stored as a vector}
 \item{L1}{list composed of the \code{\link[=CCA]{nb_fact}} first eigenvalues}
 \item{X}{list composed of the \code{\link[=CCA]{nb_fact}} first principal canonical factors}
 \item{Corr}{list composed of the \code{\link[=CCA]{nb_fact}} first correlation coefficient}
}
\references{
ACP projetée de données séquentielles, Monnez J-M., 42èmes Journées de Statistique (2010) 
}
\author{
Bar Romain, Université de Lorraine, IECL, INRIA : BIGS group
}
\seealso{ 
\code{\link{CA}},
\code{\link{MCA}},
\code{\link{MCA_iter}}
}

\keyword{correspondence analysis}
