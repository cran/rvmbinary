\name{rvmkernel}
\alias{rvmkernel}
\alias{rvmkernel.formula}
\alias{rvmkernel.default}
\alias{print.rvmkernel}
\title{Calculate Kernels for RVM Binary}

\description{
  The Relevance Vector Machine is a Bayesian model for regression and
  classification of identical functional form to the support vector
  machine. This method calculates kernels for \code{rvmbinary}
  }
\usage{

\method{rvmkernel}{formula}(formula, data=NULL, ...)

\method{rvmkernel}{default}(x,kernel= "gaus",parameters=c(0.1), ...)
\method{print}{rvmkernel}(x, ...)
}



\arguments{
 \item{formula}{Formula interface for rvmkernel}
\item{data}{For use with formula interface}
  \item{x}{The data to be fit by RVM.
    When not using a formula x can be a matrix or vector containing the training
    data}
	  
 \item{kernel}{the kernel function used in training and predicting. The AA Kernel is supplied built in which is used by setting the kernel parameter to "aa".
   The Kernlab kernels are also supplied if Kernlab is installed which provides the most popular kernel functions. These can be used by setting the kernel parameter to the following
    strings:
    \itemize{
    \item \code{rbfdot} Radial Basis kernel "Gaussian"

      \item \code{polydot} Polynomial kernel

      \item \code{vanilladot} Linear kernel 

      \item \code{tanhdot} Hyperbolic tangent kernel 

      \item \code{laplacedot} Laplacian kernel 

      \item \code{besseldot} Bessel kernel 

      \item \code{anovadot} ANOVA RBF kernel 

      \item \code{splinedot} Spline kernel 

      \item \code{stringdot} String kernel 

    }
	(default = "rbfdot")
  }
  \item{parameters}{a vector of hyper-parameters (kernel parameters).
    This is a vector which contains the parameters to be used with the
    kernel function. For valid parameters for existing kernels are :

    \itemize{
      \item \code{c(sigma)} inverse kernel width for the Radial Basis
      kernel function "rbfdot" and the Laplacian kernel "laplacedot". And the only parameter for the AA kernel (lambda) which can be set between 0.5-1.0.

      \item \code{c(degree, scale, offset)} for the Polynomial kernel "polydot"

      \item \code{c(scale, offset)} for the Hyperbolic tangent kernel
      function "tanhdot"

      \item \code{c(sigma, order, degree)} for the Bessel kernel "besseldot". 

      \item \code{c(sigma, degree)} for the ANOVA kernel "anovadot".
      
      \item \code{c(length, lambda, normalized)} for the "stringdot" kernel
      where length is the length of the strings considered, lambda the
      decay factor and normalized a logical parameter determining if the
      kernel evaluations should be normalized.
    }
    (default = c(0.1).}

\item{...}{optional parameters to be passed to the low level function
    \code{rvmkernel.default}.}

}
\details{Produces the kernel necessary for RVM.
  }
\value{
  An S4 object of class "rvmkernel" containing the calculated kernel.  
  
  \item{BASIS}{The calculated kernel}
 \item{type}{The kernel used}
  \item{parameter}{The parameter for the kernel}
  \item{data}{the training data used}  
}
\references{
Tipping, M. E. and A. C. Faul (2003). Fast marginal likelihood maximisation for sparse Bayesian models. In C. M. Bishop and B. J. Frey (Eds.), Proceedings of the Ninth International Workshop on Artificial Intelligence and Statistics, Key West, FL, Jan 3-6.

Robert Lowe, Hamse Yussuf Mussa, John Mitchell, and Robert Glen. Classifying Molecules Using a Sparse Probabilistic Kernel Binary Classifier. J. Chem. Inf. Model. DOI: 10.1021/ci200128w}
\author{ Robert Lowe \cr
\email{ral64@cam.ac.uk}}

\seealso{ \code{\link{ksvm}}}
\examples{
# create data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(401,sd=0.05)

# Calculate the kernel
tmp=rvmkernel(x,kernel="rbfdot")

}
\keyword{regression}
\keyword{classification}
\keyword{nonlinear}
