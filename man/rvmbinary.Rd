\name{rvmbinary}
\alias{rvmbinary}
\alias{rvmbinary.formula}
\alias{rvmbinary.default}
\alias{predict.rvmbinary}
\alias{print.rvmbinary}

\title{Relevance Vector Machine}

\description{
  The Relevance Vector Machine is a Bayesian model for regression and
  classification of identical functional form to the support vector
  machine.
  }
\usage{
\method{rvmbinary}{formula}(formula, data=NULL, ...)

\method{rvmbinary}{default}(x,y, kernel="rbfdot", parameters=c(0.1),iterations=100, noisevar=0.1, minmaxdiff = 1e-3, ...) 

\method{print}{rvmbinary}(x, ...)
\method{predict}{rvmbinary}(object,newdata, ...)

}



\arguments{
  \item{formula}{Formula interface for rvmbinary}
  \item{x}{The data to be fit by RVM.
    When not using a formula x can be a matrix or vector containing the training
    data or a kernel matrix of class \code{rvmkernel} of the training data.}
  \item{data}{a data frame containing the variables in the model when using the formula function}
	  
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
    a factor (for classification tasks) or a numeric vector (for
    regression).}
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


  \item{noisevar}{the initial noise variance}
  \item{iterations}{Number of iterations allowed (default: 100)} 
  \item{minmaxdiff}{termination criteria. Stop when max difference is
    equal to this parameter (default:1e-3) }
 \item{object}{The object returned from rvmbinary in order to predict}
\item{newdata}{New data for testing}
\item{...}{optional parameters to be passed to the low level function
    \code{rvmbinary.default}.}
}
\details{The Relevance Vector Machine typically leads to sparser models
  then the SVM. It is probabilistic by nature and any kernel can be used unlike SVM.
  }
\value{
  An S4 object of class "rvmbinary" containing the fitted model.
  Accessor functions can be used to access the slots of the
  object which include :
  
  
  \item{kernel}{The kernel used to produce the model}
  \item{kernelparameter}{The parameter for the kernel}
  \item{nRV}{Number of relevance vectors}
  \item{type}{Classification or Regression}
  \item{vectors}{The relevance vectors}
  \item{values}{The weights for the relevance vectors}

  
  ...
}
\references{
Tipping, M. E. and A. C. Faul (2003). Fast marginal likelihood maximisation for sparse Bayesian models. In C. M. Bishop and B. J. Frey (Eds.), Proceedings of the Ninth International Workshop on Artificial Intelligence and Statistics, Key West, FL, Jan 3-6.

Robert Lowe, Hamse Yussuf Mussa, John Mitchell, and Robert Glen. Classifying Molecules Using a Sparse Probabilistic Kernel Binary Classifier. J. Chem. Inf. Model. DOI: 10.1021/ci200128w}

\author{ Robert Lowe \cr
\email{ral64@cam.ac.uk}}

\seealso{ \code{\link{ksvm}}}
\examples{
#Classification
data(iris)
#Create training and test splits
datatest=rbind(iris[41:50,],iris[91:100,])
datatrain=rbind(iris[1:40,],iris[51:90,])

#Run model
rvm=rvmbinary(datatrain[,-dim(datatrain)[2]],datatrain[,dim(datatrain)[2]],kernel="rbfdot",parameters=0.1,1000)

#Calculate class probability for test set
y=predict(rvm,datatest[,-dim(datatest)[2]])




# create data For REGRESSION
x <- c(seq(-20,-0.1,0.1),seq(0.1,20,0.1))
y <- sin(x)/x + rnorm(400,sd=0.05)

# train relevance vector machine
foo <- rvmbinary(x, y)
foo
# print relevance vectors

# predict and plot
ytest <- predict(foo, x)
plot(x, y, type ="l")
lines(x, ytest, col="red")
}
\keyword{regression}
\keyword{classification}
\keyword{nonlinear}
