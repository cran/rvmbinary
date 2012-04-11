\name{gaus}
\alias{gaus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create the Gaussian kernel for RVM
}
\description{
  The \code{gaus} function applies the Gaussian kernel to the data}
\usage{
gaus(data, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    The data with the class label as the last column
}
  \item{lambda}{
    The lambda parameter for the Gaussian kernel
}
}

\value{
 Returns a list of BASIS and class vector. The first item of the list is the BASIS or Kernel that is calculated. This is the supplied to the rvmbinary function. The second item is the class vector.
...
}
\references{
Tipping, M. E. and A. C. Faul (2003). Fast marginal likelihood maximisation for sparse Bayesian models. In C. M. Bishop and B. J. Frey (Eds.), Proceedings of the Ninth International Workshop on Artificial Intelligence and Statistics, Key West, FL, Jan 3-6.}
\author{
Robert Lowe
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(iris)
list=gaus(iris[1:100,],0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rvm}
\keyword{
classification
binary}% __ONLY ONE__ keyword per line
