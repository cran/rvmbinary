\name{rvmbinary}
\alias{rvmbinary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binary Classification using RVM
}
\description{
  The \code{rvmbinary} function applies the RVM algorithm to a calculated kernel in order to produce a classification model.
}
\usage{
rvmbinary(BASIS, class,maxits=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASIS}{
    The BASIS or Kernel for learning as calculated by Gaussian or AA Kernel.
}
  \item{class}{
    The classification of each datapoint as a value 0 or 1. The class labels here must be coerced to double and be 0 or 1. 
}
  \item{maxits}{
    The maximum number of iterations. This should be set to a reasonably high number. Useful if the RVM algorithm cannot find an optimum value quickly enough.
}
}

\value{
 The function returns the Relevance Vector Numbers and their associated weights.
...
}
\references{
Tipping, M. E. and A. C. Faul (2003). Fast marginal likelihood maximisation for sparse Bayesian models. In C. M. Bishop and B. J. Frey (Eds.), Proceedings of the Ninth International Workshop on Artificial Intelligence and Statistics, Key West, FL, Jan 3-6.

Robert Lowe, Hamse Yussuf Mussa, John Mitchell, and Robert Glen. Classifying Molecules Using a Sparse Probabilistic Kernel Binary Classifier. J. Chem. Inf. Model. DOI: 10.1021/ci200128w
}
\author{
Robert Lowe
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(iris)
list=gaus(iris[1:100,],0.4)
rvmbinary(list[[1]],as.double(list[[2]])-1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rvm}
\keyword{
classification
binary}% __ONLY ONE__ keyword per line
