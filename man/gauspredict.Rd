\name{gaus.predict}
\alias{gaus.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create the gauss kernel for RVM
}
\description{
  The \code{gaus.predict} function predicts a test data set using the Relevance Vectors and the training set.}
\usage{
gaus.predict(data,datatest, lambda,model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    The original training data with the class label as the last column
}
  \item{datatest}{
    The test data set with the class label as the last column
}
  \item{lambda}{
    The lambda parameter for the Gaussian kernel. Must be the same as used in \code{gaus}
}
  \item{model}{
    The model as returned by \code{rvmbinary}
}

}

\value{
 Returns a vector of the mapped probabilities of the class label. The value is between 0 and 1. The class boundary is at 0.5.
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
datatrain=rbind(iris[1:40,],iris[51:90,])
datatest=rbind(iris[41:50,],iris[91:100,])
list=gaus(datatrain,0.4)
model=rvmbinary(list[[1]],as.double(list[[2]])-1)
#Produces a list of the probabilities
out=gaus.predict(datatrain,datatest,0.4,model)
#Produce the predicted class
out>0.5
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rvm}
\keyword{
classification
binary}% __ONLY ONE__ keyword per line
