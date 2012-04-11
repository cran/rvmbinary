\name{aa}
\alias{aa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create the gauss kernel for RVM
}
\description{
  The \code{aa} function applies the Aitchison-Aitken(AA) kernel to the data}
\usage{
aa(data, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    The data with the class label as the last column. The data must have binary features.
}
  \item{lambda}{
    The lambda parameter of the AA kernel. (values - $>0.5$ and $<1.0$)
}
}

\value{
 Returns a list of BASIS and class vector
...
}
\references{
Aitchison, J. and Aitken, C.G.G. (1976). Multivariate binary discrimination by kernel method. Biometrika 63, 413-420.}
\author{
Robert Lowe
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(mlbench)

#DNA database is binary features three class problem
data(DNA)
#select datapoints with class "ie" and "ei"
datapoints=c(which(DNA[,181]=="ie"),which(DNA[,181]=="ei"))
data=DNA[datapoints,]

#Create 80% training and 20% test randomly
trainpoints=sample(1:dim(data)[1],as.integer(dim(data)[1]*0.8))
testpoints=setdiff(1:dim(data)[1],trainpoints)
traindata=data[trainpoints,]
testdata=data[testpoints,]

#Compute the AA Kernel
kernel=aa(traindata,0.6)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rvm}
\keyword{
classification
binary}% __ONLY ONE__ keyword per line
