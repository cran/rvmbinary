\name{aa.predict}
\alias{aa.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Used to predict test data using the AA kernel
}
\description{
  The \code{aa.predict} function predicts a test data set using the relevance vectors and the training set.}
\usage{
aa.predict(datatrain,datatest, lambda,model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datatrain}{
    The original training data with the class label as the last column
}
  \item{datatest}{
    The test data set with the class label as the last column
}
  \item{lambda}{
    The lambda parameter of the aa kernel
}
  \item{model}{
    The model as returned by \code{rvmbinary}
}

}

\value{
 Returns a vector of the mapped probabilities of the class label. The value is between 0 and 1. 
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
print("Calculate Kernel")
kernel=aa(traindata,0.6)
#Compute the Relevance Vectors. 
#Kernel is a list of the kernel [[1]] and the associated class labels [[2]]. 
#Class labels must be 0 and 1.
model=rvmbinary(kernel[[1]],as.double(kernel[[2]])-1,100)
#Calculate the predictions for testset
y=aa.predict(traindata,testdata,0.6,model)
#Calculate classes
y>0.5
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rvm}
\keyword{
classification
binary}% __ONLY ONE__ keyword per line
