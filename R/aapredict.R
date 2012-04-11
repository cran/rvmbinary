
aa.predict <-
function (datatrain,datatest,lambda,model){
	
	vectors=model$vectors
	values=model$values
#remove the class from train and test
	class=datatrain[,length(datatrain)]
	datatrain=datatrain[,-length(datatrain)]
	datatrain=datatrain[vectors,]
	datatrain=data.matrix(datatrain)
	
	output=datatest[,length(datatest)]
	datatest=datatest[,-length(datatest)]
	datatest=data.matrix(datatest)	
	
	datarows=dim(datatrain)[1]
	testrows=dim(datatest)[1]
	cols=dim(datatrain)[2]
	BASIS=matrix(-1,dim(datatest)[1],length(vectors))

	out=.C("aapredict",as.double(datatest),as.double(datatrain),as.double(lambda),as.double(BASIS),as.integer(testrows),as.integer(cols),as.integer(datarows),PACKAGE="rvmbinary")
	
	BASIS=matrix(out[[4]],dim(datatest)[1],length(vectors))

	
	y=BASIS%*%values
	y=1/(1+exp(-y))
	return(y)
	
}
