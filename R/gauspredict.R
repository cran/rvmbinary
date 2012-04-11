gaus.predict <-
function (data,datatest,lambda,model){
	
	vectors=model$vectors
	values=model$values
	
#strip class label
	class=data[,dim(data)[2]]
	data=data[,-dim(data)[2]]
#only require relvance vectors
	data=data[vectors,]
	Y=data.matrix(data)

	trrows=dim(data)[1]
	trcols=dim(data)[2]
#strip class label	
	output=datatest[,dim(datatest)[2]]
	datatest=datatest[,-dim(datatest)[2]]
	X=data.matrix(datatest)		
	
	tsrows=dim(datatest)[1]
	tscols=dim(datatest)[2]
	
	X2=matrix(NA,tsrows,trrows)
	Y2=matrix(NA,tsrows,trrows)
	
	for(i in 1:trrows){
		sumofsquares=sum(as.double(data[i,])^2)
		Y2[,i]=rep(sumofsquares,tsrows)
	}

	for(i in 1:tsrows){
		sumofsquares=sum(as.double(datatest[i,])^2)
		X2[i,]=rep(sumofsquares,trrows)
	}	
	
	D2=X2+Y2-2*X%*%t(Y)
	
	BASIS=exp(-D2*(lambda))
	
	
	y=BASIS%*%values
	y=1/(1+exp(-y))
	return(y)
}

