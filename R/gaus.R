gaus <-
function (data,lambda){
	class=data[,length(data)]
	data=data[,-length(data)]
	data=data.matrix(data)
	X2=matrix(NA,dim(data)[1],dim(data)[1])
	Y2=matrix(NA,dim(data)[1],dim(data)[1])
	
	for(i in 1:dim(data)[1]){
		sumofsquares=sum(data[i,]^2)
		X2[i,]=rep(sumofsquares,dim(data)[1])
		Y2[,i]=rep(sumofsquares,dim(data)[1])
	}
	D2=X2+Y2-2*data%*%t(data)
	rm(X2,Y2,data)
	BASIS=exp(-D2*(lambda))
	rm(D2)
	list(BASIS,class)
}

