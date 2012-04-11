aa <-
function (data,lambda){
	class=data[,length(data)]
	data=data[,-length(data)]
	data=data.matrix(data)
	rows=dim(data)[1]
	cols=dim(data)[2]
	out=matrix(-1,rows,rows)
	
	return=.C("aakernel",as.double(data),as.double(lambda),as.double(out),as.integer(rows),as.integer(cols),PACKAGE="rvmbinary")
	
	BASIS=matrix(return[[3]],rows,rows)
	list(BASIS,class)
}