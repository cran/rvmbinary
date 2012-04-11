rvmbinary <-
function (BASIS,class,maxits=10){
	print("Calculating model...")
	relv=rep(-1,dim(BASIS)[1])
	relvals=rep(-1,dim(BASIS)[1])
	
#Check to see if class labels are 1 and 0
	if(length(intersect(which(as.double(class)!=0), which(as.double(class)!=1)))!=0){
		print("The class label supplied are not 0s and 1s. Exiting ......")
		return(0)
	}
	
	out=.C("initial",as.double(BASIS),as.double(class),as.integer(dim(BASIS)[1]),as.integer(dim(BASIS)[2]),as.integer(relv),as.double(relvals),as.integer(maxits),PACKAGE="rvmbinary")
	relv=NULL
	relvals=NULL
	for(i in 1:length(out[[5]])){
		if(out[[5]][i]!=-1){
			relv=c(relv,out[[5]][i])
			relvals=c(relvals,out[[6]][i])
		}
	}
	relevance=NULL
	relevance$vectors=(relv+1)
	relevance$values=relvals
	return(relevance)
}

