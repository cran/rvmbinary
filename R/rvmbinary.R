rvmbinary <-
function (x, ...)
UseMethod ("rvmbinary")

rvmbinary.formula <- function (formula ,data = NULL, ...){

	if (!inherits(formula, "formula"))
		stop("method is only for formula objects")

	mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)
	
	ret=rvmbinary.default(x,y, ...)
	
	ret$call <- match.call()
	ret$formula <- formula
	ret
}



rvmbinary.default <- function (x,y,kernel="rbfdot",parameters=c(0.1),iterations=100, noisevar=0.1,minmaxdiff=1e-3, ...){
	
	relevance=NULL

	formula <- inherits(x, "formula")

	#Check whethe x is an already computed kernel
	if(class(x)=="rvmkernel"){
		BASIS=x$BASIS
		training=data.matrix(x$data)
		relevance$kernel=x$type
		relevance$kernelparameter=x$parameter
	}
	else{
	#compute the kernel
		x=rvmkernel(x,kernel,parameters)
		BASIS=x$BASIS
		training=data.matrix(x$data)
		relevance$kernel=x$type
		relevance$kernelparameter=x$parameter
	}
	
	#Check y are factors and turn to 0 and 1
	if(is.factor(y)){
		y=as.factor(as.character(y))
		if(length(levels(y))>2) stop("More than 2 classes given")
	
		y=as.double(y==levels(as.factor(y))[2])
	
		print("Calculating model...")
		relv=rep(-1,dim(BASIS)[1])
		relvals=rep(-1,dim(BASIS)[1])
		
		if(length(which(is.na(BASIS)))>0) stop("NA in BASIS!")
		if(length(which(is.nan(BASIS)))>0) stop("NaN in BASIS!")
		if(length(which(BASIS==Inf))) stop("Inf in BASIS!")
		
		out=.C("initial",as.double(BASIS),as.double(y),as.integer(dim(BASIS)[1]),as.integer(dim(BASIS)[2]),as.integer(relv),as.double(relvals),as.integer(iterations),as.integer(1),as.double(noisevar),as.double(minmaxdiff),PACKAGE="rvmbinary")
		relv=NULL
		relvals=NULL
		for(i in 1:length(out[[5]])){
			if(out[[5]][i]!=-1){
				relv=c(relv,out[[5]][i])
				relvals=c(relvals,out[[6]][i])
			}
		}
		
		relevance$type="Classification"
	}
	else{
		print("Regression Case")
		print("Calculating model...")
		relv=rep(-1,dim(BASIS)[1])
		relvals=rep(-1,dim(BASIS)[1])
		out=.C("initial",as.double(BASIS),as.double(y),as.integer(dim(BASIS)[1]),as.integer(dim(BASIS)[2]),as.integer(relv),as.double(relvals),as.integer(iterations),as.integer(0),as.double(noisevar),as.double(minmaxdiff),PACKAGE="rvmbinary")
		relv=NULL
		relvals=NULL
		for(i in 1:length(out[[5]])){
			if(out[[5]][i]!=-1){
				relv=c(relv,out[[5]][i])
				relvals=c(relvals,out[[6]][i])
			}
		}
		
		relevance$type="Regression"

	}

	relevance$vectors=training[c(relv+1),]
	relevance$values=relvals
	relevance$nRV=length(relvals)
	relevance$call=match.call()
	
	class(relevance)="rvmbinary"
	
	relevance
}

predict.rvmbinary <- function(object,newdata, ...){
	
	if(length(object$values)<1) stop ("Model is empty")
	
	modeldata=as.data.frame(object$vectors)
	testdata=as.data.frame(newdata)
	if (object$nRV == 1){
		modeldata = as.data.frame(matrix(object$vectors,1))
		print("******* Warning predicting with a single relevance vector *******")
	}
	
	if(ncol(modeldata) != ncol(testdata))
		stop("Test data does not match model")
	
	kernel <- pmatch(object$kernel, c("aa","rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot"), 99) - 1

	datatrain=data.matrix(object$vectors)
	if (object$nRV ==1)
		datatrain = data.matrix(matrix(object$vectors,1))

	if(!is.null(object$formula)){
	## model has been fitted using formula interface
		datatest <- model.matrix(object$formula, newdata)
	}
	else{
		datatest=data.matrix(newdata)	
	}
	
	values=object$values
	parameters=object$kernelparameter
	

	
	if(kernel==0){
		
		datarows=dim(datatrain)[1]
		testrows=dim(datatest)[1]
		cols=dim(datatrain)[2]
		BASIS=matrix(-1,dim(datatest)[1],dim(datatrain)[1])
		
		
		out=.C("aapredict",as.double(datatest),as.double(datatrain),as.double(object$kernelparameter),as.double(BASIS),as.integer(testrows),as.integer(cols),as.integer(datarows),PACKAGE="rvmbinary")
		
		BASIS=matrix(out[[4]],dim(datatest)[1],dim(datatrain)[1])
		
		
	}
	else {
		
		if(kernel==1){
			kern <- rbfdot(sigma = parameters)
			name="rbfdot"
		}
		else if(kernel==2){
			if(length(parameters)!=3) stop("For polydot please define parameters=c(degree,scale,offset)")
			kern=polydot(degree = parameters[1], scale = parameters[2], offset = parameters[3])
			name="polydot"
		}
		else if(kernel==3){
			if(length(parameters)!=2) stop("For tanhdot please define parameters=c(scale,offset)")
			kern=tanhdot(scale = parameters[1], offset = parameters[2])
			name="tanhdot"
		}
		else if(kernel==4){
			kern=vanilladot()
			name="vanilladot"
		}
		else if(kernel==5){
			kern=laplacedot(sigma = parameters)
			name="laplacedot"
		}
		else if(kernel==6){
			if(length(parameters)!=3) stop("For besseldot please define parameters=c(sigma,order,degree)")
			kern=besseldot(sigma = parameters[1], order = parameters[2], degree = parameters[3])
			name="besseldot"
		}
		else if(kernel==7){
			if(length(parameters)!=2) stop("For anovadot please define parameters=c(sigma,degree)")
			kern=anovadot(sigma = parameters[1], degree = parameters[2])
			name="anovadot"
		}
		else if(kernel==8){
			kern=splinedot()
			name="splinedot"
		}


		BASIS=kernelMatrix(kern,datatest,datatrain)
	}
	
	y=BASIS%*%values
		
	if(object$type=="Classification"){
		y=1/(1+exp(-y))

	}

	y
}
	
	

print.rvmbinary <- function(x, ...){
	cat("Call:\n")
    print(x$call)
    cat("\nNumber of Relevance Vectors:\n")
    print(x$nRV)	
}


	
	

