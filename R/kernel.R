rvmkernel <-
function (x, ...)
UseMethod ("rvmkernel")


rvmkernel.formula <- function (formula ,data = NULL, ...){
	
	if (!inherits(formula, "formula"))
	stop("method is only for formula objects")
	
	mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
	
	ret=rvmkernel.default(x, ...)
	
	ret$call <- match.call()
	ret$formula <- formula
	ret
}



rvmkernel.default <- function (x,kernel= "gaus",parameters=c(0.1), ...){
	
	x<-as.matrix(x)
	training=x
	parameters<-as.double(parameters)

	#Get kernel type
	kernel <- pmatch(kernel, c("aa","rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot"), 99) - 1
	
	if (kernel > 10) stop("wrong kernel specification!")	
	
	#Check for NA in data
	if(length(which(is.na(x)))>0) stop("NA in data!")
	
	if(kernel==0){
		rows=dim(x)[1]
		cols=dim(x)[2]
		out=matrix(-1,rows,rows)
		
		return=.C("aakernel",as.double(x),as.double(parameters),as.double(out),as.integer(rows),as.integer(cols),PACKAGE="rvmbinary")

		out=list(BASIS=matrix(return[[3]],rows,rows))
		out$type="aa"
		out$parameter=parameters

	}
	else{
		library(kernlab)
		
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
		
		out=list(BASIS=kernelMatrix(kern,x))
		out$type=name
		out$parameter=parameters

	}
		
	
	out$data <- training
	out$call <- match.call()
	class(out) <- "rvmkernel"
	out

}

print.rvmkernel <- function(x, ...){
	cat("Call:\n")
    print(x$call)
    cat("\nKernel:\n")
    print(x$type)
    cat("\nParameter:\n")
    print(x$parameter)}