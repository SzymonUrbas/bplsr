# useful function

.chol        <- function(x, ...) tryCatch(chol(x, ...), error=function(e)    {
	d          <- nrow(x)
	eigs       <- eigen(x, symmetric = TRUE)
	eval       <- eigs$values
	evec       <- eigs$vectors
	return(chol(x + evec %*% tcrossprod(diag(pmax.int(.Machine$double.eps, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec), ...))
}
)

Id = function(q){return(diag(rep(1,q)))}

mplot = function(X,digits = 2,...){
	if(any(X<0)){breakL = -1} else{breakL=0} 
	plot(X, digits = digits,key=list(side=3, cex.axis=0.75), breaks = max(abs(X))*c(breakL,1),
		main = deparse(substitute(X)),...)
}
coverage = function(pvals, CI){
	mean(abs(pvals-0.5)<=(CI/2))
}
coverage = Vectorize(coverage, vectorize.args = 'CI')

rltrgamma    <- function(n, shape, rate = 1, trunc = 1) {
	EV         <- as.list(environment())
	if(!all(vapply(EV, is.numeric,
		logical(1L))))        stop("All arguments must be numeric",              call.=FALSE)
	if(!all(lengths(EV) == 1))           stop("All arguments must be of length 1",          call.=FALSE)
	if(any(n   != floor(n)))             stop("'n' must be coercible to an integer",        call.=FALSE)
	if(shape   <= 0)                     stop("'shape' parameter must be strictly positive",   call.=FALSE)
	if(rate    <= 0)                     stop("'rate' parameter must be strictly positive", call.=FALSE)
	# if(trunc    < 0)                     stop("'trunc' cannot be negative",                 call.=FALSE)
	if(trunc <= 0){return(rgamma(n, shape, rate))}
	U          <- stats::runif(n)
	G          <- stats::pgamma(trunc, shape, rate)
	return(stats::qgamma(G + U * (1 - G),   shape, rate))
}

dltrgamma    <- function(x, shape, rate = 1, trunc = 1) {
	EV         <- as.list(environment())
	if(!all(vapply(EV, is.numeric,
		logical(1L))))        stop("All arguments must be numeric",              call.=FALSE)
	if(!all(lengths(EV) == 1))           stop("All arguments must be of length 1",          call.=FALSE)
	if(shape   <= 0)                     stop("'shape' parameter must be strictly positive",   call.=FALSE)
	if(rate    <= 0)                     stop("'rate' parameter must be strictly positive", call.=FALSE)
	if(trunc    < 0)                     stop("'trunc' cannot be negative",                 call.=FALSE)
	if(x<trunc){return(0)}

	return(stats::dgamma(x,   shape, rate)/(1-pgamma(trunc,shape, rate)))
}
dltrgamma = Vectorize(dltrgamma,vectorize.args = 'x')


# xs = seq(0,10,0.1)
# plot(xs, dltrgamma(xs,0.5,10,1), type = 'l')

std = function(X){
	Xnew = apply(X, 2 ,function(x){(x-mean(x))/sd(x)})
	names(Xnew) = names(X)
	return(Xnew)
}

hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

empty_plot = function(...){
	plot(0,0,type= 'l',...)
}
LinLog = function(y, y0){
	if(y>=y0){
		return(y)
	} else {
		return(y0*(log(y/y0)+1))
	}

}
LinLog = Vectorize(LinLog, vectorize.args = 'y')


invLinLog = function(ytil,y0){
	if(ytil>=y0){
		return(ytil)
	} else {
		return(y0*exp(ytil/y0-1))
	}
}
invLinLog = Vectorize(invLinLog, vectorize.args = 'ytil')
