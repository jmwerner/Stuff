####################################
# Author: Jeremy Werner            #
# Assignment: STA 250 Assignment 3 #
# Due date: November 27, 2013      #
####################################
setwd("/Users/jmwerner1123/Dropbox/Skoo/F13/STA250/Assignment3")
#####################################################################################
#tick and tock functions for measuring time (circa matlab) modified to meet my needs
#source: http://stackoverflow.com/questions/1716012/stopwatch-function-in-r
tick <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")){
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}
tock <- function(){
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   #print(toc - tic)
   invisible(toc)
   return(toc - tic)
}
#####################################################################################

#Bisection algorithm for finding a root of a function
bisection = function(fun_input, interval_input, tol = 1e-6, maxiter = 1e6, verbose = FALSE){
	tick()
	if(verbose){
		cat(paste0("Begin Bisection Optimization Function...\n", 
			"Interval: (", interval_input[1],",", interval_input[2], 
			") Tol: ", tol, "\n\n"))
	}
	iter = 1
	tol_val = 1
	lower = interval_input[1]
	upper = interval_input[2]
	while((iter <= maxiter) & (tol_val > tol)){
		if(verbose & !(iter %% floor(maxiter/15))){
			cat(paste0("Iteration ", iter, "\tCurrent c: ", c, " f(c) = ",
				fun_input(c), "\n"))
		}
		
		c = (lower + upper)/2
		fun_value = fun_input(c)
		
		if(fun_value*fun_input(upper) > 0){
			upper = c
		}else{
			lower = c
		}
		
		tol_val = abs(fun_value)
		iter = iter + 1
	}
	
	if(verbose & (iter == maxiter+1)){
		cat(paste("\nERROR\n\nMaximum iterations met, current Tol value:", tol_val, "\n", 
			"  Current c: ", c, " f(c) = ", fun_input(c), "\n"))
	}
	
	time = tock()
	if(verbose){
		cat(paste("\nBisection function converged in", iter, 
			"iterations\nElapsed time:", round(time, 4), "seconds \n"))
	}
	return(c)
}#bisection

#Newton-Raphson algorithm to find the root of a function
newton_raphson = function(fun_input, derivative_input, starting_point, tol = 1e-6, maxiter = 1e6, verbose = FALSE){
	tick()
	if(verbose){
		cat(paste0("Begin Newton-Raphson Optimization Function...\n", 
			"Statring point: ", starting_point, " Tol: ", tol, "\n\n"))
	}
	iter = 1
	tol_val = 1
	xt = starting_point
	while((iter <= maxiter) & (tol_val > tol)){
		if(verbose & !(iter %% floor(maxiter/15))){
			cat(paste0("Iteration ", iter, "\tCurrent xt: ", xt, " f(xt) = ",
				fun_input(xt), "\n"))
		}
		xt = xt - fun_input(xt)/derivative_input(xt)
		fun_value = fun_input(xt)
		
		tol_val = abs(fun_value)
		iter = iter + 1
	}
	
	if(verbose & (iter == maxiter+1)){
		cat(paste("\nERROR\n\nMaximum iterations met, current f(xt) value:", 
			fun_value, "Current xt: ", xt, "\n"))
	}
	
	time = tock()
	if(verbose){
		cat(paste("\nNewton-Raphson function converged in", iter, 
			"iterations\nElapsed time:", round(time, 4), "seconds \n"))
	}
	return(xt)
}#newton_raphson


#First and second derivatives of the log likelihood function (we need the root of the
# first derivative to maximize the log likelihood)
log_likelihood_d1 = function(x){return(125/(2+x)-38/(1-x)+34/x)}
log_likelihood_d2 = function(x){return(-125/(2+x)^2-38/(1-x)^2-34/x^2)}


x_vals = seq(-4,4,.01)
pdf("first_derivative.pdf")
plot(x_vals, log_likelihood_d1(x_vals), type = "n", ylab = expression(paste("(log(L(",lambda,")))'")), xlab = expression(paste(lambda)), ylim = c(-1000,1000), main = "First Derivative of Proportional Log Likelihood Function")
abline(0,0, col = "GREY", lwd = 2)
lapply(c(-2,1,0), FUN = function(x){abline(v=x, col = "BLUE", lty = 2, lwd = 2)})
points(x_vals, log_likelihood_d1(x_vals), type = "l", lwd = 3)
dev.off()


bisection_result = bisection(log_likelihood_d1, c(0,1), verbose = TRUE)

n_r_result = newton_raphson(fun_input = log_likelihood_d1, derivative_input = log_likelihood_d2, starting_point = .5, verbose = TRUE)

