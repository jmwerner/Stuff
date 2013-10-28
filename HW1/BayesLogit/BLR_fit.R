####################################
# Author: Jeremy Werner            #
# Assignment: STA 205 HW1          #
# Due date: October 28, 2013       #
####################################

#This script was written to run as an array job on Gauss. All diagnostic tools
# were commented out, but left in to show some of my thought process throughout 
# development 

#Get command line arguments
args <- commandArgs(TRUE)
cat(paste("Command Line Arguments:", args, "\n"))


#Packages needed
library(MASS)

###
sim_start <- 1000   #starting number for jobs
length.datasets <- 200
###

if(length(args)==0){
	sinkit <- FALSE
	jobid <- sim_start + 1
	set.seed(1330931)
}else{
	# Sink output to file?
	sinkit <- TRUE
	# Decide on the job number, usually start at 1000:
	jobid <- sim_start + as.numeric(args[1])
	# Set a different random seed for every job number!!!
	set.seed(321*jobid + 1234567)
}


#Used for local fitting
#setwd("/Users/jmwerner1123/Dropbox/GitHub/Invisible/STA250/HW1")
#data = read.csv("blr_data_1001.csv")
#pars = read.csv("blr_pars_1001.csv")

#Reading data and pars files on Gauss
data_path = "/home/jwerner/Stuff/HW1/BayesLogit/data/"
data = read.csv(paste0(data_path,"blr_data_",jobid,".csv"))

#Output directory for .csv files
output_path = "/home/jwerner/Stuff/HW1/BayesLogit/results/"

##############################################################
#functions
##############################################################
##Pi Function, because pi is fun (obviously)
# beta is vector of coefficients
# y is response
# X is design matrix (size n x p)
# m is vector of n for each trial
# mu is vector of mu values
# sigma.inv is input inverse covariance matrix
"pi_fun" <- function(beta, y, X, m, mu, sigma.inv){
	return(t(beta - mu) %*% sigma.inv %*% (beta-mu) + 
		t(y) %*% (X %*% beta) - t(m) %*% log(1 + exp(X %*% beta)))
}#pi_fun
##############################################################
##bayes.logreg function
# beta is vector of starting coefficients (guesses)
# y is response
# X is design matrix (size n x p)
# m is vector of n for each trial
# Sigma.0.inv is input of inverse covariance matrix
# niter is the iterations to run the MH algorithm following the burnin
# burnin is the number of iterations to burn in the algorithm (and retune proposal variance)
# print.every is the number of iterations between update printings to the user
# retune is the number of iterations between retunings of proposal variance during burnin period
# verbose is the logical to print updates or no updates
##############################################################
"bayes.logreg" <- function(m,y,X,beta.0,Sigma.0.inv, niter=10000,burnin=1000, print.every=1000,retune=100, verbose=TRUE){
                     
	if(verbose){cat(paste("Begin bayes.logreg function with", burnin,
		"burn in iterations \n and", niter, "iterations following the burn in period\nCurrent time:", date(), "\n\n"))}
	begin_time = Sys.time()
	
	Beta_0 = beta.0
		
	#Proposal variance starting point
	v = .1
	
	acceptance_vector_burnin = c()
	beta_vector_burnin = c()
	
	#Burn-in Period
	for(i in 1:burnin){
		#proposal density of a normal variable centered at mu
		Beta_1 = mvrnorm(1,Beta_0, v*diag(length(Beta_0)))
		
		alpha = min(0,pi_fun(Beta_1, y, X, m, beta.0, Sigma.0.inv) -
					pi_fun(Beta_0, y, X, m, beta.0, Sigma.0.inv))
				
		if(alpha >= log(runif(1))){
			Beta_0 = Beta_1
			acceptance_vector_burnin[i] = TRUE
		}else{
			acceptance_vector_burnin[i] = FALSE
		}
		
		beta_vector_burnin = rbind(beta_vector_burnin, Beta_0)		
		
		#Retune
		if(!(i %% retune)){
			percent_1 = sum(acceptance_vector_burnin[(i-retune+1):i])/retune
			if(percent_1 > .7){
				v1=v
				v = v*1.4
				if(verbose){cat(paste("Retuning proposal variance to", round(v,4), "from", round(v1,4), "\n"))}
			}else{
				if(percent_1 < .4){
					v1=v
					v = v/1.4
					if(verbose){cat(paste("Retuning proposal variance to", round(v,4), "from", round(v1,4), "\n"))}
				}
			}
		}
	}#burnin loop
	
	#if(verbose){
	#	pdf(paste0("burnin_diagnostics_",jobid,".pdf"))
	#	par(mfrow = c(2,2))
	#	plot(beta_vector_burnin[,1], beta_vector_burnin[,2], type = "l")
	#	plot(beta_vector_burnin[,1], type = "l")
	#	plot(beta_vector_burnin[,2], type = "l")
	#	hist(as.numeric(acceptance_vector_burnin))
	#	dev.off()
	#}
	
	#Post burn-in Period
	if(verbose){cat("\nBurn in period complete. \n \n")}
	
	acceptance_vector = c()
	beta_vector = c()
	
	for(t in 1:niter){
		if(!(t %% print.every)){cat(paste("Algorithm is", t/niter*100, 
			"percent complete and has been running for", round(as.numeric(Sys.time()-begin_time),2),"seconds\n"))}
		Beta_1 = mvrnorm(1,Beta_0, v*diag(length(Beta_0)))
		
		alpha = min(0,pi_fun(Beta_1, y, X, m, beta.0, Sigma.0.inv) -
					pi_fun(Beta_0, y, X, m, beta.0, Sigma.0.inv))
				
		if(alpha >= log(runif(1))){
			Beta_0 = Beta_1
			acceptance_vector[t] = TRUE
		}else{
			acceptance_vector[t] = FALSE
		}	
		
		beta_vector = rbind(beta_vector, Beta_0)					
	}
	
	#if(verbose){
	#	pdf(paste0("diagnostics_",jobid,".pdf"))
	#	par(mfrow = c(2,2))
	#	plot(beta_vector[,1], beta_vector[,2], type = "l")
	#	plot(beta_vector[,1], type = "l")
	#	plot(beta_vector[,2], type = "l")
	#	hist(as.numeric(acceptance_vector))
	#	dev.off()
	#}
	
	end_time = Sys.time()
	
	if(verbose){cat(paste("\nAcceptance rate:", sum(acceptance_vector)/length(acceptance_vector), "\nElapsed time:", round(as.numeric(end_time - begin_time),3), "seconds"))} #ending verbose output
	
	return(beta_vector)
}#bayes.logreg
##############################################################

############
#Main block#
############

#Selected mu and starting sigma values made up by me
Mu = c(0,1)
starting_sigma = diag(length(Mu))
starting_b0 = mvrnorm(1,Mu,starting_sigma)


#niter, print.every, and retune defaults are acceptable and will not be changed here
betas = bayes.logreg(data$n, data$y, cbind(data$X1, data$X2), starting_b0, solve(starting_sigma),burnin=2000, verbose = FALSE)


#Creation of .csv output file of quantiles for returned beta vector
probs = (1:99)/100
percentile_table = cbind(quantile(betas[,1], probs), quantile(betas[,2], probs))
write.table(percentile_table,file=paste0(output_path,"Bayes_Logit_Percentiles_",jobid,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
	

q("no")