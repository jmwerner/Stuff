####################################
# Author: Jeremy Werner            #
# Assignment: STA 205 HW1          #
# Due date: October 28, 2013       #
####################################

#Breast cancer data analysis based off of other Bayes Logit script

#The script in its current state does not produce results that I am 100% pleased with
# However, it produces some nice results and I learned a lot while writing it

#Packages needed
library(MASS)
library(xtable)


setwd("/Users/jmwerner1123/Dropbox/GitHub/STA250/Stuff/HW1/BayesLogit/breast_cancer")
data = read.table("breast_cancer.txt", header = TRUE)

levels(data$diagnosis) = c(0,1) #slick way to reassign 1 and 0 to M and B


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
"bayes.logreg" <- function(m,y,X,beta.0,Sigma.0.inv, niter=10000,burnin=1000, print.every=1000,retune=100, tuning_loops = 1, verbose=TRUE){
                     
	if(verbose){cat(paste("Begin bayes.logreg function with", burnin,
		"burn in iterations \n and", niter, "iterations following the burn in period\nCurrent time:", date(), "\n\n"))}
	begin_time = Sys.time()
	
	Beta_0 = beta.0
		
	#Proposal variance starting point 
	v = .4

	#Burn-in Period to get a somewhat close v for all parameters
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
				v = v*1.5
				if(verbose){cat(paste("Retuning proposal variance to", round(v,4), "from", round(v1,4), "\n"))}
			}else{
				if(percent_1 < .3){
					v1=v
					v = v/1.5
					if(verbose){cat(paste("Retuning proposal variance to", round(v,4), "from", round(v1,4), "\n"))}
				}
			}
		}
	}#burnin loop #1
	

	
	
	
	acceptance_vector_burnin = c()
	beta_vector_burnin = c()
	
#################################################################################	
#burn-in period to tune each parameter's proposal variance closer to where it needs to be
	
	v_vect = rep(v, length(Beta_0)) 
	
 	for(looper in 1:tuning_loops){
		cat(paste("Starting tuning loop", looper, "of", tuning_loops, "\n"))
		for(par in 1:length(Beta_0)){
			for(k in 1:burnin){
			
				Beta_1 = Beta_0
				Beta_1[par] = rnorm(1,Beta_0[par], sqrt(v_vect[par]))
		
				alpha = min(0,pi_fun(Beta_1, y, X, m, beta.0, Sigma.0.inv) -
						pi_fun(Beta_0, y, X, m, beta.0, Sigma.0.inv))
				
				if(alpha >= log(runif(1))){
					Beta_0 = Beta_1
					acceptance_vector_burnin = c(acceptance_vector_burnin, TRUE)
				}else{
					acceptance_vector_burnin = c(acceptance_vector_burnin, FALSE)
				}
			
				beta_vector_burnin = rbind(beta_vector_burnin, Beta_0)		
			
				#Retune
				if(!(k %% retune)){
					current_index = length(acceptance_vector_burnin)
					percent_1 = sum(acceptance_vector_burnin[(current_index-retune+1):current_index])/retune
					if(percent_1 > .85){
						v_vect[par] = v_vect[par]*1.15
						if(verbose){cat(paste("Retuning proposal variance up for beta", par, "\n"))}
					}else{
						if(percent_1 < .15){
							v_vect[par] = v_vect[par]/1.15
							if(verbose){cat(paste("Retuning proposal variance down for beta", par, "\n"))}
						}
					}
				}
			}	
		}
	} #tuning looper
	

	
	v_mat = diag(v_vect)
	
	#Post burn-in Period
	if(verbose){cat("\nBurn in period complete. \n \n")}
	
	acceptance_vector = c()
	beta_vector = c()
	
	for(t in 1:niter){
		if(!(t %% print.every)){cat(paste("Algorithm is", t/niter*100, 
			"percent complete and has been running for", round(as.numeric(Sys.time()-begin_time),2),"seconds\n"))}
		Beta_1 = mvrnorm(1,Beta_0, v_mat)
		
		
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
	
	
	if(verbose){
		plot.new()
		par(mfrow = c(3,4))
		for(i in 1:(ncol(data))){
			plot(beta_vector[,i], type = "l", main = paste("Beta", i))
		}
		hist(as.numeric(acceptance_vector))
	}
	end_time = Sys.time()
	
	if(verbose){cat(paste("\nAcceptance rate:", sum(acceptance_vector)/length(acceptance_vector), "\nElapsed time:", round(as.numeric(end_time - begin_time),3), "seconds"))} #ending verbose output
	
	return(beta_vector)
}#bayes.logreg
##############################################################

############
#Main block#
############


#Standardize X columns for stability
X_in = cbind(rep(1,nrow(data)),scale(as.matrix(data[,1:10])))
#X_in = cbind(rep(1,nrow(data)),as.matrix(data[,1:10]))

starting_sigma = 1000*diag(ncol(data))
starting_b0 = apply(X_in, 2, "mean")


#niter, print.every, and retune defaults are acceptable and will not be changed here
betas = bayes.logreg(rep(1,nrow(data)), as.numeric(data$diagnosis)-1, X_in, starting_b0, solve(starting_sigma),burnin=800, niter = 15000, retune = 100, tuning_loops = 7, verbose = TRUE)


means = apply(betas, 2, mean)

autocorrs = apply(betas, 2, function(x)as.numeric(unlist(acf(x, plot=FALSE)[2])))[1,]
print(xtable(data.frame(autocorrs),digits=4), row.names = FALSE)

#Posterior predictive checks using the mean of our returned beta matrix
simulated_data = c()
for(j in 1:nrow(betas)){
	u = exp(X_in%*%betas[j,])/(1+exp(X_in%*%betas[j,]))		
	simulated_data = cbind(simulated_data,rbinom(length(u),1,u))
	if(!(j%%1000)){cat(paste(j/nrow(betas), "Percent done\n"))}
}

stats = apply(simulated_data,2, "mean")
pdf("ppc.pdf")
hist(stats, main = "Posterior Predictive Check for Mean")
abline(v=mean(as.numeric(data$diagnosis)-1), col = "red", lwd = 2, lty = 2)
legend("topleft", c("Real Data"), lty = c(2), lwd = c(2), col = c("red"))
dev.off()


stats = apply(simulated_data,2, "sd")
pdf("ppc_sd.pdf")
hist(stats, main = "Posterior Predictive Check for Standard Deviation")
abline(v=sd(as.numeric(data$diagnosis)-1), col = "red", lwd = 2, lty = 2)
legend("topleft", c("Real Data"), lty = c(2), lwd = c(2), col = c("red"))
dev.off()


