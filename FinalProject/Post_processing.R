setwd("/Users/jmwerner1123/Dropbox/GitHub/Invisible/STA250/FinalProject")

cpu = list()
gpu = list()
vegas_cpu = list()
vegas_gpu = list()
times = list()
for(i in 1:6){
	cpu[[i]] = read.csv(paste0("normal_table_cpu_",i+1,".csv"))[,-1]
	gpu[[i]] = read.csv(paste0("normal_table_gpu_",i+1,".csv"))[,-1]
	vegas_cpu[[i]] = read.csv(paste0("normal_table_vegas_cpu_",i+1,".csv"))[,-1]
	vegas_gpu[[i]] = read.csv(paste0("normal_table_vegas_gpu_",i+1,".csv"))[,-1]
	times[[i]] = read.csv(paste0("times_",i+1,".csv"))
}




max_quantile = 3
table_quantiles = (1:(max_quantile * 100)) / 100

n = length(table_quantiles)
normal_vector = pnorm(table_quantiles, lower.tail = FALSE)
normal_table = matrix(normal_vector, ncol = 10, byrow = TRUE)
	
	

cpu_sse = lapply(cpu, FUN = function(x){
	return(sum((x - normal_table)^2))
})
gpu_sse = lapply(gpu, FUN = function(x){
	return(sum((x - normal_table)^2))
})
vegas_cpu_sse = lapply(vegas_cpu, FUN = function(x){
	return(sum((x - normal_table)^2))
})
vegas_gpu_sse = lapply(vegas_gpu, FUN = function(x){
	return(sum((x - normal_table)^2))
})

sse_frame = cbind(unlist(cpu_sse), unlist(gpu_sse), unlist(vegas_cpu_sse), unlist(vegas_gpu_sse))

pdf("sse_plot.pdf")
plot(2:7, sse_frame[,1], ylim = range(sse_frame), col = "red", type = "l", xlab = "K", ylab = "SSE", lty = 2, lwd = 2)
points(2:7, sse_frame[,2], col = "blue", type = "l", lty = 3, lwd = 2)
points(2:7, sse_frame[,3], col = "green", type = "l", lty = 4, lwd = 2)
points(2:7, sse_frame[,4], col = "orange", type = "l", lty = 5, lwd = 2)
legend("topright", c("CPU", "GPU", "Vegas CPU", "Vegas GPU"), col = c("red", "blue", "green", "orange"), lty = 2:5, lwd = rep(2,4))
dev.off()






time_frame = c()
for(k in 1:6){
	time_frame = rbind(time_frame,times[[k]]$elapsed)
}
pdf("time_plot.pdf")
plot(2:7, time_frame[1:6,1], ylim = range(time_frame), col = "red", type = "l", lty = 2, lwd = 2, ylab = "Time (Seconds)", xlab = "K")
points(2:7, time_frame[1:6,2], col = "blue", type = "l", lty = 3, lwd = 2)
points(2:7, time_frame[1:6,3], col = "green", type = "l", lty = 4, lwd = 2)
points(2:7, time_frame[1:6,4], col = "orange", type = "l", lty = 5, lwd = 2)
legend("topleft", c("CPU", "GPU", "Vegas CPU", "Vegas GPU"), col = c("red", "blue", "green", "orange"), lty = 2:5, lwd = rep(2,4))
dev.off()


pdf("example_plot.pdf")
x = (1:5000)/5000
plot(x,cos(x)*sin(x)+cos(4*x)-sin(10*x) + 1.6, type = "l", ylab = "f(x)", lwd = 2)
points(runif(25), runif(25,0,2.6), pch = 8, col = "RED")
dev.off()

pdf("importance_sampling.pdf")
x = (1:3000)/1000
plot(x,dexp(x,1.5), type = "l", col = "grey", lty = 3, ylab = "Density", lwd = 2)
points(x,dnorm(x), type = "l", lwd = 2, col = "Blue")
abline(v=0, lwd = 2)
legend("topright", c("Standard Normal Density", "Exponential Density"), col = c("blue", "grey"), lty = c(1,3), lwd = c(2,2))
dev.off()
