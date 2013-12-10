setwd("/Users/jmwerner1123/Dropbox/Skoo/F13/STA250/Assignment4")

data = read.csv("time_table_all.csv")

data = data[(1:6)*-2,]
data = as.matrix(data[,-1])

pdf("speed.pdf")
plot(data[,4], type = "l", main = "Speed Comparisons", ylab = "Seconds", xlab = "K", col = "green")
points(data[,4], pch = 17, col = "green")

points(1:7,as.numeric(data[,1]) + as.numeric(data[,2]) + as.numeric(data[,3]), type = "l", col = "RED")
points(1:7,as.numeric(data[,1]) + as.numeric(data[,2]) + as.numeric(data[,3]), pch = 8, col = "RED")
points(1:7, as.numeric(data[,2]), pch = 3, col = "BLUE")
legend("topleft", c("CPU", "GPU", "GPU Computation Only"), col = c("green", "red", "blue"), pch = c(17,8,3))
dev.off()