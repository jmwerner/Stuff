setwd("/Users/jmwerner1123/Dropbox/Skoo/F13/STA250/Assignment4")

data = read.csv("all_times.csv")

cpu = as.numeric(as.character(data[c(2,5,8,11,14),4]))
gpu = as.numeric(as.character(data[c(1,4,7,10,13),4]))

pdf("speed_probit.pdf")
plot(gpu, type = "l", main = "Probit Speed Comparisons", ylab = "Seconds", xlab = "Data Set", col = "green")
points(gpu, pch = 17, col = "green")

points(1:5,cpu, type = "l", col = "RED")
points(1:5,cpu, pch = 8, col = "RED")
legend("topleft", c("CPU", "GPU"), col = c("red", "green"), pch = c(8,17))
dev.off()


