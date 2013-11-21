setwd("/Users/jmwerner1123/Dropbox/GitHub/Invisible/STA250/HW2/Hive/summary/")

data<-read.table('results.txt', header=FALSE,sep="\001")
names(data) = c("Group", "Mean", "Variance")

#try 1
col1 = "blue"
col2 = "blue"
col_fun = colorRampPalette(c(col1,col2))
plot(data$Mean, data$Variance, col = col_fun(1000), pch = "x", xlab = "Group Mean", ylab = "Group Variance")
legend("topleft", c("Group 1", "Group 1000"), col = c(col1,col2), pch = c(1,1))

#try2
smoothScatter(data$Mean, data$Variance, nbin = 300, colramp = colorRampPalette(c(col1,col2)), 
	xlab = "Group Mean", ylab = "Group Variance")

#try3
pdf("hive_graph.pdf")
plot(data$Mean, data$Variance, col = "#0030FF4A", pch = 16, xlab = "Group Mean", ylab = "Group Variance")
dev.off()

