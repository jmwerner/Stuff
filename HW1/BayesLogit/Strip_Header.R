setwd("~/Stuff/HW1/BayesLogit/results")

files = list.files()

for(i in 1:length(files)){
	temp = read.csv(files[i])
	write.table(temp, file = files[i], row.names = FALSE, col.names = FALSE, sep = ",")
}

q("no")
