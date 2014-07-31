########################################################################################################################################
## PARSE ARGS
########################################################################################################################################

## parse
args <- commandArgs(TRUE)

infile <- args[8]
folderPath <- paste(unlist(strsplit(infile,"/"))[1:length(unlist(strsplit(infile,"/")))-1],collapse="/")
outfile = paste(folderPath,"/capture_efficiency.pdf",sep="")

pdf(outfile)
data = read.table(infile,sep="\t",row.names=1)

#Switch rows and cols
data = data.frame(t(data))	

#Add Coverage "header" and reorder
data$Cov=c(0:100)
data=data[c(3,1,2)]


colors = (ncol(data))-1
plot(data[,1],data[,2],type="l",col=rainbow(colors)[1],xlab="Coverage",ylab="% of target region covered",main="Capture Efficiency",ylim=c(0,100))
for(i in 2:ncol(data)){
	lines(data[,1],data[,i],col=rainbow(colors)[i])
}

dev.off()
