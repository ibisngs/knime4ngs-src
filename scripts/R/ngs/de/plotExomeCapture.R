########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
if (!require(argparse)) {
	install.packages("argparse", repos="http://cran.rstudio.com")
	library("argparse")
}
parser <- ArgumentParser(prog="plotExomeCapture.R", description="Can be use to plot coverage distribution.")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store", dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## INFILE
parser$add_argument("--in", type="character", action="store", dest="infile", required=TRUE, help="path to input file", metavar="<path>")
parser$add_argument("--infile", type="character", action="store", dest="in", required=TRUE, help="path to input file", metavar="<path>")


## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
infile <- args$infile

folderPath <- paste(unlist(strsplit(infile,"/"))[1:length(unlist(strsplit(infile,"/")))-1],collapse="/")
outfile = paste(folderPath,"/capture_efficiency.pdf",sep="")

pdf(outfile)
data = read.table(infile,sep="\t",row.names=1)

#Switch rows and cols
data = data.frame(t(data))	

#Add Coverage "header" and reorder
data$Cov=c(0:100)
data=data[c(ncol(data),1:ncol(data)-1)]


colors = (ncol(data))-1
plot(data[,1],data[,2],type="l",col=rainbow(colors)[1],xlab="Coverage",ylab="% of target region covered",main="Capture Efficiency",ylim=c(0,100))
if(ncol(data)>2) {
	for(i in 3:ncol(data)){
		lines(data[,1],data[,i],col=rainbow(colors)[i-1])
	}
}

dev.off()
