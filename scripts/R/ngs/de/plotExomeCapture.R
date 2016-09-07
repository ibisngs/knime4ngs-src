#  Copyright (C) 2016 the Knime4NGS contributors.
#  Website: http://ibisngs.github.io/knime4ngs
#  
#  This file is part of the KNIME4NGS KNIME extension.
#  
#  The KNIME4NGS extension is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
