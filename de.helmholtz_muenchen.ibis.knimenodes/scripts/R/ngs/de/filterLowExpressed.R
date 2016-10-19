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
parser <- ArgumentParser(prog="filterLowExpressed.R", description="Can be used to remove very low expressed genes based on their counts.")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store", dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument("--countTable", type="character", action="store", dest="file.count", required=TRUE, help="path to input count file", metavar="<path>")
parser$add_argument("--annotationFile", type="character", action="store", dest="file.annotation", required=TRUE, help="path to annotation file", metavar="<path>")
parser$add_argument("--output", type="character", action="store", dest="file.output", required=TRUE, help="path to output folder", metavar="<path>")

## ARGUMENTS
parser$add_argument("--keepReads", type="integer", action="store", dest="keepReads", required=TRUE, help="normalization method for size factors", metavar="<char>")
parser$add_argument("--keepFraction", type="double", action="store", dest="keepFraction", required=TRUE, help="normalization method for CPM values", metavar="<char>")
parser$add_argument("--bothConditionsSeperate", type="integer", action="store", dest="bothConditionsSeperate", required=TRUE, help="correction method for p-values", metavar="<char>")
parser$add_argument("--filterMode", type="integer", action="store", dest="filterMode", required=TRUE, help="filter mode", metavar="<char>")


## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)


########################################################################################################################################
## LOAD COUNTS
########################################################################################################################################
exprsPath <- args$file.count
pDataPath <- args$file.annotation

pData <- read.csv3(pDataPath)
exprs <- as.matrix(read.csv3(exprsPath))

pData <- as.data.frame(pData[colnames(exprs), ]) 
rownames(pData) <- colnames(exprs)
colnames(pData) <- c("condition")


########################################################################################################################################
## PERFORM DE TEST
########################################################################################################################################
Cn1 <- unique(pData$condition)[1]
Cn2 <- unique(pData$condition)[2]

C1 <- rownames(pData)[which(pData$condition == Cn1)]
C2 <- rownames(pData)[which(pData$condition == Cn2)]
print(args$filterMode)
# minimum filter mode
if(args$filterMode == 0) {
	if(args$bothConditionsSeperate == 1) {
		keep <- rowSums(exprs[, C1]>=args$keepReads) >= args$keepFraction*length(C1) & rowSums(exprs[, C2]>=args$keepReads) >= args$keepFraction*length(C2)
	}
	if(args$bothConditionsSeperate != 1) {
		keep <- rowSums(exprs[, rownames(pData)]>=args$keepReads) >= args$keepFraction*nrow(pData) 
	}
}
if(args$filterMode == 1) {
	if(args$bothConditionsSeperate == 1) {
		keep <- rowSums(exprs[, C1])>=args$keepReads*length(C1) & rowSums(exprs[, C2])>=args$keepReads*length(C2)
	}
	if(args$bothConditionsSeperate != 1) {
		keep <- rowSums(exprs[, rownames(pData)])>=args$keepReads*nrow(pData)
		#keep <- rowSums(exprs[, C1]) >= args$keepReads | rowSums(exprs[, C2]) >= args$keepReads
	}
}

# strip results 
exprs <- exprs[keep, ]

# write results
write.table(exprs, file = args$file.output, append = FALSE, quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)