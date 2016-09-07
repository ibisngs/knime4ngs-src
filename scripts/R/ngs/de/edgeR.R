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
parser <- ArgumentParser(prog="edgeR.R", description="Differential gene expression test with edgeR")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store", dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument("--countTable", type="character", action="store", dest="file.count", required=TRUE, help="path to input count file", metavar="<path>")
parser$add_argument("--annotationFile", type="character", action="store", dest="file.annotation", required=TRUE, help="path to annotation file", metavar="<path>")
parser$add_argument("--output", type="character", action="store", dest="file.output", required=TRUE, help="path to output folder", metavar="<path>")

## ARGUMENTS
parser$add_argument("--normFactors", type="character", action="store", dest="normFactors", required=TRUE, help="normalization method for size factors", metavar="<char>")
parser$add_argument("--correctPvalue", type="character", action="store", dest="correctPvalue", required=TRUE, help="correction method for p-values", metavar="<char>")

## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
loadLib("edgeR", bioC=TRUE)

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
group <- factor(pData$condition)
y <- DGEList(counts=exprs,group=group)

# get normalization factors
y <- calcNormFactors(y, method=args$normFactors)

y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)

# prepare for output
DE <- et$table
DE$adj.PValue <- p.adjust(DE$PValue, method = args$correctPvalue)
# ensure that it is sorted by adjusted p.value
DE <- DE[order(DE$PValue), ]

# filtering
all <- DE
all$ID <- rownames(all)
all <- all[, c("ID", "logFC", "logCPM", "PValue", "adj.PValue")]
colnames(all) <- c("ID", "log2FC", "aveLog2CPM", "PValue", "adj.PValue")
rownames(all) <- seq(1, nrow(all))

# write results
write.table(all, file = args$file.output, append = FALSE, quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)