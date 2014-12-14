########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="limma.R", description="Differential gene expression test with limma")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store", dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument("--countTable", type="character", action="store", dest="file.count", required=TRUE, help="path to input count file", metavar="<path>")
parser$add_argument("--annotationFile", type="character", action="store", dest="file.annotation", required=TRUE, help="path to annotation file", metavar="<path>")
parser$add_argument("--output", type="character", action="store", dest="file.output", required=TRUE, help="path to output folder", metavar="<path>")

## ARGUMENTS
parser$add_argument("--normFactors", type="character", action="store", dest="normFactors", required=TRUE, help="normalization method for size factors", metavar="<char>")
parser$add_argument("--normCPM", type="character", action="store", dest="normCPM", required=TRUE, help="normalization method for CPM values", metavar="<char>")
parser$add_argument("--correctPvalue", type="character", action="store", dest="correctPvalue", required=TRUE, help="correction method for p-values", metavar="<char>")

## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
loadLib("Biobase", bioC=TRUE)
loadLib("limma", bioC=TRUE)
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
phenoData <- new("AnnotatedDataFrame", data=pData) 
eset <- ExpressionSet(assayData=exprs, phenoData=phenoData)

# normalize counts
nf <- calcNormFactors(eset, method=args$normFactors)
groups <- phenoData(eset)$condition
design <- model.matrix(~ groups)
y <- voom(exprs(eset), design, lib.size=colSums(exprs(eset))*nf, normalize.method=args$normCPM)

# build linear model
fit <- lmFit(y,design)
fit <- eBayes(fit)

# ensure that it is sorted by adjusted p.value
DE <- topTable(fit, coef=2, number=Inf, adjust.method=args$correctPvalue)
DE <- DE[order(DE$adj.P.Val), ]

all <- DE
all$ID <- rownames(all)
all <- all[, c("ID", "logFC", "AveExpr", "B", "t", "P.Value", "adj.P.Val")]
colnames(all) <- c("ID", "log2FC", "aveLog2CPM", "B", "t", "PValue", "adj.PValue")
all$log2FC <- -all$log2FC # invert to get HF/LF
all$t <- -all$t # invert to get HF/LF
rownames(all) <- seq(1, nrow(all))

# write results
write.table(all, file = args$file.output, append = FALSE, quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)