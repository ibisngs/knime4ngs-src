########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="gaussianGraphicalModels.R", description="The Rscript reads a datamatrix containing measured variables (columns) for several samples (rows) and calculates GGMs from the data.")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file"        , metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--data"          , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input data file (cols=variables, rows=observations)", metavar="<path>")
parser$add_argument( "--output"        , type="character", action="store"     , dest="file.edges" , required=TRUE, help="path to output file (edges list)" , metavar="<path>")

# ## ARGUMENTS
# parser$add_argument("-c", "--classes"  , type="character", action="store"     , dest="classes"      , help="comma-separated list of variable classes (one entry per variable) (default: class)" , metavar="<class.1, class.2, class.3, ..., class.ncol>")
# parser$add_argument("-r", "--ranker"   , type="character", action="store"     , dest="ranker"       , default="integer=randomForest,numeric=randomForest,factor=randomForest", help="define ranker methods used for each variable-class (as class1=method1,class2=method2,...)" , metavar="<class1=method1,class2=method2,...>")
# parser$add_argument("-p", "--param"    , type="character", action="store"     , dest="ranker.params", default="integer:;numeric:;factor:", help="define ranker methods used for each variable-class (as class1=method1,class2=method2,...)", metavar="<method1:param1=value1,param2=value2; method2:param1=value1,param2=value2,...>")
# 
# parser$add_argument("-s","--sampleNum" , type="integer"  , action='store'     , dest='sampleNum'    , default='100'  , help='sum the integers (default: find the max)')
# parser$add_argument("-z","--sampleSize", type="double"   , action='store'     , dest='samplesize'   , default='0.8'    , help='size of each subsample (default:0.8)', metavar='<double>')
# parser$add_argument("-t","--rankType"  , type="character", action='store'     , dest='rankType'     , default='local', help='shall local or global ranking be applied (global only if just one type of rankers is used!)', metavar="<local|global>")
#  
# parser$add_argument("-rs","--rseed"    , type="integer"  , action='store'     , dest='rseed'        ,                  help='random seed for reproducible random samples', metavar="<int>")
# parser$add_argument("-pc","--cores"    , type="integer"  , action='store'     , dest='parallel'     ,                  help='number of parallel threads used for calculate edgeranking (by default number of cores)', metavar="<int>")
  
## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

print(args)
# q()

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
## CRAN
loadLib("GeneNet")

##########################################################################################################################################
## READ DATA
##########################################################################################################################################
data <- read.csv3(args$file.in)

##########################################################################################################################################
## CALCULATE GAUSSIAN GRAPHICAL MODEL
##########################################################################################################################################
ggm   <- ggm.estimate.pcor(data)
edges <- ggm.test.edges(ggm, fdr=T, plot=F, direct=F)

## add node names instead of numbers
nodenames = colnames(data)
for(i in 1:nrow(edges)){
	edges[i, "node1"] = nodenames[as.numeric(edges[i, "node1"])]
	edges[i, "node2"] = nodenames[as.numeric(edges[i, "node2"])]
}

## order cols
edges = edges[, c("node1", "node2", "pcor", "pval", "qval", "prob")]
colnames(edges) = c("v1", "v2", "pcor", "p", "q", "prob")
rownames(edges) = paste(edges$v1, edges$v2, sep="~")


##########################################################################################################################################
## WRITE DATA
##########################################################################################################################################
write.csv3(edges, args$file.edges)
