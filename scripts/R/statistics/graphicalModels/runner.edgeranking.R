########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="graphicalModels.R", description="This Rscript creates graphical models from datamatrices")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file"        , metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--data"        , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input data file (cols=variables, rows=observations)", metavar="<path>")
parser$add_argument( "--output"      , type="character", action="store"     , dest="file.out1"  , required=TRUE, help="path to first output file"  , metavar="<path>")
parser$add_argument( "--output2"     , type="character", action="store"     , dest="file.out2"  , required=TRUE, help="path to second output file" , metavar="<path>")

## ARGUMENTS
parser$add_argument("-c", "--classes" , type="character", action="store"     , dest="classes"     , help="comma-separated list of variable classes (one entry per variable) (default: class)" , metavar="<class.1, class.2, class.3, ..., class.ncol>")
parser$add_argument("-r", "--ranker"  , type="character", action="store"     , dest="ranker"       ,default="integer=randomForest,numeric=randomForest,factor=randomForest", help="define ranker methods used for each variable-class (as class1=method1,class2=method2,...)" , metavar="<class1=method1,class2=method2,...>")
parser$add_argument("-p", "--param"   , type="character", action="store"     , dest="ranker.params",default="integer:;numeric:;factor:", help="define ranker methods used for each variable-class (as class1=method1,class2=method2,...)", metavar="<method1:param1=value1,param2=value2; method2:param1=value1,param2=value2,...>")

parser$add_argument("-s","--sampleNum", type="integer"  , action='store'     , dest='sampleNum', default='100' ,help='sum the integers (default: find the max)')
parser$add_argument("-t","--rankType" , type="character", action='store'     , dest='rankType' , default='local' ,help='shall local or global ranking be applied (global only if just one type of rankers is used!)', metavar="<local|global>")
 

## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

print(args)
# q()

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
## CRAN
loadLib("plyr")
#loadLib("igraph")
#loadLib("reshape2")
loadLib("kimisc")
CWD = normalizePath(if(is.null(thisfile())){getwd()}else{dirname(thisfile())})
source(paste(CWD, "/mixedGraphicalModels.R", sep=""))


##########################################################################################################################################
## READ DATA
##########################################################################################################################################
data <- read.csv3(args$file.in)

## check variable classes
if(is.null(args$classes)){
	var.classes = apply(data, MARGIN=2, FUN=class)
}else{
	var.classes = unlist(strsplit(args$classes, split="\\s*,\\s*", perl=T))
	if(length(var.classes) != ncol(data)){
		stop(paste("number of variable classes (", length(var.classes) ,") has to match number of variables (", ncol(data), ") (columns in data file)", sep=""))
	}
	names(var.classes) = colnames(data)
}


## check rankers
ranker = strsplit.named(args$ranker)
invalid.classes = var.classes[ ! var.classes %in% names(ranker)]
if(length(invalid.classes)!=0){
	stop(paste("no ranker given for the following classes: ", paste(unique(invalid.classes), collapse=","), sep=""))
}

## check ranker params
ranker.params = strsplit.named(args$ranker.params, split.1="\\s*;\\s*", split.2="\\s*:\\s*")
ranker.params = llply(ranker.params, .fun=strsplit.named, split.1="\\s*,\\s*", split.2="\\s*=\\s*")
invalid.classes = var.classes[ ! var.classes %in% names(ranker.params)]
if(length(invalid.classes)!=0){
	stop(paste("no ranker parameter given for the following classes: ", paste(unique(invalid.classes), collapse=","), sep=""))
}


##########################################################################################################################################
## CREATE MODEL
##########################################################################################################################################
model <- mixedGraficalModels(data, #
                             ranker        = ranker,
                             ranker.params = ranker.params,
                             var.classes,
                             stabSel.sampleNum = args$sampleNum, 
                             #stabSel.sampleSize , #
                             parallel.stabSel = TRUE)
                             
#print(names(model))


##########################################################################################################################################
## OUTPUT RESULTS MODEL
##########################################################################################################################################

write.csv3(model$ranks, args$file.out1)

metaData = t(data.frame(variables=paste(model$variables, collapse=","), stabSel.Sample.num=model$stabSel.sampleNum, n=model$n, p=model$p, time=model$run.time))
colnames(metaData) = "value"

write.csv3(metaData, args$file.out2)
