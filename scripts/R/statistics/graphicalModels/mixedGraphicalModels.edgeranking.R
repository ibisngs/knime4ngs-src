########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "This script samples random subsets from the datamatrix and ranks all possible edges according to the data", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global", help="path to globals file"        , metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c("--data"          ), type="character", action="store"     , dest="file.in"    , help="path to input data file (cols=variables, rows=observations)", metavar="<path>")
parser <- add_option(parser, c("--output"        ), type="character", action="store"     , dest="file.out1"  , help="path to first output file"  , metavar="<path>")
parser <- add_option(parser, c("--output2"       ), type="character", action="store"     , dest="file.out2"  , help="path to second output file" , metavar="<path>")

## ARGUMENTS
parser <- add_option(parser, c("-c", "--fclasses"), type="character", action="store"     , dest="file.classes" , help="file which contains the same columns as the input data but only one row containing the datatype of the according variable (--classes can be used instead)" , metavar="<path>")
parser <- add_option(parser, c("-fc", "--classes"), type="character", action="store"     , dest="classes"      , help="comma-separated list of variable classes (one entry per variable) (default: class) (--fclasses can be used instead)" , metavar="<class.1,class.2,...>")
parser <- add_option(parser, c("-r", "--ranker"  ), type="character", action="store"     , dest="ranker"       , default="integer=randomForest,numeric=randomForest,factor=randomForest", help="define ranker methods used for each variable-class" , metavar="<class1=method1,class2=method2,...>")
parser <- add_option(parser, c("-p", "--param"   ), type="character", action="store"     , dest="ranker.params", default="integer:;numeric:;factor:", help="define ranker methods used for each variable-class", metavar="<method1:param1=value1,param2=value2;method2:param1=value1,param2=value2,...>")

parser <- add_option(parser, c("-s","--sampleNum" ), type="integer"  , action='store'     , dest='sampleNum'    , default='100'  , help='sum the integers (default: find the max)')
parser <- add_option(parser, c("-z","--sampleSize"), type="double"   , action='store'     , dest='samplesize'   , default='0.8'  , help='size of each subsample (default:0.8)', metavar='<double>')
parser <- add_option(parser, c("-t","--ranktype"  ), type="character", action='store'     , dest='rankType'     , default='local', help='shall local or global ranking be applied (global only if just one type of rankers is used!)', metavar="<local|global>")
 
parser <- add_option(parser, c("-rs","--rseed"    ), type="integer"  , action='store'     , dest='rseed'        ,                  help='random seed for reproducible random samples', metavar="<int>")
parser <- add_option(parser, c("-pc","--cores"    ), type="integer"  , action='store'     , dest='parallel'     , default=1      , help='number of parallel threads used for calculate edgeranking (by default number of cores)', metavar="<int>")



args = parse_args(parser, args = commandArgs(trailingOnly = TRUE), print_help_and_exit = TRUE, positional_arguments = FALSE)

## mandatory args
if(is.null(args$file.global)){
	print_help(parser)
	warning("mandatory globals file (--globals) missing!")
	q(status=-1)
}
if(is.null(args$file.in)){
	print_help(parser)
	warning("mandatory input file (--data) missing!")
	q(status=-1)
}
if(is.null(args$file.out1)){
	print_help(parser)
	warning("mandatory output file (--output) missing!")
	q(status=-1)
}



########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
## CRAN
loadLib("plyr")
CWD = getCWD()
source(paste(CWD, "/mixedGraphicalModels.R", sep=""))


##########################################################################################################################################
## READ DATA
##########################################################################################################################################
data <- read.csv3(args$file.in)


## check variable classes
if(!is.null(args$file.classes)){
	var.classes <- unlist(as.list(read.csv3(args$file.classes, rownames=F, stringsAsFactors=F)))
}else if(!is.null(args$classes)){
	var.classes = unlist(strsplit(args$classes, split="\\s*,\\s*", perl=T))
	names(var.classes) = colnames(data)
}else{
	var.classes = apply(data, MARGIN=2, FUN=class)
}
if(length(var.classes) != ncol(data)){
	stop(paste("number of variable classes (", length(var.classes) ,") has to match number of variables (", ncol(data), ") (columns in data file)", sep=""))
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

## init random generator
if(is.null(args$rseed)){
	args$rseed = sample(.Machine$integer.max, 1)
}
set.seed(args$rseed)

# data = data[1:20,1:3]
# var.classes = var.classes[1:3]
# args$sampleNum=10
#args$parallel=1
##########################################################################################################################################
## CREATE MODEL
##########################################################################################################################################
model <- mixedGraficalModels(data, #
                             ranker        = ranker,
                             ranker.params = ranker.params,
                             var.classes   = var.classes,
                             stabSel.sampleNum = args$sampleNum, 
                             stabSel.sampleSize = floor(args$samplesize*nrow(data)),
                             rankType = args$rankType,
                             threads = args$parallel)
                             
#print(names(model))


##########################################################################################################################################
## OUTPUT RESULTS MODEL
##########################################################################################################################################
write.csv3(model$ranks, args$file.out1)


if(!is.null(args$file.out2)){
	metaData = data.frame(variables=paste(model$variables, collapse=","), stabSel.Sample.num=model$stabSel.sampleNum, n=model$n, p=model$p, seed=args$rseed, time=model$run.time)
	write.csv3(metaData, args$file.out2)
}

