########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="outlierRemoval.R", description="Remove outliers which deviate more than a certain treshold from mean")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input file", metavar="<path>")
parser$add_argument( "--output"      , type="character", action="store"     , dest="file.out"   , required=TRUE, help="path to first output file", metavar="<path>")
parser$add_argument( "--stats"       , type="character", action="store"     , dest="file.stats" , required=TRUE, help="path to first output file", metavar="<path>")

## ARGUMENTS
parser$add_argument("-c","--cols"    , type="character", action="store"     , dest="columns"    ,                 help="define the columns for which imputation shall be performed; default:all"  , metavar="<colnames>")
parser$add_argument("--sds"          , type="integer"  , action="store"     , dest="sds"        ,                 help="[Dev From Mean] maximal number of standard deviations a value can differ from mean", metavar="<int>")


## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))


########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
loadLib("plyr")


##############################################################################################################
## READ DATA
##############################################################################################################
data <- read.csv3(args$file.in)


## optional arguments
if(is.null(args$columns)){
	args$columns = colnames(data)
}else{
	args$columns = unlist(strsplit(args$columns, ","))
}


##############################################################################################################
## impute functions
##############################################################################################################
dev.from.mean <- function(xn){
	x = data[, xn]
	m = mean(x, na.rm=T)
	s = sd(x, na.rm=T)
	to.remove = which(!is.na(x) & abs(x-m)/s> args$sds)
	x[ to.remove ] = NA
	return(list(data=x, stats=length(to.remove)))
}


##############################################################################################################
## remove outliers
##############################################################################################################
stats = list()
for(c in args$columns){
	tmp = dev.from.mean(c)
	data[, c] = tmp$data
	stats[[c]] = tmp$stats
}
stats = data.frame(unlist(stats))
colnames(stats) = c("removed.measurements")
stats$removed.measurements.perc = stats$removed.measurements/nrow(data)
##############################################################################################################
## write output
##############################################################################################################
write.csv3(data, args$file.out)
write.csv3(stats, args$file.stats)




