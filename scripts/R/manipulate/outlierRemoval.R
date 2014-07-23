########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "Remove outliers which deviate more than a certain treshold from mean", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global", help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c( "--input"       ), type="character", action="store"     , dest="file.in"    , help="path to input file", metavar="<path>")
parser <- add_option(parser, c( "--output"      ), type="character", action="store"     , dest="file.out"   , help="path to first output file", metavar="<path>")
parser <- add_option(parser, c( "--stats"       ), type="character", action="store"     , dest="file.stats" , help="path to first output file", metavar="<path>")

## ARGUMENTS
parser <- add_option(parser, c("-c","--cols"    ), type="character", action="store"     , dest="columns"    , help="define the columns for which imputation shall be performed; default:all"  , metavar="<colnames>")
parser <- add_option(parser, c("--sds"          ), type="integer"  , action="store"     , dest="sds"        , help="[Dev From Mean] maximal number of standard deviations a value can differ from mean", metavar="<int>")


## parse
args = parse_args(parser, args = commandArgs(trailingOnly = TRUE), print_help_and_exit = TRUE, positional_arguments = FALSE)

## mandatory args
if(is.null(args$file.global)){
	print_help(parser)
	warning("mandatory globals file (--globals) missing!")
	q(status=-1)
}
if(is.null(args$file.in)){
	print_help(parser)
	warning("mandatory input file (--input) missing!")
	q(status=-1)
}
if(is.null(args$file.out)){
	print_help(parser)
	warning("mandatory output file (--output) missing!")
	q(status=-1)
}
if(is.null(args$file.stats)){
	print_help(parser)
	warning("mandatory output file (--stats) missing!")
	q(status=-1)
}
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




