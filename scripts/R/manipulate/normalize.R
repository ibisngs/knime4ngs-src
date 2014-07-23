########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "Normalize data using differenz methods", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global", help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c( "--input"       ), type="character", action="store"     , dest="file.in"    , help="path to input file", metavar="<path>")
parser <- add_option(parser, c( "--output"      ), type="character", action="store"     , dest="file.out"   , help="path to first output file", metavar="<path>")
parser <- add_option(parser, c( "--stats"       ), type="character", action="store"     , dest="file.stats" , help="path to first output file", metavar="<path>")

## ARGUMENTS
parser <- add_option(parser, c("-c","--cols"    ), type="character", action="store"     , dest="columns"    , help="define the columns which shall be normalized; default:all"  , metavar="<colnames>")
parser <- add_option(parser, c("-m","--method"  ), type="character", action="store"     , dest="method"     , help="method used for normalization"  , metavar="<'quantile normalize'>")


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

loadLib("nortest")


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
## check normality
##############################################################################################################
stats = data.frame(row.names = args$columns, normality.before=rep(NA, times=length(args$columns)), normality.after=rep(NA, times=length(args$columns)))


for(c in args$columns){
	#st = shapiro.test(data[, c])
	ad = ad.test(data[, c])
	stats[c, "normality.before"] = ad$p.value
}


##############################################################################################################
## Normalize
##############################################################################################################
if(args$method == "rank transformation"){
	loadLib("GenABEL")
	for(c in args$columns){
		data[, c] = rntransform(data[, c])
	}
}else if(args$method == "quantile normalize"){
	loadLib("preprocessCore", bioC=T)
	data[, args$columns] = normalize.quantiles(as.matrix(data[, args$columns]))
}else if(args$method == "z-score"){
	for(c in args$columns){
		data[, c] = (data[, c]-mean(data[, c], na.rm=T))/sd(data[, c], na.rm=T)
	}
}



##############################################################################################################
## check normality again
##############################################################################################################
for(c in args$columns){
	#st = shapiro.test(data[, c])
	ad = ad.test(data[, c])
	stats[c, "normality.after"] = ad$p.value
}
stats$pgain = log10(stats$normality.before/stats$normality.after)

##############################################################################################################
## write output
##############################################################################################################
write.csv3(data, args$file.out)
write.csv3(stats, args$file.stats)




