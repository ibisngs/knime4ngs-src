########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="normalize.R", description="Normalize data using differenz methods")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input file", metavar="<path>")
parser$add_argument( "--output"      , type="character", action="store"     , dest="file.out"   , required=TRUE, help="path to first output file", metavar="<path>")
parser$add_argument( "--stats"       , type="character", action="store"     , dest="file.stats" , required=TRUE, help="path to first output file", metavar="<path>")

## ARGUMENTS
parser$add_argument("-c","--cols"    , type="character", action="store"     , dest="columns"    ,                 help="define the columns which shall be normalized; default:all"  , metavar="<colnames>")
parser$add_argument("-m","--method"  , type="character", action="store"     , dest="method"     ,                 help="method used for normalization"  , metavar="<'quantile normalize'>")


## parse
#print(commandArgs(trailingOnly=TRUE))
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
#args <- parser$parse_args(args.in)

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
loadLib("GenABEL")

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
	st = shapiro.test(data[, c])
	stats[c, "normality.before"] = st$p.value
}


##############################################################################################################
## Normalize
##############################################################################################################
if(args$method == "quantile normalize"){
	for(c in args$columns){
		data[, c] = rntransform(data[, c])
	}
}



##############################################################################################################
## check normality againZ
##############################################################################################################
for(c in args$columns){
	st = shapiro.test(data[, c])
	stats[ c, "normality.after"] = st$p.value
}


##############################################################################################################
## write output
##############################################################################################################

write.csv3(data, args$file.out)
write.csv3(stats, args$file.stats)




