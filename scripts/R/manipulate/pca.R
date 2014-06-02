########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="pca.R", description="calculate principal components of given data")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input file", metavar="<path>")
parser$add_argument( "--output"      , type="character", action="store"     , dest="file.out"   , required=TRUE, help="path to first output file", metavar="<path>")
parser$add_argument( "--rotation"       , type="character", action="store"     , dest="file.rotation" , required=TRUE, help="path to first output file", metavar="<path>")
parser$add_argument( "--varexplained"       , type="character", action="store"     , dest="file.varexplained" , required=TRUE, help="path to first output file", metavar="<path>")

## ARGUMENTS
parser$add_argument("-c","--cols"    , type="character", action="store"     , dest="columns"    ,                 help="define the columns which shall be normalized; default:all"  , metavar="<colnames>")
#parser$add_argument("-m","--method"  , type="character", action="store"     , dest="method"     ,                 help="method used for normalization"  , metavar="<'quantile normalize'>")
parser$add_argument(       "--scale" , action="store_true",  dest="scale"   , help="scale data before calculating principal components" )
parser$add_argument(       "--center" , action="store_true",  dest="center" , help="center data before calculating principal components" )
parser$add_argument(     "--failOnNA" , action="store_true",  dest="fail.on.na" , help="set flag if program shall fail if data contains missing values. Only complete cases are used otherwise." )


## parse
#print(commandArgs(trailingOnly=TRUE))
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
#args <- parser$parse_args(args.in)

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)


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

x = data[ , args$columns]
x = x[ complete.cases(x), ]
if(args$fail.on.na && nrow(x) != nrow(data)){
	cat("Data contains missing values!")
	stop("Data contains missing values!", call.=F)

}
##############################################################################################################
## PCA
##############################################################################################################
pc = prcomp(x, center=args$center, scale=args$scale)


##############################################################################################################
## VARIANCE EXPLAINED
##############################################################################################################
lambda.all <- pc$sdev^2
var.expl.cum = rep(NA, length(lambda.all))
var.expl     = rep(NA, length(lambda.all))
for(i in 1:length(lambda.all)){
	var.expl.cum[i] = cumsum(lambda.all)[i]/sum(lambda.all)
	var.expl[i]     = lambda.all[i]/sum(lambda.all)
}
data.var.expl = data.frame(num=1:length(lambda.all), lambda=lambda.all, var.expl=var.expl, var.expl.cum=var.expl.cum)
rownames(data.var.expl) = colnames(pc$x)


##############################################################################################################
## write output
##############################################################################################################
output = data[ , !colnames(data) %in% args$columns ]
output[ , colnames(pc$x)] = NA
output[ rownames(pc$x), colnames(pc$x)] = pc$x

write.csv3(output, args$file.out)
write.csv3(pc$rotation, args$file.rotation)
write.csv3(data.var.expl, args$file.varexplained)




