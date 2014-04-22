########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="impute.R", description="This imputes missing data in the given data file for each variable (column)")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input file", metavar="<path>")
parser$add_argument( "--output"      , type="character", action="store"     , dest="file.out"   , required=TRUE, help="path to first output file", metavar="<path>")

## ARGUMENTS
parser$add_argument("-p","--pvalcol" , type="character", action="store"     , dest="col.pval"   ,                help="name of the column which contains the p-values"  , metavar="<colname>")
parser$add_argument("-o","--outcol"  , type="character", action="store"     , dest="col.out"    ,                help="name of the column which contains the corrected p-values"  , metavar="<name>")

parser$add_argument("-m","--method"  , type="character", action="store"     , dest="method"     , required=TRUE, help="define which method should be used for p-value correction", metavar="<holm, hochberg, hommel, bonferroni, BH, BY, fdr, none>")
parser$add_argument("--n"            , type="integer"  , action="store"     , dest="n"          ,                help="number of tests correct for (by default number of rows)", metavar="<int>")

parser$add_argument( "--format"      , type="character", action="store"     , dest="format"     , required=TRUE, help="output format append corrected p-values (a), replace p-value column (r) or return corrected p-values only (o)", metavar="<a,r,o>")

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

##############################################################################################################
## adjust P-values
##############################################################################################################
corrected.pvals = p.adjust(data[, args$col.pval], method=args$method, n=args$n)

if(args$format == "a"){
	data[, args$col.out] = corrected.pvals
}else if(args$format == "r"){
	id = which(colnames(data) == args$col.pval)
	data[, id] = corrected.pvals
	colnames(data)[id] = args$col.out
}else if(args$format == "o"){
	data[, args$col.out] = corrected.pvals
	data = data[ , args$col.out, drop=F]
} else{
	stop("Don't know what to return! Format must be one of 'a', 'r', 'o'")
}

##############################################################################################################
## write output
##############################################################################################################
write.csv3(data, args$file.out)




