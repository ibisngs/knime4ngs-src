########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "This script adjusts p-values for multiple testing", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global"              , help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c( "--input"       ), type="character", action="store"     , dest="file.in"                  , help="path to input file", metavar="<path>")
parser <- add_option(parser, c( "--output"      ), type="character", action="store"     , dest="file.out"                 , help="path to first output file", metavar="<path>")

## ARGUMENTS
parser <- add_option(parser, c("-p","--pvalcol" ), type="character", action="store"     , dest="col.pval"                 , help="name of the column which contains the p-values"  , metavar="<colname>")
parser <- add_option(parser, c("-o","--outcol"  ), type="character", action="store"     , dest="col.out"                  , help="name of the column which contains the corrected p-values"  , metavar="<name>")

parser <- add_option(parser, c("-m","--method"  ), type="character", action="store"     , dest="method"     , default="BH", help="define which method should be used for p-value correction", metavar="<holm, hochberg, hommel, bonferroni, BH, BY, fdr, none>")
parser <- add_option(parser, c("--n"            ), type="integer"  , action="store"     , dest="n"                        , help="number of tests correct for (by default number of rows)", metavar="<int>")

parser <- add_option(parser, c( "--format"      ), type="character", action="store"     , dest="format"     , default="a" , help="output format append corrected p-values (a), replace p-value column (r) or return corrected p-values only (o)", metavar="<a,r,o>")

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
if(is.null(args$file.in)){
	print_help(parser)
	warning("mandatory output file (--output) missing!")
	q(status=-1)
}


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




