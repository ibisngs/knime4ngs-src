#  Copyright (C) 2016 the Knime4NGS contributors.
#  Website: http://ibisngs.github.io/knime4ngs
#  
#  This file is part of the KNIME4NGS KNIME extension.
#  
#  The KNIME4NGS extension is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "calculate principal components of given data", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global", help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c( "--data"        ), type="character", action="store"     , dest="file.in"           , help="path to input file", metavar="<path>")
parser <- add_option(parser, c( "--output"      ), type="character", action="store"     , dest="file.out"          , help="path to first output file", metavar="<path>")
parser <- add_option(parser, c( "--rotation"    ), type="character", action="store"     , dest="file.rotation"     , help="path to first output file", metavar="<path>")
parser <- add_option(parser, c( "--varexplained"), type="character", action="store"     , dest="file.varexplained" , help="path to first output file", metavar="<path>")

## ARGUMENTS
parser <- add_option(parser, c("-c","--cols"    ), type="character", action="store"     , dest="columns"                           , help="define the columns (as comma-separated list or perl-regex) which shall be normalized; default:all"  , metavar="<colnames>")
#parser <- add_option(parser, c("-m","--method" ), type="character", action="store"     , dest="method"                            , help="method used for normalization"  , metavar="<'quantile normalize'>")
parser <- add_option(parser, c(       "--scale" )                  , action="store_true", dest="scale"              , default=FALSE, help="scale data before calculating principal components" )
parser <- add_option(parser, c(       "--center")                  , action="store_true", dest="center"             , default=FALSE, help="center data before calculating principal components" )
parser <- add_option(parser, c(     "--failOnNA")                  , action="store_true", dest="fail.on.na"         , default=FALSE, help="set flag if program shall fail if data contains missing values. Only complete cases are used otherwise." )


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
	if(length(args$columns)==1){
		args$columns = colnames(data)[grepl(args$columns, colnames(data), perl=T)]
	}
}


x = data[ , args$columns]
x = x[ complete.cases(x), ]
if(args$fail.on.na && (nrow(x) != nrow(data)) ){
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
if(!is.null(args$file.rotation)){
	write.csv3(pc$rotation, args$file.rotation)
}
if(!is.null(args$file.varexplained)){
	write.csv3(data.var.expl, args$file.varexplained)
}



