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
parser <- OptionParser(usage = "usage: %prog [options]", description = "This script imputes missing data in the given data file for each variable (column)", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global", help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c( "--input"       ), type="character", action="store"     , dest="file.in"                         , help="path to input file", metavar="<path>")
parser <- add_option(parser, c( "--output"      ), type="character", action="store"     , dest="file.out"                        , help="path to first output file", metavar="<path>")

## ARGUMENTS
parser <- add_option(parser, c("-m","--method"  ), type="character", action="store"     , dest="method"     , default = "random" , help="define which method should be used for imputation (default:random)", metavar="<random, min, max, mean, mice, knn, SVD, lm>")
parser <- add_option(parser, c("-c","--cols"    ), type="character", action="store"     , dest="columns"                         , help="define the columns for which imputation shall be performed; default:all"  , metavar="<colnames>")

# knn imputation
parser <- add_option(parser, c("--knn"          ), type="integer"  , action="store"     , dest="knn"                             , help="[KNN-Imputation] the number of neighbors to use for imputation", metavar="<int>")
parser <- add_option(parser, c("--dist"         ), type="character", action="store"     , dest="dist"                            , help="[KNN-Imputation] distance measure used for determination of neighbours", metavar="<corr, euclidean, maximum, manhattan, canberra, binary, minkowski>")

# SVD imputation
parser <- add_option(parser, c("--rank.k"       ), type="integer"  , action="store"     , dest="rank.k"                          , help="[SVD-Imputation] the rank-k approximation to use for x", metavar="<int>")
parser <- add_option(parser, c("--num.iters"    ), type="integer"  , action="store"     , dest="num.iters"                       , help="[SVD-Imputation] the number of times to compute the rank-k approximation and impute the missing data", metavar="<int>")


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
varorder = colnames(data)
method   = args$method
## optional arguments
if(is.null(args$columns)){
	args$columns = colnames(data)
}else{
	args$columns = unlist(strsplit(args$columns, ","))
}
# data = data.frame(matrix(rnorm(15), ncol=3, nrow=5))
# data[,1] = as.factor(data[,3])
# data[,2] = as.factor(data[,3])
# data[,3] = as.factor(data[,3])
# rownames(data) = paste("row.", c(1:5), sep="")
# colnames(data) = paste("col.", c(1:3), sep="")
# data[2,3] <- NA
# args <- list()
# args$columns=colnames(data)

##############################################################################################################
## impute functions
##############################################################################################################
impute.simple <- function(x, method="min"){
	n = colnames(x)
	x<- x[,1]
	fun=get(method)
	imp = fun(x, na.rm=T)
	x[is.na(x)] = imp
	print(paste("[", n,"] impute using ",method,  ": ", imp, sep=""))
	return(x)
}
impute.random <- function(x, cat.cutoff=10){
	if(length(unique(x))>cat.cutoff & is.numeric(x)){
		x.mean = mean(x, na.rm=T)
		x.sd   = sd(x, na.rm=T)
		#cat("mean: ", x.mean, "  sd:", x.sd)
		x[is.na(x)] <- rnorm(sum(is.na(x)), mean=x.mean, sd=x.sd)
	}else{
		x[is.na(x)] <- sample(x[!is.na(x)], size=sum(is.na(x)), replace=T)
	}
	return(x)
}

##############################################################################################################
## impute
##############################################################################################################
if(length(args$columns)>0){
	if(method == "knn"){
		loadLib("imputation")
		dist.mat = NULL
		if(args$dist == "corr"){
			loadLib("Hmisc")
			dist.mat = rcorr(data)
		}else if(args$dist %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
			dist.mat = dist(data, method=args$dist, diag=FALSE, upper=TRUE)
		}
		imputed = kNNImpute(data, k=args$knn, x.dist=dist.mat, verbose=T) ## TODO:distance measure?
		imputed = imputed$x
		
	}else if(method == "lm"){
		loadLib("imputation")
		imputed = lmImpute(data)
	
	}else if(method == "SVD"){
		loadLib("imputation")
		imputed = SVDImpute(data, args$rank.k, num.iters=args$num.iters, verbose=T)
	
	}else if(method == "min" || method == "max" || method == "mean"){
		loadLib("plyr")
		imputed = t(aaply(.data=data[,args$columns], .margins=2, .fun=impute.simple, method=method))
		rownames(imputed) = rownames(data)
		
	}else if(method == "mice"){ ## TODO additional arguments
		loadLib("mice", bioC=T)
		imputed = mice(data)
		print(imputed$loggedEvents)
		imputed = imputed$data
		
	}else if(method == "random"){ ## TODO variable cutoff for categorical
		imputed = apply(data[,args$columns], MARGIN=2, FUN=impute.random, cat.cutoff=10)
		rownames(imputed) = rownames(data)
		colnames(imputed) = args$columns
	}else{
		warning("Method unknown")
	q()
	}
	
	## replace original data with imputed data
	data[, args$columns] = imputed[, args$columns]

}

##############################################################################################################
## write output
##############################################################################################################

data = data[, varorder]
write.csv3(data, args$file.out)




