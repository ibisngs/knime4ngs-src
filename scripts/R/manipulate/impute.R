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
parser$add_argument("-m","--method"  , type="character", action="store"     , dest="method"     , required=TRUE , help="define which method should be used for imputation", metavar="<random, min, max, mean, mice, knn, SVD, lm>")
parser$add_argument("-c","--cols"    , type="character", action="store"     , dest="columns"    ,                 help="define the columns for which imputation shall be performed; default:all"  , metavar="<colnames>")

# knn imputation
parser$add_argument("--knn"          , type="integer"  , action="store"     , dest="knn"        ,                 help="[KNN-Imputation] the number of neighbors to use for imputation", metavar="<int>")
parser$add_argument("--dist"         , type="character", action="store"     , dest="dist"       ,                 help="[KNN-Imputation] distance measure used for determination of neighbours", metavar="<corr, euclidean, maximum, manhattan, canberra, binary, minkowski>")

# SVD imputation
parser$add_argument("--rank.k"       , type="integer"  , action="store"     , dest="rank.k"     ,                 help="[SVD-Imputation] the rank-k approximation to use for x", metavar="<int>")
parser$add_argument("--num.iters"    , type="integer"  , action="store"     , dest="num.iters"  ,                 help="[SVD-Imputation] the number of times to compute the rank-k approximation and impute the missing data", metavar="<int>")


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
		x.mean = mean(x, rm=T)
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




