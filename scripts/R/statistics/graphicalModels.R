########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
require(kimisc)
CWD = normalizePath(if(is.null(thisfile())){getwd()}else{dirname(thisfile())})
source(paste(CWD, "/../utils/GLOBALS.R", sep=""))

## CRAN
loadLib("plyr")
#loadLib("igraph")
#loadLib("reshape2")

## rankers
loadLib("randomForest")
loadLib("glmnet")
source(paste(CWD, "/graphicalModels_Ranker.R", sep=""))

########################################################################################################################################
## FUNCTIONS
########################################################################################################################################
# X			data matrix with columns containing variables and rows containing  (p*n)
# method.*		ranker method to be used for certain type of data
# stabSel.sampleNum	number of subsets to be drawn for stability selection
# stabSel.sampleSize	size of subsets for stability selection (by default floor(n/2))
# parallel.stabSel	execute subsamples of stability selection in parallel
mixedGraficalModels <- function(X, #
                                ranker        = list(numeric="randomForest", double="randomForest" , factor="randomForest"),
                                ranker.params = list(numeric=list()        , double=list()         , factor=list()),
                                var.classes,
                                stabSel.sampleNum = 100, 
                                stabSel.sampleSize, #
                                parallel.stabSel = TRUE,#
                                ...){

	if(parallel.stabSel){
		loadLib("foreach")
		loadLib("doMC")
		registerDoMC(cores=detectCores())
	}
	
	## measure time of the entire procedure
	time.start <- Sys.time()
    	n = nrow(X)
	p = ncol(X)
	
	## default subset size n/2
	if(missing(stabSel.sampleSize) || is.na(stabSel.sampleSize)){
		stabSel.sampleSize = floor(n/2)
	}

	if(missing(var.classes)){
		var.classes = apply(data, MARGIN=2, FUN=class)
	}
	
	## sample subsets for stability selection
	# TODO: set seed /random initialization...
	if( stabSel.sampleNum == 1){
		samples = list(c(1:n))
	}else{
		samples = llply(.data=c(1:stabSel.sampleNum),.fun=function(i){return(sample(rownames(X),stabSel.sampleSize,replace=F))})
	}
	## order of edges (from an adjacency matrix) as vector
	edges.matrix  <- matrix(1:(p*p),ncol=p, nrow=p, dimnames=list(colnames(X), colnames(X))) #, 
	edges.indices <- as.vector(edges.matrix[upper.tri(edges.matrix)])

	## rank edges in each subsample
	ranks   = laply(.data=samples, 
			.fun = function(subset){rankEdges(X[subset,], var.classes=var.classes, ranker=ranker, ranker.params=ranker.params, edges.indices=edges.indices)},
			.parallel = parallel.stabSel, .inform=F
			)
	rownames(ranks) = paste("ranks.subset", c(1:stabSel.sampleNum), sep=".")
	colnames(ranks) = make.edgesNames(colnames(X))

	## create output
	time.end <- Sys.time()
	MODEL = list(ranks=ranks, variables=colnames(X), stabSel.sampleNum=stabSel.sampleNum, n=n, p=p, run.time=as.numeric(time.end-time.start, units="secs"))
	return(invisible(MODEL))
}

##########################################################################################################################################
## FUNCTIONS
##########################################################################################################################################
# stabSel.inclusionPerc percentage of subsamples in which a edge must be present to be included in final graph
# stabSel.edges 	number of edges to choose in each subsample (leave blank if desired E.v is given)
# E.v			desired maximal number of False positive edges; only guaranteed if stabSel.edges is not given manually
getGraph <- function(model, 
                     stabSel.sampleNumedges, #
                     stabSel.inclusionPerc=0.8, #
                     E.v=20,
                     mode="adjacency"){
	## if not given by the user, calculate stabSel.sampleNumedges to meet certain FWER (for details Fellinghauer et al. 2013)
	if(missing(stabSel.sampleNumedges)){
		stabSel.sampleNumedges = floor(sqrt((2*stabSel.inclusionPerc-1) * E.v * model$p * (model$p-1) / 2))
	}
	
	edge.ranks = model$ranks
	n.edges    = ncol(edge.ranks)
	n.subsets  = nrow(edge.ranks)
	
	if(mode == "adjacency"){
		## per stability selection sample choose x best edges
		edges.included            = aaply(.data=edge.ranks, .margins=1, .fun=function(x){a<-rep(0, times=n.edges); a[which(rank(x, ties.method='min') < stabSel.sampleNumedges)]<-1;return(a)})
		
		## count number (and percentage) of samples in which a certain edge is included
		edges.included.percentage = aaply(.data=edges.included, .margins=2, .fun=sum)
		names(edges.included.percentage) = colnames(edge.ranks)
		edges.included.percentage = edges.included.percentage/n.subsets
		
		## select edges which meet frequency cutoff
		edges.selected            = edges.included.percentage[edges.included.percentage >= stabSel.inclusionPerc]
		
		## generate adjacency matrix from selected edges
		adjacency.matrix          = matrix(0, ncol=model$p, nrow=model$p)
		rownames(adjacency.matrix)= colnames(adjacency.matrix) = model$variables
		adjacency.matrix[as.numeric(names(edges.selected))] = edges.selected
		adjacency.matrix = adjacency.matrix+t(adjacency.matrix)

		return(adjacency.matrix)
	} else if( mode == "edge.ranks") {
		edges.ranks.mean          = aaply(.data=edge.ranks, .margins=1, .fun=function(x){return(rank(x, ties.method='min'))})
		edges.ranks.mean          = aaply(.data=edges.ranks.mean, .margins=2, .fun=mean)
		#edges.ranks.mean          = rank(edges.ranks.mean)
		return(edges.ranks.mean)
	}
}

