########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
require(kimisc)
CWD = normalizePath(if(is.null(thisfile())){getwd()}else{dirname(thisfile())})
source(paste(CWD, "/../../utils/GLOBALS.R", sep=""))

## CRAN
loadLib("plyr")
#loadLib("igraph")
#loadLib("reshape2")

## rankers
loadLib("randomForest")
loadLib("glmnet")
loadLib("party")
source(paste(CWD, "/Ranker.R", sep=""))
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
                                rank.type = "local",
                                ...){

	if(parallel.stabSel){
		loadLib("foreach")
		loadLib("doMC")
		registerDoMC(cores=detectCores())
		warning("USING ", detectCores(), " cores for comuputations!")
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
			.fun = function(subset){rankEdges(X[subset,], var.classes=var.classes, ranker=ranker, ranker.params=ranker.params, edges.indices=edges.indices, rank.type=rank.type)},
			.parallel = parallel.stabSel
			)
	ranks = as.data.frame(matrix(ranks, nrow=stabSel.sampleNum))
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
getGraph <- function(edge.ranks, 
                     stabSel.sampleNumedges, #
                     stabSel.inclusionPerc=0.8, #
                     E.v=20){
	loadLib("igraph")
	
	variables = unique(unlist(strsplit(colnames(edge.ranks), "~")))
	p         = length(variables)
	n.edges    = ncol(edge.ranks)
	n.subsets  = nrow(edge.ranks)
	
	## if not given by the user, calculate stabSel.sampleNumedges to meet certain FWER (for details Fellinghauer et al. 2013)	
	if(missing(stabSel.sampleNumedges)){
		stabSel.sampleNumedges = floor(sqrt((2*stabSel.inclusionPerc-1) * E.v * p * (p-1) / 2))
		cat("Including top ", stabSel.sampleNumedges, "edges from each sample")
	}

	
	## per stability selection sample choose x best edges ## TODO HIGHER VALUE MEANS HIGHER PROBABILITY OF INCLUSION
	edges.included            = aaply(.data=edge.ranks, .margins=1, .fun=function(x){a<-rep(0, times=n.edges); a[which(rank(-x) < stabSel.sampleNumedges)]<-1;return(a)}, .expand=FALSE)

	#colnames(edges.included) = colnames(edge.ranks)
	#srownames(edges.included) = rownames(edge.ranks)
	
	## count number (and percentage) of samples in which a certain edge is included
	edges           = data.frame(do.call(rbind, strsplit(colnames(edge.ranks), split = "~")))
	colnames(edges) = c("v1", "v2")
	rownames(edges) = colnames(edge.ranks)
	edges$incl.n    = aaply(.data=edges.included, .margins=2, .fun=sum)
	edges$mean.rank = aaply(.data=edge.ranks, .margins=2, .fun = function(x){mean(x[,1])})
	edges$incl.perc = edges$incl.n/n.subsets
	edges$selected  = edges$incl.perc >= stabSel.inclusionPerc
	
	## generate adjacency matrix from selected edges
	g = graph.data.frame(edges[, c("v1", "v2", "selected")], directed=F)
	adjacency = get.adjacency(g,sparse=FALSE)
	adjacency = adjacency[variables,variables]

	return(list(edges=edges, adjacency=adjacency, graph=g))
}

