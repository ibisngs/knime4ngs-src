# ########################################################################################################################################
# ## LOAD ARGUMENTS
# ########################################################################################################################################
# require(argparse)
# parser <- ArgumentParser(prog="mixedGraphicalModels.R", description="This script is not intendet to be run! Use the wrapper scripts!")
# 
# ## GLOBALS 
# parser$add_argument("-g", "--globals", type="character", action="store", dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")
# 
# ## parse
# args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
# 
# ########################################################################################################################################
# ## LOAD LIBRARIES
# ########################################################################################################################################
# source(args$file.glob)

## CRAN
loadLib("plyr")
## rankers
loadLib("randomForest")
loadLib("glmnet")
loadLib("party")
rankerpath = paste(getCWD(), "/Ranker.R", sep="")
source(rankerpath)

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
                                threads = NULL,#
                                rank.type = "local",
                                ...){

	## parallelization
	parallel.stabSel = FALSE
        if(is.null(threads) | threads<1){
		loadLib("doMC")
		threads = detectCores()
        }
	if(threads > 1){
		loadLib("foreach")
		loadLib("doMC")
		registerDoMC(cores=threads)
		warning("USING ", threads, " threads for computations!")
		parallel.stabSel = TRUE
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
	samples = llply(.data=c(1:stabSel.sampleNum),.fun=function(i){return(sample(rownames(X),stabSel.sampleSize,replace=F))})

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
getEdgesList <- function(edge.ranks, 
                         stabSel.sampleNumedges, #
                         stabSel.inclusionPerc=0.8){
	
	variables = unique(unlist(strsplit(colnames(edge.ranks), "~")))
	p         = length(variables)
	n.edges    = ncol(edge.ranks)
	n.subsets  = nrow(edge.ranks)

		
	## per stability selection sample choose x best edges ## TODO HIGHER VALUE MEANS HIGHER PROBABILITY OF INCLUSION
	edges.included            = aaply(.data=edge.ranks, .margins=1, .fun=function(x){a<-rep(0, times=n.edges); a[which(rank(-x) < stabSel.sampleNumedges)]<-1;return(a)}, .expand=FALSE)
	#colnames(edges.included) = colnames(edge.ranks)
	#srownames(edges.included)= rownames(edge.ranks)
	

	## count number (and percentage) of samples in which a certain edge is included
	edges           = data.frame(do.call(rbind, strsplit(colnames(edge.ranks), split = "~")))
	colnames(edges) = c("v1", "v2")
	rownames(edges) = colnames(edge.ranks)
	edges$incl.n    = aaply(.data=edges.included, .margins=2, .fun=sum)
	edges$mean.rank = aaply(.data=edge.ranks, .margins=2, .fun = function(x){mean(x[,1])})
	edges$incl.perc = edges$incl.n/n.subsets

	return(edges)
}


getAdjacency <- function(edges, variables, v1="v1", v2="v2", feature="selected"){
	loadLib("igraph")
	g = graph.data.frame(edges[, c(v1, v2, feature)], directed=F)
	adjacency = get.adjacency(g,sparse=FALSE, attr=feature)
	adjacency = adjacency[variables,variables]
}


# stabSel.inclusionPerc percentage of subsamples in which a edge must be present to be included in final graph
# stabSel.edges 	number of edges to choose in each subsample (leave blank if desired E.v is given)
# E.v			desired maximal number of False positive edges; only guaranteed if stabSel.edges is not given manually
getGraph.FWER <- function(edge.ranks, 
                          stabSel.inclusionPerc=0.8, #
                          E.v=20){

	## if not given by the user, calculate stabSel.sampleNumedges to meet certain FWER (for details Fellinghauer et al. 2013)	
	variables = unique(unlist(strsplit(colnames(edge.ranks), "~")))
	p         = length(variables)
	stabSel.sampleNumedges = floor(sqrt((2*stabSel.inclusionPerc-1) * E.v * p * (p-1) / 2))
	#cat("Including top ", stabSel.sampleNumedges, "edges from each sample")

	edges = getEdgesList(edge.ranks, stabSel.sampleNumedges=stabSel.sampleNumedges, stabSel.inclusionPerc=stabSel.inclusionPerc)
	edges$selected  = edges$incl.perc >= stabSel.inclusionPerc
	
	## generate adjacency matrix from selected edges
	adjacency = getAdjacency(edges, variables)

	## filter edgeslist
	edges = edges[edges$selected, ]
	edges = edges[, colnames(edges) != "selected"]
	
	return(list(edges=edges, adjacency=adjacency))
}



# stabSel.sampleNumedges: number of edges to be selected in each subset
# stabSel.inclusionPerc percentage of subsamples in which a edge must be present to be included in final graph
getGraph.empiricalPvalues <- function(edge.ranks, 
                                      edge.ranks.background,
                                      stabSel.sampleNumedges,
                                      stabSel.inclusionPerc=0.8){

	variables = unique(unlist(strsplit(colnames(edge.ranks), "~")))
	p         = length(variables)
	n.edges    = ncol(edge.ranks)
	n.subsets  = nrow(edge.ranks)
	
	## BACKGROUND DISTRIBUTION
	ep.samples =  sub(".*\\.", "", rownames(edge.ranks.background))
	background.samples = split(edge.ranks.background, ep.samples)
	background.inclusionPerc = data.frame()
	background.inclusionNum  = data.frame()
	for(s in 1:length(background.samples)){
		edges <- getEdgesList(background.samples[[s]], stabSel.sampleNumedges=stabSel.sampleNumedges, stabSel.inclusionPerc=stabSel.inclusionPerc)
		background.inclusionPerc = rbind(background.inclusionPerc, edges$incl.perc)
	}
	colnames(background.inclusionPerc) <- colnames(background.samples[[1]])

	## COMPARE TO REAL DATA
	edges <- getEdgesList(edge.ranks, stabSel.sampleNumedges=stabSel.sampleNumedges, stabSel.inclusionPerc=stabSel.inclusionPerc)
	edges$p.value = NA
	for(e in rownames(edges)){
		## compare inclusion percentage with background data
		edges[e, "p.value"] = pnorm(edges[e, "incl.perc"], mean=mean(background.inclusionPerc[, e]), sd=sd(background.inclusionPerc[, e]), lower.tail=F)
	}
	
	## generate adjacency matrix from selected edges
	adjacency = getAdjacency(edges, variables, feature="p.value")
	
	return(list(edges=edges, adjacency=adjacency))
}

