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
                                pairedSamples = FALSE,
                                ...){

	## parallelization
	parallel.stabSel = FALSE
        if(is.null(threads) | threads<1){
		loadLib("doMC")
		threads = detectCores()-1
        }
	if(threads > 1){
		loadLib("doMC")
		registerDoMC(cores=threads)
		warning("USING ", threads, " threads for computations!")
		parallel.stabSel = TRUE
	}
	cat("USING ", threads, " threads for computations!\n")
	
	## measure time of the entire procedure
	time.start <- Sys.time()
    	n = nrow(X)
	p = ncol(X)
	
	## default subset size n/2
	if(missing(stabSel.sampleSize) || is.na(stabSel.sampleSize)){
		stabSel.sampleSize = floor(n/2)
	}

	if(missing(var.classes)){
		var.classes = apply(X, MARGIN=2, FUN=class)
	}
	
	## sample subsets for stability selection
	sampleNum = stabSel.sampleNum
	if(pairedSamples){
		if(stabSel.sampleNum %% 2 != 0){
			warning("Paired subsample drawing only possible for even sample numbers!")
			q(status=-1)
		}
		sampleNum = stabSel.sampleNum/2
	}
	samples = unlist(llply(.data=c(1:sampleNum),.fun=drawSamples, X, stabSel.sampleSize=stabSel.sampleSize, paired=T), recursive=F)

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
	
	cat("Finished after ", as.numeric(time.end-time.start, units="hours"), " hours\n")
	return(invisible(MODEL))
}

drawSamples = function(i, X, stabSel.sampleSize=0.5 , paired = F){
	
	if(!paired){
		return(list(sample(rownames(X),stabSel.sampleSize,replace=F)))
	}else{
		sample.1 = sample(rownames(X),floor(nrow(X)/2),replace=F)
		sample.2 = rownames(X)[!rownames(X) %in% sample.1]
		# if uneven sample size -> remove one from sample.2
		if(length(sample.2) > length(sample.1)){
			sample.2 = sample.2[-sample(c(1:length(sample.2)), 1)]
		}
		return(list(sample.1,sample.2))
	}
}
##########################################################################################################################################
## FUNCTIONS
##########################################################################################################################################
# stabSel.inclusionPerc percentage of subsamples in which a edge must be present to be included in final graph
# stabSel.edges 	number of edges to choose in each subsample (leave blank if desired E.v is given)
getEdgesList <- function(edge.ranks, 
                         stabSel.sampleNumedges, #
                         stabSel.inclusionPerc = 0.8,
                         parallel = FALSE){

	n.edges    = ncol(edge.ranks)
	n.subsets  = nrow(edge.ranks)

	## per stability selection sample choose x best edges ## TODO HIGHER VALUE MEANS HIGHER PROBABILITY OF INCLUSION
	edges.included            = aaply(.data=edge.ranks, .margins=1, .fun=function(x){a<-rep(0, times=n.edges); a[which(rank(-x) < stabSel.sampleNumedges)]<-1;return(a)}, .expand=FALSE, .parallel=parallel)
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
	adjacency = 1*adjacency # to numeric
	
	return(adjacency)
}


# stabSel.inclusionPerc percentage of subsamples in which a edge must be present to be included in final graph
# stabSel.edges 	number of edges to choose in each subsample (leave blank if desired E.v is given)
# E.v			desired maximal number of False positive edges; only guaranteed if stabSel.edges is not given manually
getGraph.FWER <- function(edge.ranks, 
                          stabSel.inclusionPerc=0.8, #
                          E.v=20,
                          parallel = FALSE){

	## if not given by the user, calculate stabSel.sampleNumedges to meet certain FWER (for details Fellinghauer et al. 2013)	
	variables = unique(unlist(strsplit(colnames(edge.ranks), "~")))
	p         = length(variables)
	stabSel.sampleNumedges = floor(sqrt((2*stabSel.inclusionPerc-1) * E.v * p * (p-1) / 2))
	#cat("Including top ", stabSel.sampleNumedges, "edges from each sample")

	edges = getEdgesList(edge.ranks, stabSel.sampleNumedges=stabSel.sampleNumedges, stabSel.inclusionPerc=stabSel.inclusionPerc, parallel=parallel)
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
                                      stabSel.inclusionPerc=0.8,
                                      parallel=parallel){

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
		edges <- getEdgesList(background.samples[[s]], stabSel.sampleNumedges=stabSel.sampleNumedges, stabSel.inclusionPerc=stabSel.inclusionPerc, parallel=parallel)
		background.inclusionPerc = rbind(background.inclusionPerc, edges$incl.perc)
	}
	colnames(background.inclusionPerc) <- colnames(background.samples[[1]])

	## COMPARE TO REAL DATA
	edges <- getEdgesList(edge.ranks, stabSel.sampleNumedges=stabSel.sampleNumedges, stabSel.inclusionPerc=stabSel.inclusionPerc, parallel=parallel)
	edges$p.value = NA
	for(e in rownames(edges)){
		## compare inclusion percentage with background data
		edges[e, "p.value"] = pnorm(edges[e, "incl.perc"], mean=mean(background.inclusionPerc[, e]), sd=sd(background.inclusionPerc[, e])+0.000000000001, lower.tail=F)
	}
	
	## generate adjacency matrix from selected edges
	adjacency = getAdjacency(edges, variables, feature="p.value")
	
	return(list(edges=edges, adjacency=adjacency))
}







##########################################################################################################################################
## FUNCTIONS
##########################################################################################################################################
## from http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R
## see "Variable Selection with error control: another look at stability selection" (Shah, 2013)
r.TailProbs <- function(eta, B, r) {
  # TailProbs returns a vector with the tail probability for each \tau = ceil{B*2\eta}/B + 1/B,...,1
  # We return 1 for all \tau = 0, 1/B, ... , ceil{B*2\eta}/B
  MAXa <- 100000
  MINa <- 0.00001
  s <- -1/r
  etaB <- eta * B
  k_start <- (ceiling(2 * etaB) + 1)
  output <- rep(1, B)
  if(k_start > B) return(output)
  
  a_vec <- rep(MAXa,B)
  
  Find.a <- function(prev_a) uniroot(Calc.a, lower = MINa, upper = prev_a, tol = .Machine$double.eps^0.75)$root
  Calc.a <- function(a) {
    denom <- sum((a + 0:k)^(-s))
    num <- sum((0:k) * (a + 0:k)^(-s))
    num / denom - etaB
  }
  
  for(k in k_start:B) a_vec[k] <- Find.a(a_vec[k-1])
  
  OptimInt <- function(a) {
    num <- (k + 1 - etaB) * sum((a + 0:(t-1))^(-s))
    denom <- sum((k + 1 - (0:k)) * (a + 0:k)^(-s))
    1 - num / denom
  }
  # NB this function makes use of several gloabl variables
  
  prev_k <- k_start
  for(t in k_start:B) {
    cur_optim <- rep(0, B)
    cur_optim[B] <- OptimInt(a_vec[B])
    if (prev_k <= (B-1)) {
      for (k in prev_k:(B-1))
        cur_optim[k] <- optimize(f=OptimInt,lower = a_vec[k+1], upper = a_vec[k], maximum  = TRUE)$objective
    }
    output[t] <- max(cur_optim)
    prev_k <- which.max(cur_optim)
  }
  return(output)
}

minD <- function(theta, B, r = c(-1/2, -1/4)) {
  pmin(c(rep(1, B), r.TailProbs(theta^2, B, r[1])), r.TailProbs(theta, 2*B, r[2]))
}
