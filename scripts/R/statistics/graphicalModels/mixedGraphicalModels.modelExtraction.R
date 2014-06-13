########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="mixedGraphicalModels.modelExtraction.R", description="Extract mixed graphical models from sampled edgerankings")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file"        , metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input data file with edge ranks(cols=variables, rows=observations)", metavar="<path>")
parser$add_argument( "--background"  , type="character", action="store"     , dest="file.backgr",                help="path to input data file with edge ranks from randomly shuffled data for background distribution estimation(cols=variables, rows=observations). If this parameter is used, empirical p-values are calculated and the --ev parameter is ignored.", metavar="<path>")
parser$add_argument( "--edgeslist"   , type="character", action="store"     , dest="file.el"    ,                help="path to edges list output file"  , metavar="<path>")
parser$add_argument( "--adjacency"   , type="character", action="store"     , dest="file.adj"   ,                help="path to adjacency matrix output file"  , metavar="<path>")

## ARGUMENTS
parser$add_argument("-p", "--percIncl", type="double"    , action="store"    , dest="percIncl"  , default=0.8  , help="percentage of stability selection samples in which the edge has to be contained to be included in the final model" , metavar="<double>")

parser$add_argument("-e", "--ev"      , type="integer"   , action="store"    , dest="ev"        ,                help="(FWER control) upper limit for false positive edges in the resulting graph" , metavar="<int>")
parser$add_argument("-n", "--nedges"  , type="integer"   , action="store"    , dest="nEdges"    ,                help="(empirical p-values) number of edges to be selected in each sampled model" , metavar="<integer>")
parser$add_argument("-c", "--cores"   , type="integer"   , action="store"    , dest="cores"     ,                help="number of cores used for calculations" , metavar="<integer>")



## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)

## CRAN
loadLib("plyr")
CWD = getCWD()
source(paste(CWD, "/mixedGraphicalModels.R", sep=""))


##########################################################################################################################################
## READ DATA
##########################################################################################################################################
edge.ranks <- read.csv3(args$file.in)

## parallelization
parallel = FALSE
if(is.null(args$cores) | args$cores<1){
	loadLib("parallel")
	args$cores = detectCores()-1
}
if(args$cores > 1){
	loadLib("doMC")
	registerDoMC(cores=args$cores)
	warning("USING ", args$cores, " threads for computations!")
	parallel = TRUE
}

##########################################################################################################################################
## CREATE MODEL
##########################################################################################################################################
## EMPIRICAL P_VALUES
if(!is.null(args$file.backgr)){
	if(is.null(args$nEdges)){
		stop(paste("Number of Edges to be selected in each subsample has to be given for empirical p-value calculation (--nedges parameter)\n", parser$print_help(), sep=""))
	}
	edge.ranks.background <- read.csv3(args$file.backgr)
	model = getGraph.empiricalPvalues(edge.ranks,
	                                  edge.ranks.background,
	                                  stabSel.sampleNumedges = args$nEdges,
	                                  stabSel.inclusionPerc=args$percIncl,
	                                  parallel = parallel)


## FWER CONTROL
}else{
	if(is.null(args$ev)){
		stop(paste("Upper limit for false positives must be given for FWER control(--ev parameter)\n", parser$print_help(), sep=""))
	}
	model = getGraph.FWER(edge.ranks,
	                      stabSel.inclusionPerc=args$percIncl, #
	                      E.v=args$ev,
	                      parallel = parallel)
}


##########################################################################################################################################
## OUTPUT RESULTS MODEL
##########################################################################################################################################
if(! is.null(args$file.el)){
	write.csv3(model$edges, args$file.el)
}

if(! is.null(args$file.adj)){
	write.csv3(model$adjacency, args$file.adj)
}
#write.graph(model$g, "~/g.graphml",format="graphml")