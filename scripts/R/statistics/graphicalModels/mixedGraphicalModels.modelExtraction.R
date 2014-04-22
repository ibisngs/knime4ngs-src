########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="graphicalModels.R", description="This Rscript creates graphical models from datamatrices")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file"        , metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input data file (cols=variables, rows=observations)", metavar="<path>")
parser$add_argument( "--edgeslist"   , type="character", action="store"     , dest="file.el"    ,                help="path to edges list output file"  , metavar="<path>")
parser$add_argument( "--adjacency"   , type="character", action="store"     , dest="file.adj"   ,                help="path to adjacency matrix output file"  , metavar="<path>")

## ARGUMENTS
parser$add_argument("-e", "--ev"      , type="integer"   , action="store"    , dest="ev"      , required=TRUE, help="upper limit for false positive edges in the resulting graph (FWER control)" , metavar="<int>")
parser$add_argument("-p", "--percIncl", type="double"    , action="store"    , dest="percIncl", default=0.8  , help="percentage of stability selection samples in which the edge has to be contained to be included in the final model" , metavar="<double>")

## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)

## CRAN
loadLib("plyr")
require(kimisc)
CWD = normalizePath(if(is.null(thisfile())){getwd()}else{dirname(thisfile())})
source(paste(CWD, "/mixedGraphicalModels.R", sep=""))


##########################################################################################################################################
## READ DATA
##########################################################################################################################################
edge.ranks <- read.csv3(args$file.in)


##########################################################################################################################################
## CREATE MODEL
##########################################################################################################################################
model = getGraph(edge.ranks,
                 stabSel.inclusionPerc=args$percIncl, #
                 E.v=args$ev)
                     


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