########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "The Rscript reads a datamatrix containing measured variables (columns) for several samples (rows) and calculates GGMs from the data.", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global", help="path to globals file"        , metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c("--data"         ), type="character", action="store"     , dest="file.in"    , help="path to input data file (cols=variables, rows=observations)", metavar="<path>")
parser <- add_option(parser, c("--output"       ), type="character", action="store"     , dest="file.edges" , help="path to output file (edges list)" , metavar="<path>")

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
	warning("mandatory input file (--data) missing!")
	q(status=-1)
}
if(is.null(args$file.out1)){
	print_help(parser)
	warning("mandatory output file (--output) missing!")
	q(status=-1)
}




########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
## CRAN
loadLib("GeneNet")

##########################################################################################################################################
## READ DATA
##########################################################################################################################################
data <- read.csv3(args$file.in)

##########################################################################################################################################
## CALCULATE GAUSSIAN GRAPHICAL MODEL
##########################################################################################################################################
ggm   <- ggm.estimate.pcor(data)
edges <- ggm.test.edges(ggm, fdr=T, plot=F, direct=F)

## add node names instead of numbers
nodenames = colnames(data)
for(i in 1:nrow(edges)){
	edges[i, "node1"] = nodenames[as.numeric(edges[i, "node1"])]
	edges[i, "node2"] = nodenames[as.numeric(edges[i, "node2"])]
}

## order cols
edges = edges[, c("node1", "node2", "pcor", "pval", "qval", "prob")]
colnames(edges) = c("v1", "v2", "pcor", "p", "q", "prob")
rownames(edges) = paste(edges$v1, edges$v2, sep="~")


##########################################################################################################################################
## WRITE DATA
##########################################################################################################################################
write.csv3(edges, args$file.edges)





# require(huge)
# L = huge.generator(n = 100, d = 30, graph = "hub", g = 6)
# #graph path estimation using mb
# out1 = huge(L$data)
# out1
# plot(out1) #Not aligned
# plot(out1, align = TRUE) #Aligned
# huge.plot(out1$path[[3]])
# select = huge.select(out1)