initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(paste(sep="/", script.basename, "plotting.R"))

########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="histogram.R", description="Create Histogram from Data")

plotting.addArgs(parser)

parser$add_argument("-b" , "--binwidth"  , type="double"  , action="store"      , dest="binwidth"         , help="width of image"  , metavar="<double>")
parser$add_argument(       "--dens"      , action="store_true",  dest="dens"     ,help="plot density instead of counts" )
parser$add_argument(       "--densCurve" , action="store_true",  dest="dens.cur" , help="plot density curve layer" )


## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

##############################################################################################################
## LOAD PACKAGES
##############################################################################################################
source(args$file.global)
loadLib("ggplot2")


##############################################################################################################
## READ FILES
##############################################################################################################
data <- plotting.readData(args)

##############################################################################################################
## PROCESS PARAMS
##############################################################################################################
args = plotting.checkArgs(args)

## bin width 
if(is.null(args$binwidth) || args$binwidth<=0){
	args$binwidth = NULL
}

##############################################################################################################
## PLOTTING
##############################################################################################################   
## create plot
p <- ggplot(data, aes_string(x=args$col.x, color=args$col.color, fill=args$col.fill))


## histogram
if(args$dens){
	p <- p + geom_histogram(aes(y=..density..))
	if(args$dens.cur){
		p <- p +  geom_density(alpha=.2)
	}
}else{
	p <- p + geom_histogram(binwidth=args$binwidth)
}


## add facets
p = plotting.addFacets(p, args)

## add scales, labels and guides
p = plotting.addScalesAndLabelsAndGuides(p, args)

## change layout
p <- p + geom_default(legend=="vertical")


##############################################################################################################
## WRITE OUTPUT
##############################################################################################################  
plotting.print(p, args)
	