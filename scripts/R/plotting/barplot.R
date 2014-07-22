initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(paste(sep="/", script.basename, "plotting.R"))

########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = Create Barplot from Data"", epilogue = "(c) Jonas Zierer")

## add global plotting args
parser = plotting.addArgs(parser)

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


##############################################################################################################
## PLOTTING
##############################################################################################################   
## create plot
p <- ggplot(data, aes_string(x=args$col.x, y=args$col.y, color=args$col.color, fill=args$col.fill)) + 
	geom_bar(stat = "identity")


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
	