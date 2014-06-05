initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(paste(sep="/", script.basename, "plotting.R"))

########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="barplot.R", description="Create Barplot from Data")

plotting.addArgs(parser)


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
	