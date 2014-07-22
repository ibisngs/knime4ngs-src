initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(paste(sep="/", script.basename, "plotting.R"))

########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "Create Boxplot from given data", epilogue = "(c) Jonas Zierer")

## add global plotting args
parser = plotting.addArgs(parser)

## add specific args
parser <- add_option(parser, c("-p" , "--points"), type="character", action="store"   , dest="points" , default="outliers", help="define if seperate points be plotted over the box"                         , metavar="<outliers, no, all, all jittered>")
parser <- add_option(parser, c("--columns"      ), type="character", action="store"   , dest="columns"                    , help="if one box per column shall be plotted (instead of using col.x and col.y)" , metavar="<colname1,colname2, ...>")

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
if(!is.null(args$columns) && args$columns!=""){
	loadLib("reshape2")
	args$columns = unlist(strsplit(args$columns, ","))
	data.m = melt(data[, args$columns])
	data = data.frame(data[, !colnames(data)%in%args$columns], data.m, row.names=c(1:nrow(data.m)))
	args$col.x="variable"
	args$col.y="value"
}

args = plotting.checkArgs(args)


##############################################################################################################
## PLOTTING
##############################################################################################################   
## create plot
p <- ggplot(data, aes_string(x=args$col.x, y=args$col.y, fill=args$col.fill))


## create boxplot
#<outliers, no, all, all jittered>
if(args$points == "outliers"){
	p <- p + geom_boxplot() ## todo add shape and or color
}else{
	p <- p + geom_boxplot(outlier.shape = NA) 
	if(args$points == "all" | args$points == "all jittered"){
		position="identity"
		if(args$points == "all jittered"){
			position = "jittered"
		}
		p <- p + geom_point(aes_string(color=args$col.color, shape=args$col.shape), position=position) ## todo add she and color
	}
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
