initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(paste(sep="/", script.basename, "plotting.R"))

########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "create histogram from given data", epilogue = "(c) Jonas Zierer")

## add global plotting args
parser = plotting.addArgs(parser)

## add specific args
parser <- add_option(parser, c("-b" , "--binwidth" ), type="double", action="store"     , dest="binwidth"  , help="width of image"                , metavar="<double>")
parser <- add_option(parser, c(       "--dens"     ),                action="store_true", dest="dens"      , help="plot density instead of counts" )
parser <- add_option(parser, c(       "--densCurve"),                action="store_true", dest="dens.cur"  , help="plot density curve layer" )
parser <- add_option(parser, c(       "--densColor"),                action="store"     , dest="dens.color", help="color of density curve"        , metavar="<color>")

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
p <- plotting.histogram(data, args)


##############################################################################################################
## WRITE OUTPUT
##############################################################################################################  
plotting.print(p, args)
	