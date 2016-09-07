#  Copyright (C) 2016 the Knime4NGS contributors.
#  Website: http://ibisngs.github.io/knime4ngs
#  
#  This file is part of the KNIME4NGS KNIME extension.
#  
#  The KNIME4NGS extension is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.basename <- dirname(script.name)
source(paste(sep="/", script.basename, "plotting.R"))

########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "create scatterplot from given data", epilogue = "(c) Jonas Zierer")

## add global plotting args
parser = plotting.addArgs(parser)

## add specific args
parser <- add_option(parser, c("-m" , "--matrix"), type="character", action="store", dest="matrix"             , help="give a comma separated list of columns for scatterplot matrix (instead of x and y)", metavar="<int>")
parser <- add_option(parser, c("--alpha"        ), type="double"   , action="store", dest="alpha" , default=0.6, help="alpha value of points"                                                             , metavar="<double>" )
parser <- add_option(parser, c("--size"         ), type="double"   , action="store", dest="size"  , default=1  , help="pointsize"                                                                         , metavar="<col1,col2,...>")

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

if(!is.null(args$matrix)){
	args$matrix = unlist(strsplit(args$matrix, ","))
	data.pairs  = plotting.makePairs(data[, args$matrix])
	data        = data.frame(data[, !colnames(data) %in% args$matrix],
						data.pairs$all)
	
	args$col.x = "x"
	args$col.y = "y"
	args$col.facet.x = "xvar"
	args$col.facet.y = "yvar"
	
}

if(is.null(args$col.x) || is.null(args$col.y)){
	warning("x and y col missing!")
}
##############################################################################################################
## PROCESS PARAMS
##############################################################################################################
args = plotting.checkArgs(args)

##############################################################################################################
## PLOTTING
##############################################################################################################   
p = plotting.scatterplot(data, args)


##############################################################################################################
## WRITE OUTPUT
##############################################################################################################  
plotting.print(p, args)
	
