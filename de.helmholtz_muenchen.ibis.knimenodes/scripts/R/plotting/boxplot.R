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
parser <- OptionParser(usage = "usage: %prog [options]", description = "Create Boxplot from given data", epilogue = "(c) Jonas Zierer")

## add global plotting args
parser = plotting.addArgs(parser)

## add specific args
parser <- add_option(parser, c("-p" , "--points"), type="character", action="store"   , dest="points" , default="outliers", help="define if seperate points be plotted over the box"                         , metavar="<outliers, no, all, all jittered>")
parser <- add_option(parser, c("--columns"      ), type="character", action="store"   , dest="columns"                    , help="if one box per column shall be plotted either specify a list of columns or a (perl) regex to define boxes" , metavar="<colname1,colname2, ...>")

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
	if(length(args$columns)==1){
		args$columns = colnames(data)[grepl(args$columns, colnames(data), perl=T)]
	}
	data.m = melt(data[, args$columns])
	data = data.frame(data[, !colnames(data)%in%args$columns], data.m, row.names=c(1:nrow(data.m)))
	args$col.x="variable"
	args$col.y="value"
}

args = plotting.checkArgs(args)


##############################################################################################################
## PLOTTING
##############################################################################################################   
p = plotting.boxplot(data, args)

##############################################################################################################
## WRITE OUTPUT
##############################################################################################################  
plotting.print(p, args)
