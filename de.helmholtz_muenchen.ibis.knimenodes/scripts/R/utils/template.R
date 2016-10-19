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

########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(optparse)
parser <- OptionParser(usage = "usage: %prog [options]", description = "this is a dummy template script", epilogue = "(c) Your Name")

## GLOBALS 
# read globals.R which contains function to read/write files in a uniform format, ...
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global", help="path to globals file"        , metavar="<path>")

## ARGUMENTS
# add arguments to parser, various types are possibl. check the optparse 
parser <- add_option(parser, c( "--input"       ), type="character", action="store"     , dest="file.in"    , help="path to input file"         , metavar="<path>")
parser <- add_option(parser, c( "--output"      ), type="character", action="store"     , dest="file.out1"  , help="path to first output file"  , metavar="<path>")
parser <- add_option(parser, c( "--output2"     ), type="character", action="store"     , dest="file.out2"  , help="path to second output file" , metavar="<path>")
parser <- add_option(parser, c("-i", "--int"    ), type="integer"  , action="store"     , dest="int"        , help="an integer"                 , metavar="<integer>")
parser <- add_option(parser, c("-d","--double"  ), type="double"   , action="store"     , dest="double"     , help="a double number"            , metavar="<double>")
parser <- add_option(parser, c("-s","--string"  ), type="character", action="store"     , dest="string"     , help="a string"                   , metavar="<string>")
parser <- add_option(parser, c("-l","--list"    ), type="character", action="store"     , dest="list"       , help="a comma-separated list"     , metavar="<comma-separated list>")
parser <- add_option(parser, c("-b","--bool"    ),                   action="store_true", dest="bool"       , help="a boolean flag")
 
 
## parse
# parse the commandline arguments
args = parse_args(parser, args = commandArgs(trailingOnly = TRUE), print_help_and_exit = TRUE, positional_arguments = FALSE)

## mandatory args
# check if some mandatory args are missing
if(is.null(args$file.global)){
	print_help(parser)
	warning("mandatory globals file (--globals) missing!")
	q(status=-1)
}
if(is.null(args$file.in)){
	print_help(parser)
	warning("mandatory input file (--input) missing!")
	q(status=-1)
}

########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
loadLib("plyr") # use loadLib which installs the package automatically if it is missing


########################################################################################################################################
## READ INPUT
########################################################################################################################################
data = read.csv3(args$file.in) # use read.csv3 and write.csv3 (defined in GLOBALS.R) to write in a format that KNIME can read

#Sys.sleep(2000000000)

########################################################################################################################################
## DO INTELLIGENT STUFF
########################################################################################################################################
columns <- unlist(strsplit(args$list, ","))
data <- data[, columns]


data.args = NULL
for(i in 1:length(args)){
	data.args = rbind(data.args, data.frame(names(args)[i], args[[i]]))
}
colnames(data.args) = c("argument", "value")
rownames(data.args) = c(1:nrow(data.args))


########################################################################################################################################
## WRITE OUTPUT
########################################################################################################################################
write.csv3(data, args$file.out1)
write.csv3(data.args, args$file.out2)


