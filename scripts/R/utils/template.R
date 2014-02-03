
########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="template.R", description="This is a dummy template R script that does nothing meaningful")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file"        , metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input file"         , metavar="<path>")
parser$add_argument( "--output"      , type="character", action="store"     , dest="file.out1"  , required=TRUE, help="path to first output file"  , metavar="<path>")
parser$add_argument( "--output2"     , type="character", action="store"     , dest="file.out2"  , required=TRUE, help="path to second output file" , metavar="<path>")

## ARGUMENTS
parser$add_argument("-i", "--int"    , type="integer"  , action="store"     , dest="int"        , required=TRUE, help="an integer"                 , metavar="<integer>")
parser$add_argument("-d","--double"  , type="double"   , action="store"     , dest="double"     , required=TRUE, help="a double number"            , metavar="<double>")
parser$add_argument("-s","--string"  , type="character", action="store"     , dest="string"     , required=TRUE, help="a string"                   , metavar="<string>")
parser$add_argument("-l","--list"    , type="character", action="store"     , dest="list"       , required=TRUE, help="a comma-separated list"     , metavar="<comma-separated list>")
parser$add_argument("-b","--bool"    ,                   action="store_true", dest="bool"       ,                help="a boolean flag")
 
## parse
print(commandArgs(trailingOnly=TRUE))
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))


########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)
loadLib("plyr")


########################################################################################################################################
## READ INPUT
########################################################################################################################################
data = read.csv3(args$file.in)

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


