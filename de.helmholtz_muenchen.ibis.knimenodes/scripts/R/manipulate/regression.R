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
parser <- ArgumentParser(prog="normalize.R", description="")
parser <- OptionParser(usage = "usage: %prog [options]", description = "Perform regression", epilogue = "(c) Jonas Zierer")

## GLOBALS 
parser <- add_option(parser, c("-g", "--globals"), type="character", action="store"     , dest="file.global"                        , help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser <- add_option(parser, c( "--input"            ), type="character", action="store"     , dest="file.in"                       , help="path to input file", metavar="<path>")
parser <- add_option(parser, c( "--coef"             ), type="character", action="store"     , dest="file.coef"                     , help="path to coefficients output file", metavar="<path>")
parser <- add_option(parser, c( "--residuals"        ), type="character", action="store"     , dest="file.output"                   , help="path to output file", metavar="<path>")


## ARGUMENTS
parser <- add_option(parser, c("-t","--target"       ), type="character", action="store"     , dest="target"                        , help="target variable"  , metavar="<colname>")
parser <- add_option(parser, c("-f","--target.family"), type="character", action="store"     , dest="target.fam", default="gaussian", help="family of target variable"  , metavar="<colname>")
parser <- add_option(parser, c("-c","--confounders"  ), type="character", action="store"     , dest="confounders"                   , help="Independent Variables"  , metavar="<var1,var2,...>")
parser <- add_option(parser, c(       "--q"          )                  , action="store_true", dest="fail.on.na"                    , help="fail if datatable contains NAs" )


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
	warning("mandatory input file (--input) missing!")
	q(status=-1)
}
if(is.null(args$file.coef)){
	print_help(parser)
	warning("mandatory output file (--coef) missing!")
	q(status=-1)
}
if(is.null(args$file.output)){
	print_help(parser)
	warning("mandatory output file (--residuals) missing!")
	q(status=-1)
}

if(is.null(args$target)){
	print_help(parser)
	warning("mandatory parameter 'target' (--target) missing!")
	q(status=-1)
}
if(is.null(args$confounders)){
	print_help(parser)
	warning("mandatory parameter 'confounder' (--confounders) missing!")
	q(status=-1)
}
########################################################################################################################################
## LOAD LIBRARIES
########################################################################################################################################
source(args$file.glob)


##############################################################################################################
## READ DATA
##############################################################################################################
data <- read.csv3(args$file.in)

args$confounders = gsub(",", " + ", args$confounders)
na.action = "na.omit"
if(args$fail.on.na){
	na.action="na.fail"
}
##############################################################################################################
## CALCULATE MODEL
##############################################################################################################
glm = glm(as.formula(paste(args$target, args$confounders, sep="~")), data=data, family=args$target.fam, na.action=na.action)

glm.0 <- glm(as.formula(paste(args$target, 1, sep="~"))            , data=data, family=args$target.fam, na.action=na.action)
anov = anova(glm, glm.0, test="F")

## coefficients
coefficients = as.data.frame(summary(glm)$coefficients)
coefficients = rbind(coefficients, anova=c(NA,NA,NA,anov[2, "Pr(>F)"]))

## residuals
output = data.frame(residuals = glm$residuals, predictors=glm$linear.predictors, fitted=glm$fitted.values, effects=glm$effects, weights=glm$weights, prior.weights=glm$prior.weights)
rownames(output) = rownames(data)
#colnames(residuals) = args$target

##############################################################################################################
## write output
##############################################################################################################
write.csv3(output      , args$file.output)
write.csv3(coefficients, args$file.coef)





