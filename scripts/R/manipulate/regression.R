########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="normalize.R", description="Normalize data using differenz methods")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input file", metavar="<path>")
parser$add_argument( "--coef"        , type="character", action="store"     , dest="file.coef"  , required=TRUE, help="path to coefficients output file", metavar="<path>")
parser$add_argument( "--residuals"   , type="character", action="store"     , dest="file.output"   , required=TRUE, help="path to output file", metavar="<path>")


## ARGUMENTS
parser$add_argument("-t","--target"       , type="character", action="store"     , dest="target"    ,required=TRUE, help="target variable"  , metavar="<colname>")
parser$add_argument("-f","--target.family", type="character", action="store"     , dest="target.fam",default="gaussian", help="family of target variable"  , metavar="<colname>")
parser$add_argument("-c","--confounders"  , type="character", action="store"     , dest="confounders",required=TRUE, help="Independent Variables"  , metavar="<var1,var2,...>")
parser$add_argument(       "--q"   , action="store_true", dest="fail.on.na"                                  , help="fail if datatable contains NAs" )


## parse
args <- parser$parse_args(commandArgs(trailingOnly=TRUE))


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





