########################################################################################################################################
## PARSE ARGS
########################################################################################################################################
require(argparse)
parser <- ArgumentParser(prog="regressionBatch.R", description="Fit Regression Models for multiple Variables")

## GLOBALS 
parser$add_argument("-g", "--globals", type="character", action="store"     , dest="file.global", required=TRUE, help="path to globals file", metavar="<path>")

## IN- AND OUTPUT-FILES
parser$add_argument( "--input"       , type="character", action="store"     , dest="file.in"    , required=TRUE, help="path to input file", metavar="<path>")
parser$add_argument( "--residuals"        , type="character", action="store"     , dest="file.residuals"  , required=TRUE, help="path to residuals output file", metavar="<path>")
parser$add_argument( "--output"   , type="character", action="store"     , dest="file.output"   , required=TRUE, help="path to output file", metavar="<path>")


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

args$target      = unlist(strsplit(args$target,","))
args$confounders = gsub(",", " + ", args$confounders)

na.action = "na.omit"
if(args$fail.on.na){
	na.action="na.fail"
}
##############################################################################################################
## CALCULATE MODEL
##############################################################################################################
residuals = data.frame(row.names=rownames(data))
results   = data.frame()
for(t in args$target){
	glm = glm(as.formula(paste(t, args$confounders, sep="~")), data=data, family=args$target.fam, na.action=na.action)
	
	glm.0 <- glm(as.formula(paste(t, 1, sep="~"))            , data=data, family=args$target.fam, na.action=na.action)
	anov = anova(glm, glm.0, test="F")

	## coefficients
	results = rbind(results, data.frame(target=t, p=anov[2, "Pr(>F)"], n=sum(complete.cases(data[, c(t,args$confounders)]))))

	
	## residuals
	residuals[names(glm$residuals), t] = glm$residuals
}
##############################################################################################################
## write output
##############################################################################################################
data[, args$target] = residuals[, args$target]
rownames(results) = results$target

write.csv3(results  , args$file.output)
write.csv3(data     , args$file.residuals)





