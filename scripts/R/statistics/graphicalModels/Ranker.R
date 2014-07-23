########################################################################################################################################
## RANKER
########################################################################################################################################
## INPUT 
##	dep.var = name (column name) of dependent variable
##	D       = data matrix containing variables in cols and samples in rows
##	...     = further arguments which are passed to the ranking function
## OUTPUT      
##	vector named and ordered as the columns of D 
##	containing the node ranks (low=dep.var depends on this variable, high=no relationship with dep.var)



############################################################################################################
## RANDOM FOREST RANKER
############################################################################################################
## calculate variable importance for all variables for the dep.var
## ADDITIONAL PARAMS
##	type  = how to asses variable importance? (1=mean decrease in accuracy, 2=mean decrease in node impurity)
##	...   = passed to randomForest function
##
##
ranker.randomForest <- function(dep.var, D, type=1, ...){
	dep.var.idx = which(dep.var == names(D))
	## preprocess outcome and cancel if variable has too little variation
	y = D[, dep.var]
	if(is.numeric(y)){
		if(length(unique(y)) < 6){
			return(rep(NA, ncol(D)))
		}
	}else if(is.factor(y)){
		y = droplevels(y)
		if(length(levels(y))==1){
			return(rep(NA, ncol(D)))
		}
	}

	var.imp = importance(randomForest(x=D[,-dep.var.idx], y=y, importance=TRUE, proximity=FALSE, keep.forest=FALSE, ...), type=type)[,1]
	var.imp[dep.var] = NA
	return(var.imp[colnames(D)])
}
############################################################################################################
## RANDOM FOREST RANKER
############################################################################################################
## calculate variable importance for all variables for the dep.var
## ADDITIONAL PARAMS
##	type  = how to asses variable importance? (1=mean decrease in accuracy, 2=mean decrease in node impurity)
##	...   = passed to cforest_unbiased (control) function
##
##
ranker.randomForest.party <- function(dep.var, D, type="AUC", mtry=NA, ...){
	y= D[, dep.var]
	## preprocess outcome and cancel if variable has too little variation
	if(is.factor(y)){
		y = droplevels(y)
		if(length(levels(y))==1){
			return(rep(NA, times=ncol(D)))
		}
	}
	if(is.na(mtry)){
		mtry = if (!is.null(y) && !is.factor(y)) max(floor(ncol(D)/3), 1) else floor(sqrt(ncol(D)))
	}
	
	forest = cforest(as.formula(paste(dep.var, "~ .")), data = D, control = cforest_unbiased(mtry=mtry, ...))
	if(type=="AUC"){
		var.imp = varimpAUC(forest)
	}else{
		var.imp = varimp(forest)
	}

	var.imp[dep.var] = NA
	#ranking[names(var.imp)] = rank(-var.imp)
	return(var.imp[colnames(D)])
}


############################################################################################################
## LASSO RANKER
############################################################################################################
## run lasso regression using glmnet on depend variable
## rank all variables according to the lambda for which the according beta becomes >0
## ADDITIONAL PARAMS
##	...  = passed to glmnet function (particularly 'family' has to be given if not gaussian)
## 
##
## suppress pmax warnings and not reached maxit warnings
suprWarning.glmnet <- function(w){
	if( any( grepl( "exceeds pmax=", w) || grepl("not reached after maxit=", w ) )){
		#warning("SUPRESS!")#
		invokeRestart( "muffleWarning" )
	}
}
ranker.lasso <- function(dep.var, D, nlambda=100, ...){
	dep.var.idx = which(dep.var == names(D))
	
	## preprocess outcome and cancel if variable has too little variation
	if(is.factor(D[, dep.var])){
		D[, dep.var] = droplevels(D[, dep.var])
		if(length(levels(D[, dep.var]))==1){
			return(rep(NA, times=ncol(D)))
		}
	}
	
	withCallingHandlers(
		lasso.cv  <- glmnet(data.matrix(D[, -dep.var.idx]), D[, dep.var], standardize=T, nlambda=nlambda, ...), #pmax=2*round((9)/2)-1#, # family="gaussian", 
 		warning = suprWarning.glmnet
 	)
	## get beta coefficients for all lamda parameters
	betas <- lasso.cv$beta
	# in case of multinomial regression one beta matrix for each level -> take lowest possible lambda from all matrices
	if(is.list(betas)){
		betas = Reduce(function(x,y){abs(x)+abs(y)}, betas)
	}
	#betas <- betas[rowSums(betas)!=0,]                # ignore variables which are never used
	bestLambda = aaply(.data=betas, .margins=1, .fun=function(x){a = which(x!=0); if(length(a)==0)return(-Inf); return(min(a))})
	bestLambda[dep.var] = NA
	
	return(nlambda-bestLambda[colnames(D)])
}


############################################################################################################
## REGRESSION RANKER
############################################################################################################
## run regression using glm on depend variable
## run an anova analysis and use the explained Deviance as ranking criterion
## ADDITIONAL PARAMS
##	...  = passed to glmnet function (particularly 'family' has to be given if not gaussian)
## 
##
ranker.regression <- function(dep.var, D, ...){
	## preprocess outcome and cancel if variable has too little variation
	if(is.factor(D[, dep.var])){
		D[, dep.var] = droplevels(D[, dep.var])
		if(length(levels(D[, dep.var]))==1){
			return(rep(NA, times=ncol(D)))
		}
	}
	
	model <- glm(as.formula(paste(dep.var, "~ .")), data=D, ...)
	ano = anova(model)
	ranking = ano[, "Deviance"]
	names(ranking) = rownames(ano)
	ranking[dep.var] = NA
# 	betas = model$coefficients[names(model$coefficients) != "(Intercept)"]
# 	names(betas) = colnames(D)[colnames(D)!=dep.var]
# 	betas[dep.var] = NA
# 	betas = betas[names(betas) != "(Intercept)"]
# 	return(betas[colnames(D)])
	return(ranking[colnames(D)])
}





############################################################################################################
## RANK EDGES
############################################################################################################
## INPUT
##	D             = data matrix containing variables in cols and samples in rows
##	var.classes   = class of variables contained in D (e.g. numeric, factor,...)
##	ranker        = list containing one entry foreach class in var.classes describing the ranker used for this class
##	ranker.params = list containing one entry foreach class in var.classes. Each entry is a list of parameters for this ranker
##	rank.type     = local or global ranking?
##	edges.indices = order of edges (how to cast adjacency matrix to vector). Should be the indices of 'upper.triangle' by default
##
## OUTPUT
##	vector of edges (ordered as inedges.indices). One rank for each edge is returned, the smaller the more likely it is a true edge
##
##
rankEdges <- function(D, var.classes, 
                      ranker        = list(numeric="randomForest", double="randomForest" , factor="randomForest"),
                      ranker.params = list(numeric=list() , double=list()  , factor=list()),
                      rank.type = "local",
                      edges.indices, 
                      parallel.ranking=FALSE){
	p = ncol(D)
	## ranking for all types of variables
	ranks = laply(.data=c(1:p), .parallel=parallel.ranking, .inform=F,
		.fun=function(i){
			do.call(paste("ranker", ranker[[var.classes[i]]], sep="."), c(list(colnames(D)[i], D) , ranker.params[[var.classes[i]]]))
		})

	#rownames(ranks) = colnames(ranks)
	ranks[is.na(ranks)] = -10e10
	if(rank.type == "local"){
		# ranks row-wise
		ranks = aaply(ranks, .margins=1, .fun=function(x){(rank(x)-1)/(p-1)})
	}else if(rank.type == "global"){
		ranks = matrix((rank(ranks)-1)/(p*p-1), nrow=p, ncol=p)#, dimnames=list(colnames(D), colnames(D))
	}else{
		stop(paste("Unknown rank type", rank.type))
	}
	
	## combine local rankings
	ranks <- get.max.upper.triangle(ranks)
	ranks <- c(ranks[edges.indices])

	return(ranks)
}


############################################################################################################
## GET MAX UPPER TRIANGLE
############################################################################################################
## INPUT
##	A             = symmetric matrix
##
## OUTPUT
##	A matrix with only the upper trianlge filled. Each entry is the maximum of the to symmetric entries
##	Maximum of each pair A[i,j], A[[j,i] in upper triangle of matrix
##
get.max.upper.triangle <- function(A){
    if (dim(A)[1]!=dim(A)[2])
        stop("Not a quadratic matrix!")

    p <- dim(A)[1]
    ## True/False check: Is entry in upper triangle larger than the one in lower triangle
    DIFF <- A[upper.tri(A)] - t(A)[upper.tri(A)] >= 0
    A.up.tri.max <- matrix(0,p,p)
    ## save all values where _upper_ tri is larger
    A.up.tri.max[upper.tri(A)][DIFF] <- A[upper.tri(A)][DIFF]
    ##  save all values where _lower__ tri is larger
    A.up.tri.max[upper.tri(A)][!DIFF] <- t(A)[upper.tri(A)][!DIFF]
    return(A.up.tri.max)
}

############################################################################################################
## GET MIN UPPER TRIANGLE
############################################################################################################
## INPUT
##	A             = symmetric matrix
##
## OUTPUT
##	A matrix with only the upper trianlge filled. Each entry is the maximum of the to symmetric entries
##	Minimum of each pair A[i,j], A[[j,i] in upper triangle of matrix
##
get.min.upper.triangle <- function(A){
    if (dim(A)[1]!=dim(A)[2])
        stop("Not a quadratic matrix!")

    p <- dim(A)[1]
    ## True/False check: Is entry in upper triangle larger than the one in lower triangle
    DIFF <- A[upper.tri(A)] - t(A)[upper.tri(A)] < 0
    A.up.tri.min <- matrix(0,p,p)
    ## save all values where _upper_ tri is larger
    A.up.tri.min[upper.tri(A)][DIFF] <- A[upper.tri(A)][DIFF]
    ##  save all values where _lower__ tri is larger
    A.up.tri.min[upper.tri(A)][!DIFF] <- t(A)[upper.tri(A)][!DIFF]
    return(A.up.tri.min)
}

############################################################################################################
## GET MAX UPPER TRIANGLE
############################################################################################################
## INPUT
##	A  = vector of variable names (e.g. colnames of adjacency matrix)
##
## OUTPUT
##	vector of (pasted) pairs of variable names ordered as upper.triangle (e.g. of the adjacency matrix)
##
make.edgesNames <- function(x){
	do.call(c, llply(.data=2:length(x), .fun=function(j){paste(x[1:(j-1)], x[j], sep="~")}))
}