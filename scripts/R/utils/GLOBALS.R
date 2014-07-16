################################################
### VARS
################################################
## FILES
FILE_SEP <- ";"
XLS_DATE_ORIGIN = "1899-12-30"

################################################
### LIBRARIES
################################################
loadLib <- function(x, bioC=FALSE){
	isInstalled =  is.element(x, installed.packages()[,1])
	
	if(!isInstalled){
		cat("Trying to install ", x, "...")
		if(bioC){
			source("http://bioconductor.org/biocLite.R")
			biocLite(x)
		}else{
			success = install.packages(x, repos="http://cran.us.r-project.org")
			cat(success)
		}
	}
	require(x, character.only=T)
	cat("... done\n")
}

################################################
### FUNCTIONS
################################################
## write ; separated csv
write.csv3 <- function(x, file = "", row.names.name="row.names", quote=TRUE, sep=FILE_SEP, row.names = TRUE,col.names = TRUE, ...){
	options("scipen" = 100) # prevent scientific notation
	append=FALSE
	if(row.names & col.names){
		write.table(row.names.name, file=file, append=FALSE, row.names=FALSE, col.names=FALSE, eol=FILE_SEP, sep=sep, ...)
		append=TRUE
	}
	write.table(x, file=file, append=append, row.names=row.names,col.names=col.names, sep=sep)
	options("scipen" = 0) # back to default
}

read.csv3 <- function(filepath, header=T, rownames=T, sep=FILE_SEP, check.names=FALSE, fill=T, ...){
	data <- read.table(file=filepath, header=header, sep=sep, check.names=check.names, fill=fill, ...)
	if(rownames){
		row.names(data) <- data[,1]
		data <- data[,2:ncol(data), drop=F]
	}
	return(data)
}

## get parent directory of file
parent <- function(x){
	return(sub("/+[^/]*$", "", x, perl=T))
}

## clear screen
clear = function(){
	system("clear")
}


## merge two data.frames by rownames
merge.byrownames = function(x,y, ...){
  d = merge(x, y, by=0, ...)
  rownames(d) = d$Row.names
  d = d[, !grepl("Row.names", colnames(d))]
  
  return(d)
}
merge.all = function(l, by...){
  if(length(l)==0){
	return(null)
  }else if(length(l)==1){
	return(l[[1]])
  }else{
	result = l[[1]]
	for(i in 2:length(l)){
		result = merge(result, l[[i]], by=by, all=T, suffixes=F)
	}
	return(result)
  }
}

## split string in named list "name1=value1,name2=value2"
strsplit.named <- function(l, split.1="\\s*,\\s*", split.2="\\s*=\\s*", perl=T){
	if(is.na(l) | is.null(l) | l==""){
		return(list())
	}
	split=as.list(unlist(strsplit(l, split=split.1, perl=perl)))

	strsplit.inner=function(x){
		y=unlist(strsplit(x, split=split.2, perl=perl))
		if(length(y)==2){
			x=list(y[2])
			names(x)=y[1]
			return(x)
		} else if(length(y)==1){
			x=list("")
			names(x) = y[1]
			return(x)
		}else{
			warning("wrong format of string!")
		} 
	}
	l = as.list(unlist(lapply(split, FUN=strsplit.inner)))
	return(l)
}








################################################
### CODING OF VARIABLES
################################################
YES = "yes"
NO  = "no" 

#### recode
## x : variable
## levels: list of new levels of the variable; each level contains a list of old codings which should be replaced
recode = function(x, levels=list(), others=NA){
	for(l in names(levels)){
		x[ x %in% levels[[l]] ] = l
	}
	x[ !x %in% names(levels) ] = others
	return(x)
}
recode.binary = function(x, yes=c("y", "Y", "yes"), no=c("n", "N", "No", "NO", "no"), others=NA){
	recode(x, list('yes'=yes, 'no'=no), others=others)
}


################################################
### FILEPATHS
################################################
thisfile <- function() {
  if (!is.null(res <- thisfile_source())) res
  else if (!is.null(res <- thisfile_rscript())) res
  else NULL
}

# Helper functions
thisfile_source <- function() {
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source))
      return (normalizePath(sys.frame(i)$ofile))
  }

  NULL
}

thisfile_rscript <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  cmdArgsTrailing <- commandArgs(trailingOnly = TRUE)
  cmdArgs <- cmdArgs[seq.int(from=1, length.out=length(cmdArgs) - length(cmdArgsTrailing))]
  res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmdArgs)

  # If multiple --file arguments are given, R uses the last one
  res <- tail(res[res != ""], 1)
  if (length(res) > 0)
    return (res)

  NULL
}

getCWD <- function(){
	return(normalizePath(if(is.null(thisfile())){getwd()}else{dirname(thisfile())}))
}

