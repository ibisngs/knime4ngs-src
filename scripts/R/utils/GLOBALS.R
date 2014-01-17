################################################
### VARS
################################################
## FILES
FILE_SEP <- ";"


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
		data <- data[,2:ncol(data)]
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