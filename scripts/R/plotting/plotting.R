##############################################################################################################
## PLOTTING ARGS
##############################################################################################################
plotting.addArgs = function(parser){

	## GLOBALS 
	parser <- add_option(parser, c("-g", "--globals"      ), type="character", action="store"      , dest="file.global"                 , help="path to globals file"                                              , metavar="<path>")
	## IN- AND OUTPUT-FILES
	parser <- add_option(parser, c( "--data"              ), type="character", action="store"      , dest="file.in"                     , help="path to input data file (cols=variables, rows=observations)"       , metavar="<path>")

	## ARGUMENTS
	parser <- add_option(parser, c("-x" , "--col.x"       ), type="character", action="store"      , dest="col.x"                       , help="name of column to be plotted on x-axis (defining different boxes)" , metavar="<column name>")
	parser <- add_option(parser, c("-y" , "--col.y"       ), type="character", action="store"      , dest="col.y"                       , help="name of column to be plotted on y-axis (defining different boxes)" , metavar="<column name>")
	# color
	parser <- add_option(parser, c("-c" , "--col.color"   ), type="character", action="store"      , dest="col.color"                   , help="name of column which defines color"                                , metavar="<column name>")
	parser <- add_option(parser, c("-cl", "--lab.color"   ), type="character", action="store"      , dest="lab.color"                   , help="label for color legend"                                            , metavar="<column name>")
	parser <- add_option(parser, c("-cs", "--scale.color" ), type="character", action="store"      , dest="scale.color"                 , help="scale palette of color"                                            , metavar="<string>")
	# fill color
	parser <- add_option(parser, c("-f" , "--col.fill"    ), type="character", action="store"      , dest="col.fill"                    , help="name of column which defines fill"                                 , metavar="<column name>")
	parser <- add_option(parser, c("-fl", "--lab.fill"    ), type="character", action="store"      , dest="lab.fill"                    , help="label for fill color legend"                                       , metavar="<column name>")
	parser <- add_option(parser, c("-fs", "--scale.fill"  ), type="character", action="store"      , dest="scale.fill"                  , help="scale palette of fill color"                                       , metavar="<string>")
	# shape
	parser <- add_option(parser, c("-s" , "--col.shape"   ), type="character", action="store"      , dest="col.shape"                   , help="name of column which defines shape"                                , metavar="<column name>")
	parser <- add_option(parser, c("-sl", "--lab.shape"   ), type="character", action="store"      , dest="lab.shape"                   , help="label for shape legend"                                            , metavar="<string>")
	# facet
	parser <- add_option(parser, c("-fx", "--facet.x"     ), type="character", action="store"      , dest="col.facet.x"                 , help="name of column which defines x axis of facet grid"                 , metavar="<column name>")
	parser <- add_option(parser, c("-fy", "--facet.y"     ), type="character", action="store"      , dest="col.facet.y"                 , help="name of column which defines x axis of facet grid"                 , metavar="<column name>")
	# labels
	parser <- add_option(parser, c("-t" , "--title"       ), type="character", action="store"      , dest="lab.plot"                    , help="Title of plot"                                                     , metavar="<String>")
	parser <- add_option(parser, c("-xl", "--lab.x"       ), type="character", action="store"      , dest="lab.x"                       , help="Label of x axis"                                                   , metavar="<String>")
	parser <- add_option(parser, c("-yl", "--lab.y"       ), type="character", action="store"      , dest="lab.y"                       , help="Label of y axis"                                                   , metavar="<String>")
	# output
	parser <- add_option(parser, c("-w" , "--width"       ), type="integer"  , action="store"      , dest="width"         , default=800 , help="width of image"                                                    , metavar="<int>")
	parser <- add_option(parser, c(       "--height"      ), type="integer"  , action="store"      , dest="height"        , default=600 , help="height of image"                                                   , metavar="<int>")
	parser <- add_option(parser, c("-i" , "--image"       ), type="character", action="store"      , dest="image"         ,             , help="name of png file"                                                  , metavar="<path>")
	
	return(parser)
}



##############################################################################################################
## PLOT IMAGE
##############################################################################################################
plotting.print = function(p, args){
	loadLib("gridExtra")
	png(args$image, width=args$width, height=args$height)
		grid.draw(multiLegendAlign(p))
	dev.off()
}



##############################################################################################################
## MAKE PAIRS FOR SCATTERPLOT MATRIX
##############################################################################################################
# from http://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
plotting.makePairs <- function(data){
	grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
	grid <- subset(grid, x != y)
	all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
		xcol <- grid[i, "x"]
		ycol <- grid[i, "y"]
		data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], x = data[, xcol], y = data[, ycol], data)
	}))
	all$xvar <- factor(all$xvar, levels = names(data))
	all$yvar <- factor(all$yvar, levels = names(data))
	densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
		data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
	}))
	list(all=all, densities=densities)
}


##############################################################################################################
## DEFAULT THEME
##############################################################################################################
geom_default = function(colour=NA, shape=NA, legend=NA){
	loadLib("gridExtra")
	props = theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

	## LEGEND
	if( !is.na(legend)){
		if(legend=="vertical"){
			props <- props + theme(legend.position=c(0,1), legend.justification=c(0,1), legend.direction="vertical", legend.box="vertical")
		}
	}
	return(props)
}


##############################################################################################################
## ALIGN MULTIPLE LEGENDS
##############################################################################################################
multiLegendAlign <- function(p, align="left"){
	loadLib("gridExtra")
	data <- ggplot_build(p)
	gtable <- ggplot_gtable(data)

	# Determining index of legends table
	lbox <- which(sapply(gtable$grobs, paste) == "gtable[guide-box]")
	if(length(lbox)==1){
		# Each legend has several parts, wdth contains total widths for each legend
		# Determining narrower legend
		if(length(gtable$grobs[[lbox]])==1){
			wdth <- with(gtable$grobs[[lbox]], c(sum(as.vector(grobs[[1]]$widths)),sum(as.vector(grobs[[1]]$widths))))
			id <- 1
		}else{
			wdth <- with(gtable$grobs[[lbox]], c(sum(as.vector(grobs[[1]]$widths)),sum(as.vector(grobs[[2]]$widths))))
			id <- which.min(wdth)
		}
		# Adding a new empty column of abs(diff(wdth)) mm width on the right of 
		# the smaller legend box
		if(align=="right"){
			pos = 0
		}else if(align=="left"){
			pos = -1
		}
		gtable$grobs[[lbox]]$grobs[[id]] <- gtable::gtable_add_cols(gtable$grobs[[lbox]]$grobs[[id]], unit(abs(diff(wdth)), "mm"), pos=pos)
	}

	# Plotting
	return(gtable)
}


##############################################################################################################
## CHECK ARGS
##############################################################################################################
replaceInvalidChars = function(x){
	 gsub("[^0-9A-Za-z\\.]", "...", x) 
}

plotting.checkArgs = function(args){
	## x
	args$col.x     = replaceInvalidChars(args$col.x)
	
	## y
	args$col.y     = replaceInvalidChars(args$col.y)
	
	## color
	if( is.null(args$col.color) || is.na(args$col.color)|| args$col.color=="NULL" || args$col.color=="" || args$col.color=="NA" ){
		args$col.color <- NULL
		args$lab.color <- NULL
	}else{
		args$col.color = replaceInvalidChars(args$col.color)
		
		if(is.null(args$lab.color) || args$lab.color == "NULL"){
			args$lab.color = NULL
		}else if(args$lab.color == "NA" || args$lab.color == ""){
			args$lab.color = NA
		}
	}

	## fill
	if( is.null(args$col.fill) || is.na(args$col.fill)|| args$col.fill=="NULL" || args$col.fill=="" || args$col.fill=="NA" ){
		args$col.fill <- NULL
		args$lab.fill <- NULL
	}else{
		args$col.fill = replaceInvalidChars(args$col.fill)
		
		if(is.null(args$lab.fill) || args$lab.fill == "NULL"){
			args$lab.fill = NULL
		}else if(args$lab.fill == "NA" || args$lab.fill == ""){
			args$lab.fill = NA
		}
	}
	
	## shape
	if( is.null(args$col.shape) || is.na(args$col.shape)|| args$col.shape=="NULL" || args$col.shape=="" || args$col.shape=="NA" ){
		args$col.shape <- NULL
		args$lab.shape <- NULL
	}else{
		args$col.shape = replaceInvalidChars(args$col.shape)
		
		if(is.null(args$lab.shape) || args$lab.shape == "NULL"){
			args$lab.shape = NULL
		}else if(args$lab.shape == "NA" || args$lab.shape == ""){
			args$lab.shape = NA
		}
	}

	## X Label
	if( is.null(args$lab.x) || args$lab.x == "" || args$lab.x == "NA" || args$lab.x == "NULL" ){
		args$lab.x = NULL
	}

	## Y Label
	if( is.null(args$lab.y) || args$lab.y == "" || args$lab.y == "NA" || args$lab.y == "NULL" ){
		args$lab.y = NULL
	}

	## Main Label
	if( is.null(args$lab.plot) || args$lab.plot == "" || args$lab.plot == "NA" || args$lab.plot == "NULL" ){
		args$lab.plot = NULL
	}
	
	## Facets
	if( is.null(args$facet.x) || args$facet.x == "" || args$facet.x == "NA" || args$facet.x == "NULL" ){
		args$facet.x = NULL
	}
	if( is.null(args$facet.y) || args$facet.y == "" || args$facet.y == "NA" || args$facet.y == "NULL" ){
		args$facet.y = NULL
	}
	
	## Scales
	if(!is.null(args$scale.color) && args$scale.color!="default"){
		## todo
	}
	
	
	return(args)
}

plotting.readData <- function(args){
	data <- read.csv3(args$file.in)
	data$ROWID = rownames(data)
	colnames(data) = replaceInvalidChars(colnames(data))
	return(data) 
}

##############################################################################################################
## ADD LABEL AND GUIDES
##############################################################################################################
plotting.addScalesAndLabelsAndGuides = function(p, args){
	## scales
	if(!is.null(args$scale.color) && args$scale.color!="default"){
		p <- p + scale_color_brewer(palette=args$scale.color)
	}
	if(!is.null(args$scale.fill) && args$scale.fill!="default"){
		p <- p + scale_fill_brewer(palette=args$scale.fill)
	}
	
	## labels
	p <- p + labs(x=args$lab.x, y=args$lab.y, title=args$lab.plot, color=args$lab.color, shape=args$lab.shape, fill=args$lab.fill)

	## guides
	if(is.null(args$lab.color)){
		p <- p + guides(color="none")
	}
	if(is.null(args$lab.fill)){
		p <- p + guides(fill="none")
	}
	if(is.null(args$lab.shape)){
		p <- p + guides(shape="none")
	}
	return(p)
}



##############################################################################################################
## ADD FACETS
##############################################################################################################
plotting.addFacets = function(p, args){
	if( !is.null(args$col.facet.x) || !is.null(args$col.facet.y)){
		if(is.null(args$col.facet.x)){
			args$col.facet.x = "."
		}
		if(is.null(args$col.facet.y)){
			args$col.facet.y = "."
		}
		if(is.null(args$facet.scales)){
			args$facet.scales = "fixed"
		}
		p <- p + facet_grid(as.formula(paste(args$col.facet.x,args$col.facet.y, sep=" ~ ")), scales = args$facet.scales)
	}
	return(p)
}



##############################################################################################################
## HISTOGRAM
##############################################################################################################
plotting.histogram = function(data, args){
	## bin width 
	if(is.null(args$binwidth) || args$binwidth<=0){
		args$binwidth = NULL
	}

	## create plot
	p <- ggplot(data, aes_string(x=args$col.x, color=args$col.color, fill=args$col.fill))

	## histogram
	if(args$dens){
		p <- p + geom_histogram(aes(y=..density..))
		if(args$dens.cur){
			if(is.null(args$dens.color)){
				args$dens.color = "black"
			}
			p <- p +  geom_density(alpha=.2, color=args$dens.color)
		}
	}else{
		p <- p + geom_histogram(binwidth=args$binwidth)
	}

	## add facets
	p = plotting.addFacets(p, args)

	## add scales, labels and guides
	p = plotting.addScalesAndLabelsAndGuides(p, args)

	## change layout
	p <- p + geom_default(legend=="vertical")
	
	return(p)
}

##############################################################################################################
## SCATTERPLOT
##############################################################################################################
plotting.scatterplot = function(data, args){
	## create plot
	p <- ggplot(data, aes_string(x=args$col.x, y=args$col.y, fill=args$col.fill, color=args$col.color, shape=args$col.shape)) +
		geom_point(na.rm = TRUE, alpha=args$alpha, size=args$size) 

	## add facets
	p = plotting.addFacets(p, args)
	## add density if matrix
	if( !is.null(args$matrix)){
		p <- p + stat_density(aes(x = x, y = ..scaled.. * diff(range(x)) + min(x)), data = data.pairs$densities, position = "identity",colour = "grey20", geom = "line")
	}

	## add scales, labels and guides
	p = plotting.addScalesAndLabelsAndGuides(p, args)

	## change layout
	p <- p + geom_default(legend=="vertical")
	
	return(p)
}

##############################################################################################################
## BOXPLOT
##############################################################################################################
plotting.boxplot = function(data, args){
	## create plot
	p <- ggplot(data, aes_string(x=args$col.x, y=args$col.y, fill=args$col.fill))

	## create boxplot
	#<outliers, no, all, all jittered>
	if(args$points == "outliers"){
		p <- p + geom_boxplot() ## todo add shape and or color
	}else{
		p <- p + geom_boxplot(outlier.shape = NA) 
		if(args$points == "all" | args$points == "all jittered"){
			position="identity"
			if(args$points == "all jittered"){
				position = "jittered"
			}
			p <- p + geom_point(aes_string(color=args$col.color, shape=args$col.shape), position=position) ## todo add she and color
		}
	}

	## add facets
	p = plotting.addFacets(p, args)

	## add scales, labels and guides
	p = plotting.addScalesAndLabelsAndGuides(p, args)

	## change layout
	p <- p + geom_default(legend=="vertical")

	
	return(p)
}



##############################################################################################################
## BARPLOT
##############################################################################################################
plotting.barplot = function(data, args){
	## create plot	
	p <- ggplot(data, aes_string(x=args$col.x, y=args$col.y, color=args$col.color, fill=args$col.fill)) + 
		geom_bar(stat = "identity")


	## add facets
	p = plotting.addFacets(p, args)

	## add scales, labels and guides
	p = plotting.addScalesAndLabelsAndGuides(p, args)

	## change layout
	p <- p + geom_default(legend=="vertical")
	return(p)
}