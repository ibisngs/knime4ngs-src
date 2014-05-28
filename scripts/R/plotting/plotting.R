##############################################################################################################
## PLOTTING ARGS
##############################################################################################################
plotting.addArgs = function(parser){

	## GLOBALS 
	parser$add_argument("-g", "--globals"      , type="character", action="store"      , dest="file.global" , required=TRUE, help="path to globals file"        , metavar="<path>")
	## IN- AND OUTPUT-FILES
	parser$add_argument( "--data"              , type="character", action="store"      , dest="file.in"     , required=TRUE, help="path to input data file (cols=variables, rows=observations)", metavar="<path>")

	## ARGUMENTS
	parser$add_argument("-x" , "--col.x"       , type="character", action="store"      , dest="col.x"                          , help="name of column to be plotted on x-axis (defining different boxes)" , metavar="<column name>")
	parser$add_argument("-y" , "--col.y"       , type="character", action="store"      , dest="col.y"                          , help="name of column to be plotted on y-axis (defining different boxes)" , metavar="<column name>")
	# color
	parser$add_argument("-c" , "--col.color"   , type="character", action="store"      , dest="col.color"                      , help="name of column which defines color"         , metavar="<column name>")
	parser$add_argument("-cl", "--lab.color"   , type="character", action="store"      , dest="label.color"                    , help="label for color legend"         , metavar="<column name>")
	parser$add_argument("-cs", "--scale.color" , type="character", action="store"      , dest="scale.color"                    , help="scale palette of color"         , metavar="<string>")
	# fill color
	parser$add_argument("-f" , "--col.fill"    , type="character", action="store"      , dest="col.fill"                       , help="name of column which defines fill"          , metavar="<column name>")
	parser$add_argument("-fl", "--lab.fill"    , type="character", action="store"      , dest="label.fill"                     , help="label for fill color legend"          , metavar="<column name>")
	parser$add_argument("-fs", "--scale.fill"  , type="character", action="store"      , dest="scale.fill"                     , help="scale palette of fill color"         , metavar="<string>")
	# shape
	parser$add_argument("-s" , "--col.shape"   , type="character", action="store"      , dest="col.shape"                      , help="name of column which defines shape"          , metavar="<column name>")
	parser$add_argument("-sl", "--lab.shape"   , type="character", action="store"      , dest="label.shape"                    , help="label for shape legend"          , metavar="<string>")
	# facet
	parser$add_argument("-fx", "--facet.x"     , type="character", action="store"      , dest="col.facet.x"                    , help="name of column which defines x axis of facet grid"         , metavar="<column name>")
	parser$add_argument("-fy", "--facet.y"     , type="character", action="store"      , dest="col.facet.y"                    , help="name of column which defines x axis of facet grid"         , metavar="<column name>")
	# labels
	parser$add_argument("-t" , "--title"       , type="character", action="store"      , dest="label.plot"                     , help="Title of plot"   , metavar="<String>")
	parser$add_argument("-xl", "--lab.x"       , type="character", action="store"      , dest="label.x"                        , help="Label of x axis" , metavar="<String>")
	parser$add_argument("-yl", "--lab.y"       , type="character", action="store"      , dest="label.y"                        , help="Label of y axis" , metavar="<String>")
	# output
	parser$add_argument("-w" , "--width"       , type="integer"  , action="store"      , dest="width"         , default=800    , help="width of image"  , metavar="<int>")
	parser$add_argument(       "--height"      , type="integer"  , action="store"      , dest="height"        , default=600    , help="height of image" , metavar="<int>")
	parser$add_argument("-i" , "--image"       , type="character", action="store"      , dest="image"         , required=TRUE  , help="name of png file" , metavar="<path>")
}


    


##############################################################################################################
## PLOT IMAGE
##############################################################################################################
plotting.print = function(p, args){
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
	loadLib("ggplot2")
	loadLib("grid")
	
	props = theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

	## LEGEND
	if( !is.na(legend)){
		if(legend=="vertical"){
			props <- props + theme(legend.position=c(0,1), legend.justification=c(0,1), legend.direction="vertical", legend.box="vertical")
		}
	}
# 	
	return(props)
}


##############################################################################################################
## ALIGN MULTIPLE LEGENDS
##############################################################################################################
multiLegendAlign <- function(p, align="left"){
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
plotting.checkArgs = function(args){
	## color
	if( is.null(args$col.color) || is.na(args$col.color)|| args$col.color=="NULL" || args$col.color=="" || args$col.color=="NA" ){
		args$col.color <- NULL
		args$lab.color <- NULL
	}else{
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
		if(is.null(args$lab.shape) || args$lab.shape == "NULL"){
			args$lab.shape = NULL
		}else if(args$lab.shape == "NA" || args$lab.shape == ""){
			args$lab.shape = NA
		}
	}

	## X Label
	if( is.null(args$label.x) || args$label.x == "" || args$label.x == "NA" || args$label.x == "NULL" ){
		args$label.x = NULL
	}

	## Y Label
	if( is.null(args$label.y) || args$label.y == "" || args$label.y == "NA" || args$label.y == "NULL" ){
		args$label.y = NULL
	}

	## Main Label
	if( is.null(args$label.plot) || args$label.plot == "" || args$label.plot == "NA" || args$label.plot == "NULL" ){
		args$label.plot = NULL
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
	p <- p + labs(x=args$label.x, y=args$label.y, title=args$label.plot, color=args$label.color, shape=args$label.shape, fill=args$label.fill)

	## guides
	if(is.null(args$label.color)){
		p <- p + guides(color="none")
	}
	if(is.null(args$label.fill)){
		p <- p + guides(fill="none")
	}
	if(is.null(args$label.shape)){
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
