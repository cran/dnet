#' Function to visualise input data matrix using advanced heatmap
#'
#' \code{visHeatmapAdv} is supposed to visualise input data matrix using advanced heatmap. It allows for adding multiple sidecolors in both columns and rows. Besides, the sidecolor can be automatically added via cutting histogram into groups. Note: this heatmap displays matrix in a top-to-bottom direction
#'
#' @param data an input gene-sample data matrix used for heatmap
#' @param scale a character indicating when the input matrix should be centered and scaled. It can be one of "none" (no scaling), "row" (being scaled in the row direction), "column" (being scaled in the column direction)
#' @param Rowv determines if and how the row dendrogram should be reordered. By default, it is TRUE, which implies dendrogram is computed and reordered based on row means. If NULL or FALSE, then no dendrogram is computed and no reordering is done. If a dendrogram, then it is used "as-is", ie without any reordering. If a vector of integers, then dendrogram is computed and reordered based on the order of the vector
#' @param Colv determines if and how the column dendrogram should be reordered.	Has the options as the Rowv argument above and additionally when x is a square matrix, Colv = "Rowv" means that columns should be treated identically to the rows
#' @param dendrogram character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms. Defaults to 'both'. However, if Rowv (or Colv) is FALSE or NULL and dendrogram is 'both', then a warning is issued and Rowv (or Colv) arguments are honoured
#' @param dist.metric distance metric used to calculate the distance metric between columns (or rows). It can be one of "none" (i.e. no dendrogram between rows), "pearson", "spearman", "kendall", "euclidean", "manhattan", "cos" and "mi". See details at  \url{http://suprahex.r-forge.r-project.org/sDistance.html}
#' @param linkage.method the agglomeration method used to cluster/linkages columns (or rows). This should be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". See 'Note' below for details
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param RowSideColors NULL or a matrix of "numRowsidebars" X nrow(x), where "numRowsidebars" stands for the number of sidebars annotating rows of x. This matrix contains the color names for vertical sidebars. By default, it sets to NULL. In this case, sidebars in rows can still be enabled by cutting the row dendrogram into several clusters (see the next two parameters)
#' @param row.cutree an integer scalar specifying the desired number of groups being cut from the row dendrogram. Note, this optional is only enabled when the ColSideColors is NULL
#' @param row.colormap short name for the colormap to color-code the row groups (i.e. sidebar colors used to annotate the rows)
#' @param ColSideColors NULL or a matrix of ncol(x) X "numColsidebars", where "numColsidebars" stands for the number of sidebars annotating the columns of x. This matrix contains the color names for horizontal sidebars. By default, it sets to NULL. In this case, sidebars in columns can still be enabled by cutting the column dendrogram into several clusters (see the next two parameters)
#' @param column.cutree an integer scalar specifying the desired number of groups being cut from the column dendrogram. Note, this optional is only enabled when the column dengrogram is built
#' @param column.colormap short name for the colormap to color-code the column groups (i.e. sidebar colors used to annotate the columns)
#' @param ... additional graphic parameters. For most parameters, please refer to \url{http://www.inside-r.org/packages/cran/gplots/docs/heatmap.2}. For example, the parameters "srtRow" and "srtCol" to control the angle of row/column labels (in degrees from horizontal: 45 degrees for the column, 0 degrees for the row, by default), i.e. string rotation. The parameters "offsetRow" and "offsetCol" to indicate the number of character-width spaces to place between row/column labels and the edge of the plotting region. Unique to this function, there are two parameters "RowSideWidth" and RowSideLabelLocation, to respectively indicate the fraction of the row side width and the location (either bottom or top) of the row side labelling; the other two parameters "ColSideHeight" and "ColSideLabelLocation" for the column side height and the location (either left or right) of the column side labelling; and two parameters "RowSideBox" and "ColSideBox" to indicate whether there are boxes outside. 
#' @return 
#' invisible
#' @note The clustering/linkage methods are provided:
#' \itemize{
#' \item{"ward": Ward's minimum variance method aims at finding compact, spherical clusters}
#' \item{"single": The single linkage method (which is closely related to the minimal spanning tree) adopts a 'friends of friends' clustering strategy}
#' \item{"complete": The complete linkage method finds similar clusters}
#' \item{"average","mcquitty","median","centroid": These methods can be regarded as aiming for clusters with characteristics somewhere between the single and complete link methods. Two methods "median" and "centroid" are not leading to a monotone distance measure, or equivalently the resulting dendrograms can have so called inversions (which are hard to interpret)}
#' }
#' @export
#' @seealso \code{\link{visHeatmapAdv}}
#' @include visHeatmapAdv.r
#' @examples
#' # 1) generate data with an iid matrix of 100 x 9
#' data <- cbind(matrix(rnorm(100*3,mean=0,sd=1), nrow=100, ncol=3), 
#' matrix(rnorm(100*3,mean=0.5,sd=1), nrow=100, ncol=3), 
#' matrix(rnorm(100*3,mean=-0.5,sd=1), nrow=100, ncol=3))
#' colnames(data) <- c("S1_R1","S1_R2","S1_R3","S2_R1","S2_R2","S2_R3","S3_R1","S3_R2","S3_R3")
#'
#' # 2) heatmap after clustering both rows and columns
#' # 2a) shown with row and column dendrograms
#' visHeatmapAdv(data, dendrogram="both", colormap="gbr", zlim=c(-2,2), KeyValueName="log2(Ratio)",
#' add.expr=abline(v=(1:(ncol(data)+1))-0.5,col="white"), 
#' lmat=rbind(c(4,3), c(2,1)), lhei=c(1,5), lwid=c(1,3))
#' # 2b) shown with row dendrogram only
#' visHeatmapAdv(data, dendrogram="row", colormap="gbr", zlim=c(-2,2))
#' # 2c) shown with column dendrogram only
#' visHeatmapAdv(data, dendrogram="column", colormap="gbr", zlim=c(-2,2))
#'
#' # 3) heatmap after only clustering rows (with 2 color-coded groups)
#' visHeatmapAdv(data, Colv=FALSE, colormap="gbr", zlim=c(-2,2), 
#' row.cutree=2, row.colormap="jet", labRow=NA)
#'
#' # 4) prepare colors for the column sidebar
#' # color for stages (S1-S3)
#' stages <- sub("_.*","",colnames(data))
#' lvs <- unique(stages)
#' lvs_color <- visColormap(colormap="rainbow")(length(lvs))
#' col_stages <- sapply(stages, function(x) lvs_color[x==lvs])
#' # color for replicates (R1-R3)
#' replicates <- sub(".*_","",colnames(data))
#' lvs <- unique(replicates)
#' lvs_color <- visColormap(colormap="rainbow")(length(lvs))
#' col_replicates <- sapply(replicates, function(x) lvs_color[x==lvs])
#' # combine both color vectors
#' ColSideColors <- cbind(col_stages,col_replicates)
#' colnames(ColSideColors) <- c("Stages","Replicates")
#'
#' # 5) heatmap without clustering on rows and columns but with the two sidebars in columns
#' visHeatmapAdv(data, Rowv=FALSE, Colv=FALSE, colormap="gbr", zlim=c(-2,2), 
#' density.info="density", tracecol="yellow", ColSideColors=ColSideColors, 
#' ColSideHeight=0.5, ColSideLabelLocation="right")

visHeatmapAdv <- function (data, scale=c("none","row","column"), Rowv=T, Colv=T, dendrogram=c("both","row","column","none"), dist.metric=c("euclidean","pearson","spearman","kendall","manhattan","cos","mi"), linkage.method=c("complete","ward","single","average","mcquitty","median","centroid"), 
colormap=c("bwr","jet","gbr","wyr","br","yr","rainbow","wb"), ncolors=64, zlim=NULL, RowSideColors=NULL, row.cutree=NULL, row.colormap=c("jet"), ColSideColors=NULL, column.cutree=NULL, column.colormap=c("jet"), ...)
{

    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
    }else if(is.null(data)){
        stop("The input data must be matrix.\n")
    }
    
    scale <- match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    dist.metric <- match.arg(dist.metric)
    linkage.method <- match.arg(linkage.method)
    
    ## make sure RowSideColors and ColSideColors are always matrix
    if (is.vector(RowSideColors)){
        RowSideColors <- matrix(RowSideColors, nrow=1)
    }else if(is.matrix(RowSideColors) | is.data.frame(RowSideColors)){
        RowSideColors <- as.matrix(RowSideColors)
    }
    if (is.vector(ColSideColors)){
        ColSideColors <- matrix(ColSideColors, ncol=1)
    }else if(is.matrix(ColSideColors) | is.data.frame(ColSideColors)){
        ColSideColors <- as.matrix(ColSideColors)
    }
    
    ## for rows
    if(Rowv==TRUE){
        if(!is.null(row.cutree)){
            if(is.null(RowSideColors)){
                if(row.cutree==as.integer(row.cutree) & row.cutree>=2 & row.cutree<=nrow(data)){
                    distance <- as.dist(sDistance(data, metric=dist.metric))
                    cluster <- hclust(distance, method=linkage.method)
                    row_cutree <- cutree(cluster, k=row.cutree)
                    row_color <- visColormap(colormap=row.colormap)(length(unique(row_cutree)))
                    RowSideColors <- matrix(row_color[as.vector(row_cutree)], nrow=1)
                }
            }
        }
    }
    ## for columns
    if(Colv==TRUE){
        if(!is.null(column.cutree)){
            if(is.null(ColSideColors)){
                if(column.cutree==as.integer(column.cutree) & column.cutree>=2 & column.cutree<=ncol(data)){
                    distance <- as.dist(sDistance(t(data), metric=dist.metric))
                    cluster <- hclust(distance, method=linkage.method)
                    column_cutree <- cutree(cluster, k=column.cutree)
                    column_color <- visColormap(colormap=column.colormap)(length(unique(column_cutree)))
                    ColSideColors <- matrix(column_color[as.vector(column_cutree)], ncol=1)
                }
            }
        }
    }

    ## determine the color range
    vmin <- floor(quantile(data, 0.05))
    vmax <- ceiling(quantile(data, 0.95))
    if(vmin < 0 & vmax > 0){
        vsym <- abs(min(vmin, vmax))
        vmin <- -1*vsym
        vmax <- vsym
    }
    if(!is.null(zlim)){
        if(zlim[1] < floor(min(data)) | zlim[2] > ceiling(max(data))){
            #zlim <- c(vmin,vmax)
        }
    }else{
        zlim <- c(vmin,vmax)
    }

    ## Defining breaks for the color scale
    myBreaks <- seq(zlim[1], zlim[2], length.out=ncolors+1)


################################################################################################
################################################################################################
heatmap.2 <- function(x,
                      
                      # dendrogram control
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      
                      # data scaling
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      
                      # image plot
                      revC = identical(Colv,"Rowv"),
                      
                      # expression that will be evaluated after the call to image. Can be used to add components to the plot.
                      add.expr, # add.expr=abline(h=(1:(nrow(x)+1))-0.5,v=(1:(ncol(x)+1))-0.5, col="white", lwd=3, lty=1)
                      
                      # mapping data to colors
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      
                      # block sepration: vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      
                      # cell labeling: (optional) matrix of character strings which will be placed within each color cell, e.g. p-value symbols
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      
                      # level trace
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(7,5), # controling the margins (1st for the bottom and 2nd for the right)
                      
                      # Row/Column Labeling/annotating
                      RowSideColors,
                      ColSideColors,
                      RowSideBox = TRUE, # whether adding box for row sides
                      ColSideBox = TRUE, # whether adding box for column sides
                      RowSideWidth = 0.25, # fraction of the width of row sides
                      ColSideHeight = 0.25, # fraction of the height of column sides
                      RowSideLabelLocation = c("top","bottom"), # location of column side labelling (either at bottom or top)
                      ColSideLabelLocation = c("right","left"), # location of row side labelling (either at left or right)
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      srtRow = 0, # angle of row labels, in degrees from horizontal
                      srtCol = 45, # angle of column labels, in degrees from horizontal
                      offsetRow = 0, # Number of character-width spaces to place between row labels and the edge of the plotting region
                      offsetCol = 0, # Number of character-width spaces to place between column labels and the edge of the plotting region
                      adjRow = c(0,NA), # 2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation)
                      adjCol = c(NA,0), # 2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation)
                      
                      # color key + density info
                      ## logical indicating whether a color-key should be shown
                      key = TRUE,
                      ## numeric value indicating the size of the key
                      keysize = 1.5,
                      ## character string indicating whether to superimpose "histogram","density","none"
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      ## Boolean indicating whether the color key should be made symmetric about 0. Defaults to TRUE if the data includes negative values, and to FALSE otherwise
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      ## Numeric scaling value for tuning the kernel width when a density plot is drawn on the color key
                      densadj = 0.25,
                      
                      # plot labels
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      
                      # visual layout: position matrix, column height, column width
                      # always 1 for heatmap, 2 for row dendrogram, 3 for col dendrogram, 4 for color key. See: layout.show(layout(lmat, widths=lwid, heights=lhei))
                      lmat = NULL, # lmat=rbind(c(4,3,0), c(2,1,0), c(0,0,0))
                      lhei = NULL, # lhei=c(1,2,1)
                      lwid = NULL, # lwid=c(1,2,1)
                      
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName = "Value",
                      ...)
{
    ########## Add RowSideLabelLocation and ColSideLabelLocation for the location of labelling
    RowSideLabelLocation <- match.arg(RowSideLabelLocation)
    ColSideLabelLocation <- match.arg(ColSideLabelLocation)
    ##########
     
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }
 
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
 
        if (!missing(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], ColSideHeight*NumColSideColors, lhei[2])
        }
 
        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], RowSideWidth*NumRowSideColors, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
 
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
 
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
 
    if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors)){
        if (nrow(RowSideColors)==0){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
                ## add box
                if(RowSideBox==TRUE) box(lwd=1,col="black")
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            ## add box
            if(RowSideBox==TRUE){
                if(nrow(RowSideColors) >= 2){
                    tmp <- dim(rsc)[2]-1
                    abline(v=(((0:(tmp+1))-0.5)/max(1,tmp)), lwd=1, col="black")
                }
                box(lwd=1, col="black")
            }
            if (length(rownames(RowSideColors)) > 0) {
                if(RowSideLabelLocation=="bottom"){
                    axis(1, 0:(dim(rsc)[2] - 1)/(max(1,dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE, line=-0.8)
                }else if(RowSideLabelLocation=="top"){
                    axis(3, 0:(dim(rsc)[2] - 1)/(max(1,dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE, line=-0.8)
                }
            }
        }
    }
 
    if (!missing(ColSideColors)) {
 
        if (ncol(ColSideColors)==0){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
            ## add box
            if(ColSideBox==TRUE) box(lwd=1,col="black")
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            ##### add box
            if(ColSideBox==TRUE){
                if(ncol(ColSideColors) >= 2){
                    tmp <- dim(csc)[2]-1
                    abline(h=(((0:(tmp+1))-0.5)/max(1,tmp)), lwd=1, col="black")
                }
                box(lwd=1, col="black")
            }
            #####
            if (length(colnames(ColSideColors)) > 0) {
                if(ColSideLabelLocation=="left"){
                    axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE, line=-0.8)
                }else if(ColSideLabelLocation=="right"){
                    axis(4, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE, line=-0.8)
                }
            }
        }
    }
 
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    ###########
    box(lwd=1, col="black")
    ###########
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",col = na.color, add = TRUE)
        ###########
        box(lwd=1, col="black")
        ###########
    }
    
    ############ Add srtCol, offsetCol and adjCol according to heatmap.2 from gplots
    if (is.null(srtCol)) 
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
            padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    ############
    ##axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
    if (!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25)
    
    ############ Add srtRow, offsetRow and adjRow according to heatmap.2 from gplots
    if (is.null(srtRow)) {
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
            tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
                srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    ############
    ##axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
    if (!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25)
    
    if (!missing(add.expr))
        #eval(substitute(add.expr))
        #eval(parse(text=add.expr))
        eval(add.expr)
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
 
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}
################################################################################################

    if(dendrogram=="both"){
        if(Rowv==TRUE & Colv==FALSE){
            dendrogram <- "row"
        }else if(Rowv==FALSE & Colv==TRUE){
            dendrogram <- "column"
        }else if(Rowv==FALSE & Colv==FALSE){
            dendrogram <- "none"
        }
    }
    
    if(!is.null(RowSideColors) & !is.null(ColSideColors)){
        heatmap.2(as.matrix(data), scale=scale, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, distfun=function(x) as.dist(sDistance(x, metric=dist.metric)), hclustfun=function(x) hclust(x, method=linkage.method), col=visColormap(colormap=colormap)(ncolors), breaks=myBreaks, zlim=zlim,RowSideColors=RowSideColors,ColSideColors=ColSideColors,...)
    }else if(!is.null(RowSideColors) & is.null(ColSideColors)){
        heatmap.2(as.matrix(data), scale=scale, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, distfun=function(x) as.dist(sDistance(x, metric=dist.metric)), hclustfun=function(x) hclust(x, method=linkage.method), col=visColormap(colormap=colormap)(ncolors), breaks=myBreaks, zlim=zlim,RowSideColors=RowSideColors,...)
    }else if(is.null(RowSideColors) & !is.null(ColSideColors)){
        heatmap.2(as.matrix(data), scale=scale, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, distfun=function(x) as.dist(sDistance(x, metric=dist.metric)), hclustfun=function(x) hclust(x, method=linkage.method), col=visColormap(colormap=colormap)(ncolors), breaks=myBreaks, zlim=zlim,ColSideColors=ColSideColors,...)
    }else{
        heatmap.2(as.matrix(data), scale=scale, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, distfun=function(x) as.dist(sDistance(x, metric=dist.metric)), hclustfun=function(x) hclust(x, method=linkage.method), col=visColormap(colormap=colormap)(ncolors), breaks=myBreaks, zlim=zlim,...)
    }
    
    invisible()
}

