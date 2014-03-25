#' Function to visualise input data matrix using heatmap
#'
#' \code{visHeatmap} is supposed to visualise input data matrix using heatmap. Note: this heatmap displays matrix in a bottom-to-top direction
#'
#' @param data an input gene-sample data matrix used for heatmap
#' @param scale a character indicating when the input matrix should be centered and scaled. It can be one of "none" (no scaling), "row" (being scaled in the row direction), "column" (being scaled in the column direction)
#' @param row.metric distance metric used to calculate the distance metric between rows. It can be one of "none" (i.e. no dendrogram between rows), "pearson", "spearman", "kendall", "euclidean", "manhattan", "cos" and "mi". See details at  \url{http://suprahex.r-forge.r-project.org/sDistance.html}
#' @param row.method the agglomeration method used to cluster rows. This should be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". See 'Note' below for details
#' @param column.metric distance metric used to calculate the distance metric between columns. It can be one of "none" (i.e. no dendrogram between rows), "pearson", "spearman", "kendall", "euclidean", "manhattan", "cos" and "mi". See details at  \url{http://suprahex.r-forge.r-project.org/sDistance.html}
#' @param column.method the agglomeration method used to cluster columns. This should be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". See 'Note' below for details
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum z/patttern values for which colors should be plotted, defaulting to the range of the finite values of z. Each of the given colors will be used to color an equispaced interval of this range. The midpoints of the intervals cover the range, so that values just outside the range will be plotted
#' @param row.cutree an integer scalar specifying the desired number of groups being cut from the row dendrogram. Note, this optional is only enabled when the row dengrogram is built
#' @param row.colormap short name for the colormap to color-code the row groups (i.e. sidebar colors used to annotate the rows)
#' @param column.cutree an integer scalar specifying the desired number of groups being cut from the column dendrogram. Note, this optional is only enabled when the column dengrogram is built
#' @param column.colormap short name for the colormap to color-code the column groups (i.e. sidebar colors used to annotate the columns)
#' @param ... additional graphic parameters. Type ?heatmap for the complete list.
#' @return 
#' invisible
#' @note The clustering methods are provided:
#' \itemize{
#' \item{"ward": Ward's minimum variance method aims at finding compact, spherical clusters}
#' \item{"single": The single linkage method (which is closely related to the minimal spanning tree) adopts a 'friends of friends' clustering strategy}
#' \item{"complete": The complete linkage method finds similar clusters}
#' \item{"average","mcquitty","median","centroid": These methods can be regarded as aiming for clusters with characteristics somewhere between the single and complete link methods. Two methods "median" and "centroid" are not leading to a monotone distance measure, or equivalently the resulting dendrograms can have so called inversions (which are hard to interpret)}
#' }
#' @export
#' @seealso \code{\link{visHeatmap}}
#' @include visHeatmap.r
#' @examples
#' # 1) generate data with an iid matrix of 100 x 9
#' data <- cbind(matrix(rnorm(100*3,mean=0,sd=1), nrow=100, ncol=3), 
#' matrix(rnorm(100*3,mean=0.5,sd=1), nrow=100, ncol=3), 
#' matrix(rnorm(100*3,mean=-0.5,sd=1), nrow=100, ncol=3))
#' colnames(data) <- c("S1","S1","S1","S2","S2","S2","S3","S3","S3")
#'
#' # 2) prepare colors for the column sidebar
#' lvs <- unique(colnames(data))
#' lvs_color <- visColormap(colormap="rainbow")(length(lvs))
#' my_ColSideColors <- sapply(colnames(data), function(x) lvs_color[x==lvs])
#'
#' # 3) heatmap with row dendrogram (with 10 color-coded groups)
#' visHeatmap(data, row.metric="euclidean", row.method="average", colormap="gbr", zlim=c(-2,2), 
#' ColSideColors=my_ColSideColors, row.cutree=10, row.colormap="jet", labRow=NA)

visHeatmap <- function (data, scale=c("none","row","column"), row.metric=c("none","pearson","spearman","kendall","euclidean","manhattan","cos","mi"), row.method=c("ward","single","complete","average","mcquitty","median","centroid"), column.metric=c("none","pearson","spearman","kendall","euclidean","manhattan","cos","mi"), column.method=c("ward","single","complete","average","mcquitty","median","centroid"),
colormap=c("bwr","jet","gbr","wyr","br","yr","rainbow","wb"), ncolors=64, zlim=NULL, row.cutree=NULL, row.colormap=c("rainbow"), column.cutree=NULL, column.colormap=c("rainbow"), ...)
{

    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
    }else if(is.null(data)){
        stop("The input data must be matrix.\n")
    }
    
    row.metric <- match.arg(row.metric)
    row.method <- match.arg(row.method)
    column.metric <- match.arg(column.metric)
    column.method <- match.arg(column.method)
    scale <- match.arg(scale)
    
    ## for rows
    my_RowSideColors <- NULL
    if(row.metric!="none"){
        row_distance <- as.dist(sDistance(data, metric=row.metric))
        row_cluster <- hclust(row_distance, method=row.method)
        my_Rowv <- as.dendrogram(row_cluster)
        
        if(!is.null(row.cutree)){
            if(row.cutree==as.integer(row.cutree) & row.cutree>=2 & row.cutree<=nrow(data)){
                row_cutree <- cutree(row_cluster, k=row.cutree)
                row_color <- visColormap(colormap=row.colormap)(length(unique(row_cutree)))
                my_RowSideColors <- row_color[as.vector(row_cutree)]
            }
        }
    }else{
        my_Rowv <- NA
    }
    ## for columns 
    my_ColSideColors <- NULL
    if(column.metric!="none"){
        column_distance <- as.dist(sDistance(t(data), metric=column.metric))
        column_cluster <- hclust(column_distance, method=column.method)
        my_Colv <- as.dendrogram(column_cluster)
        
        if(!is.null(column.cutree)){
            if(column.cutree==as.integer(column.cutree) & column.cutree>=2 & column.cutree<=ncol(data)){
                column_cutree <- cutree(column_cluster, k=column.cutree)
                column_color <- visColormap(colormap=column.colormap)(length(unique(column_cutree)))
                my_ColSideColors <- column_color[as.vector(column_cutree)]
            }
        }
    }else{
        my_Colv <- NA
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

    mat_data <- data
    mat_data[mat_data < zlim[1]] <- zlim[1]
    mat_data[mat_data > zlim[2]] <- zlim[2]
    
    if(!is.null(my_RowSideColors) & !is.null(my_ColSideColors)){
        heatmap(as.matrix(mat_data),col=visColormap(colormap=colormap)(ncolors),zlim=zlim,Rowv=my_Rowv,Colv=my_Colv, scale=scale, RowSideColors=my_RowSideColors, ColSideColors=my_ColSideColors,  ...)
    }else if(!is.null(my_RowSideColors) & is.null(my_ColSideColors)){
        heatmap(as.matrix(mat_data),col=visColormap(colormap=colormap)(ncolors),zlim=zlim,Rowv=my_Rowv,Colv=my_Colv, scale=scale, RowSideColors=my_RowSideColors,  ...)
    }else if(is.null(my_RowSideColors) & !is.null(my_ColSideColors)){
        heatmap(as.matrix(mat_data),col=visColormap(colormap=colormap)(ncolors),zlim=zlim,Rowv=my_Rowv,Colv=my_Colv, scale=scale, ColSideColors=my_ColSideColors,  ...)
    }else{
        heatmap(as.matrix(mat_data),col=visColormap(colormap=colormap)(ncolors),zlim=zlim,Rowv=my_Rowv,Colv=my_Colv, scale=scale,  ...)
    }
    
    invisible()
}