#' Function to find the tip node(s) of a direct acyclic graph (DAG)
#'
#' \code{dDAGtip} is supposed to find the tip node(s) of a direct acyclic graph (DAG; an ontology). It return the name (i.e Term ID) of the tip node(s). 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @return 
#' \itemize{
#'  \item{\code{tip}: the tip name (i.e. Term ID)}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dDAGtip}}
#' @include dDAGtip.r
#' @examples
#' # 1) load HPPA as igraph object
#' data(ig.HPPA)
#' g <- ig.HPPA
#'
#' # 2) find tips
#' tips <- dDAGtip(g)
#' tips

dDAGtip <- function (g)
{

    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    igr <- dDAGreverse(ig)
    tips <- dDAGroot(igr)
    
    return(tips)
}