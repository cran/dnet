#' Function to view enrichment results of dEnricher
#'
#' \code{dEnricherView} is supposed to view results of enrichment analysis by dEnricher. 
#'
#' @param eTerm an object of class "eTerm"
#' @param top_num the maximum number of gene sets (terms) will be viewed
#' @param sortBy which statistics will be used for sorting and viewing gene sets (terms). It can be "adjp" for adjusted p value, "pvalue" for p value, "zscore" for enrichment z-score, "nSet" for the number of sets (terms), "nOverlap" for the number in overlaps, and "none" for ordering according to ID of gene sets (terms)
#' @param decreasing logical to indicate whether to sort in a decreasing order. If it is null, it would be true for "zscore", "nSet" or "nOverlap"; otherwise it would be false
#' @param details logical to indicate whether the detail information of gene sets (terms) is also viewed. By default, it sets to false for no inclusion
#' @return
#' a data frame with following components:
#' \itemize{
#'  \item{\code{setID}: term ID}
#'  \item{\code{nSet}: number of sets (terms)}
#'  \item{\code{nOverlap}: number in overlaps}
#'  \item{\code{zscore}: enrichment z-score}
#'  \item{\code{pvalue}: nominal p value}
#'  \item{\code{adjp}: adjusted p value}
#'  \item{\code{name}: term name; optional, it is only appended when "details" is true}
#'  \item{\code{namespace}: term namespace; optional, it is only appended when "details" is true}
#'  \item{\code{distance}: term distance; optional, it is only appended when "details" is true}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dEnricher}}
#' @include dEnricherView.r
#' @examples
#' \dontrun{
#' dEnricherView(eTerm, top_num=10, sortBy="adjp", decreasing=FALSE, details=TRUE)
#' dEnricherView(eTerm, top_num=NULL, sortBy="none")
#' }

dEnricherView <- function(eTerm, top_num=10, sortBy=c("adjp","pvalue","zscore","nSet","nOverlap","none"), decreasing=NULL, details=F) 
{

    if (class(eTerm) != "eTerm" ){
        stop("The function must apply to a 'eTerm' object.\n")
    }
    
    sortBy <- match.arg(sortBy)
    
    if( is.null(top_num) | length(top_num > length(eTerm$set_info$setID))==0 ){
        top_num <- length(eTerm$set_info$setID)
    }
    
    tab <- data.frame(setID         = eTerm$set_info$setID,
                       nSet         = sapply(eTerm$gs,length),
                       nOverlap     = sapply(eTerm$overlap,length),
                       zscore       = eTerm$zscore,
                       pvalue       = eTerm$pvalue,
                       adjp         = eTerm$adjp,
                       name         = eTerm$set_info$name,
                       namespace    = eTerm$set_info$namespace,
                       distance     = eTerm$set_info$distance
                      )
    
    
    if(details == T){
        res <- tab[,c(1:9)]
    }else{
        res <- tab[,c(1:6)]
    }
    
    if(is.null(decreasing)){
        if(sortBy=="zscore" | sortBy=="nSet" | sortBy=="nOverlap"){
            decreasing <- T
        }else{
            decreasing <- F
        }
    }
    
    switch(sortBy, 
        adjp={res <- res[order(res[,6], decreasing=decreasing)[1:top_num],]},
        pvalue={res <- res[order(res[,5], decreasing=decreasing)[1:top_num],]},
        zscore={res <- res[order(res[,4], decreasing=decreasing)[1:top_num],]},
        nSet={res <- res[order(res[,2], decreasing=decreasing)[1:top_num],]},
        nOverlap={res <- res[order(res[,3], decreasing=decreasing)[1:top_num],]},
        none={res <- res[order(res[,1], decreasing=decreasing)[1:top_num],]}
    )
    
    res
}