#' Function to setup the pipeline for finding maximum-scoring module from an input graph and the signficance imposed on its nodes
#'
#' \code{dNetPipeline} is supposed to finish ab inito maximum-scoring module identification for the input graph with the node information on the significance (p-values). It returns an object of class "igraph" or "graphNEL". 
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param pval a vector containing input p-values. For each element, it must have the name that could be mapped onto the input graph. Also, the names in input "pval" should contain all those in the input graph "g", but the reverse is not necessary
#' @param method the method used for the transformation. It can be either "pdf" for the method based on the probability density function of the fitted model, or "cdf" for the method based on the cumulative distribution function of the fitted model
#' @param fdr the given FDR threshold. By default, it is set to NULL, meaning there is no constraint. If given, those p-values with the FDR below this are considered significant and thus scored positively. Instead, those p-values with the FDR above this given FDR are considered insigificant and thus scored negatively
#' @param nsize the desired number of nodes constrained to the resulting module. It is not nulll, a wide range of FDR will be scanned to find the FDR threshold leading to the desired number of nodes in the resulting module. Notably, the given FDR threshold will be overwritten. 
#' @param plot logical to indicate whether the histogram plot, contour plot and scatter plot should be drawn. By default, it sets to false for no plotting
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' a module with a maximum score, an object of class "igraph" or "graphNEL"
#' @note The pipeline sequentially consists of: 
#' \itemize{
#' \item{i) \code{\link{dBUMfit}} used to fit the p-value distribution under beta-uniform mixture model.}
#' \item{ii) if there is the desired number of nodes constrained to the resulting module, a wide range of FDR (including rough stage with large intervals, and finetune stage with smaller intervals) will be scanned to find the FDR threshold to meet the desired number of nodes.}
#' \item{iii) \code{\link{dBUMscore}} used to calculate the scores according to the fitted BUM and FDR threshold.}
#' \item{iv) \code{\link{dNetFind}} used to find maximum-scoring module from the input graph and scores imposed on its nodes.}
#' }
#' @export
#' @seealso \code{\link{dBUMfit}}, \code{\link{dBUMscore}}, \code{\link{dNetFind}}
#' @include dBUMfit.r
#' @examples
#' # 1) generate an vector consisting of random values from beta distribution
#' x <- rbeta(1000, shape1=0.5, shape2=1)
#' names(x) <- as.character(1:length(x))
#'
#' # 2) generate a random graph according to the ER model
#' g <- erdos.renyi.game(1000, 1/100)
#'
#' # 3) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#'
#' # 4) find maximum-scoring module based on fdr=0.1 threshold
#' module <- dNetPipeline(g=subg, pval=x, fdr=0.1)
#' 
#' \dontrun{
#' # 5) find maximum-scoring module with the desired node number nsize=20
#' # module <- dNetPipeline(g=subg, pval=x, nsize=20)
#' }

dNetPipeline <- function(g, pval, method=c("pdf","cdf","fdr"), fdr=NULL, nsize=NULL, plot=F, verbose=T)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    method <- match.arg(method)
    ####################
    
    if(method!="fdr"){
        if(verbose){
            message(sprintf("First, fit the p-value distribution under beta-uniform mixture model..."), appendLF=T)
        }
        ## fit a p-value distribution under beta-uniform mixture model
        if(plot){
            fit <- dBUMfit(pval, ntry=1, hist.bum=T, contour.bum=T, verbose=verbose)
        }else{
            fit <- dBUMfit(pval, ntry=1, hist.bum=F, contour.bum=F, verbose=verbose)
        }
    }else if(method=="fdr"){
        if(verbose){
            message(sprintf("First, consider the fdr distribution"), appendLF=T)
        }
    }
    
    if(verbose){
        message(sprintf("Second, determine FDR threshold..."), appendLF=T)
    }
    ## Determine the final fdr threshold
    if(is.null(nsize)){
        fdr_final <- fdr
    }else{
        fdr_final <- NULL
        ####################################
        if(verbose){
            message(sprintf("\tVia constraint on the size of subnetwork to be identified (%d nodes)", nsize), appendLF=T)
        }
        ## Constraint on the network size
        
        if(verbose){
            message(sprintf("\tScanning FDR at rough stage..."), appendLF=T)
        }
        ## at rough phase
        fdr_rough <- NULL
        nsize_rough <- 0
        for(i in seq(from=floor(log10(min(pval[pval!=0]))),to=0)){
            fdr_test <- 10^i
            
            if(method!="fdr"){
                scores_test <- dBUMscore(fit=fit, method=method, fdr=fdr_test, scatter.bum=F)
            }else if(method=="fdr"){
                scores_test <- dFDRscore(pval, fdr.threshold=fdr_test, scatter=F)
            }
    
            module_test <- suppressWarnings(dNetFind(g, scores_test))
            nsize_test <- vcount(module_test)
            
            if(verbose){
                message(sprintf("\t\tFDR: %1.2e, corresponding to the network size (%d nodes)", fdr_test, nsize_test), appendLF=T)
            }
            
            if(nsize_test >= nsize){
                fdr_rough <- fdr_test
                nsize_rough <- nsize_test
                break
            }
        }
        
        if(nsize_rough==nsize){
            fdr_final <- fdr_rough
        }else{
            if(verbose){
                message(sprintf("\tScanning FDR at finetune stage..."), appendLF=T)
            }
            ## at finetune phase
            fdr_final <- NULL
            for(fdr_test in seq(from=fdr_rough/10+fdr_rough/20,to=fdr_rough-fdr_rough/20,by=fdr_rough/20)){
            
                if(method!="fdr"){
                    scores_test <- dBUMscore(fit=fit, method=method, fdr=fdr_test, scatter.bum=F)
                }else if(method=="fdr"){
                    scores_test <- dFDRscore(pval, fdr.threshold=fdr_test, scatter=F)
                }
            
                module_test <- suppressWarnings(dNetFind(g, scores_test))
                nsize_test <- vcount(module_test)
            
                if(verbose){
                    message(sprintf("\t\tFDR: %1.2e, corresponding to the network size (%d nodes)", fdr_test, nsize_test), appendLF=T)
                }
            
                if(nsize_test >= nsize){
                    fdr_final <- fdr_test
                    break
                }
            }
        }
        
    }
    
    if(is.null(fdr_final)){
        if(verbose){
            message(sprintf("\tNo FDR threshold required"), appendLF=T)
        }
    }else{
        if(verbose){
            message(sprintf("\tFDR threshold: %1.2e", fdr_final), appendLF=T)
        }
    }
    
    if(method!="fdr"){
        if(verbose){
            message(sprintf("Third, calculate the scores according to the fitted BUM and FDR threshold (if any)..."), appendLF=T)
        }
        ## calculate the scores according to the fitted BUM and fdr threshold (fdr_final) 
        if(plot){
            scores <- dBUMscore(fit=fit, method=method, fdr=fdr_final, scatter.bum=T)
        }else{
            scores <- dBUMscore(fit=fit, method=method, fdr=fdr_final, scatter.bum=F)
        }
    }else if(method=="fdr"){
        if(verbose){
            message(sprintf("Third, calculate the scores according to the given fdr and the threshold (if any)..."), appendLF=T)
        }
        ## calculate the scores according to the fitted BUM and fdr threshold (fdr_final) 
        if(plot){
            scores <- dFDRscore(pval, fdr.threshold=fdr_final, scatter=T)
        }else{
            scores <- dFDRscore(pval, fdr.threshold=fdr_final, scatter=F)
        }
    }
    ####################
    
    if(verbose){
        message(sprintf("\tAmongst %d scores, there are %d positives.", length(scores), sum(scores>0)), appendLF=T)
        message(sprintf("Finally, find the subnetwork from the input network with %d nodes and %d edges...", vcount(g), ecount(g)), appendLF=T)
    }
    ## find the module with the maximum score
    module <- dNetFind(g, scores)
    if(verbose){
        message(sprintf("\tSize of the module: %d nodes and %d edges", vcount(module), ecount(module)), appendLF=T)
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    return(module)

}