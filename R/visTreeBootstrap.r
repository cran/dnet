#' Function to build and visualise the bootstrapped tree
#'
#' \code{visTreeBootstrap} is supposed to build the tree, perform bootstrap analysis and visualise the bootstrapped tree. It returns an object of class "phylo". For easy downstream analysis, the bootstrapped tree is rerooted either at the internal node with the miminum bootstrap/confidence value or at any customised internal node. 
#'
#' @param data an input data matrix used to build the tree. The built tree describes the relationships between rows of input matrix
#' @param algorithm the tree-building algorithm. It can be one of "nj" for the neighbor-joining tree estimation, "fastme.ols" for the minimum evolution algorithm with ordinary least-squares (OLS) fitting of a metric to a tree structure, and "fastme.bal" for the minimum evolution algorithm under a balanced (BAL) weighting scheme
#' @param metric distance metric used to calculate a distance matrix between rows of input matrix. It can be: "pearson" for pearson correlation, "spearman" for spearman rank correlation, "kendall" for kendall tau rank correlation, "euclidean" for euclidean distance, "manhattan" for cityblock distance, "cos" for cosine similarity, "mi" for mutual information
#' @param num.bootstrap an integer specifying the number of bootstrap replicates
#' @param consensus logical to indicate whether to return the consensus tree. By default, it sets to false for not doing so. Note: if true, there will be no visualisation of the bootstrapped tree
#' @param consensus.majority a numeric value between 0.5 and 1 (or between 50 and 100) giving the proportion for a clade to be represented in the consensus tree
#' @param reroot determines if and how the bootstrapped tree should be rerooted. By default, it is "min.bootstrap", which implies that the bootstrapped tree will be rerooted at the internal node with the miminum bootstrap/confidence value. If it is an integer between 1 and the number of internal nodes, the tree will be rerooted at the internal node with this index value
#' @param plot.phylo.arg a list of main parameters used in the function "ape::plot.phylo" \url{http://www.inside-r.org/packages/cran/ape/docs/plot.phylo}. See 'Note' below for details on the parameters
#' @param nodelabels.arg a list of main parameters used in the function "ape::nodelabels" \url{http://www.inside-r.org/packages/cran/ape/docs/nodelabels}. See 'Note' below for details on the parameters
#' @param visTree logical to indicate whether the bootstrap tree will be visualised. By default, it sets to true for display. Note, the consensus tree can not be enabled for visualisation
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param ... additional "ape::plot.phylo" parameters
#' @return 
#' an object of class "phylo". It can return a bootstrapped tree or a consensus tree (if enabled):
#' When a bootstrapped tree is returned (also visualised by default), the "phylo" object has a list with following components:
#' \itemize{
#'  \item{\code{Nnode}: the number of internal nodes}
#'  \item{\code{node.label}: the labels for internal nodes. Here, each internal node is associated with the bootstrap value}
#'  \item{\code{tip.label}: the lables for tip nodes. Tip labels come from the row names of the input matrix, but are not necessarily the same order as they appear in the input matrix}
#'  \item{\code{edge}: a two-column matrix describing the links between tree nodes (including internal and tip nodes)}
#'  \item{\code{edge.length}: a vector indicating the edge length in the 'edge'}
#'  \item{Note: the tree structure is indexed with 1:Ntip for tip nodes, and (\eqn{Ntip}+1):(\eqn{Ntip}+\eqn{Nnode}) for internal nodes, where \eqn{Ntip} is the number of tip nodes and \eqn{Nnode} for the number of internal nodes. Moreover, \eqn{nrow(data)=Ntip=Nnode-2}.}
#' }
#' When a consensus tree is returned (no visualisation), the "phylo" object has a list with following components:
#' \itemize{
#'  \item{\code{Nnode}: the number of internal nodes}
#'  \item{\code{tip.label}: the lables for tip nodes. Tip labels come from the row names of the input matrix, but are not necessarily the same order as they appear in the input matrix}
#'  \item{\code{edge}: a two-column matrix describing the links between tree nodes (including internal and tip nodes)}
#' }
#' @note 
#' A list of main parameters used in the function "ape::plot.phylo":
#' \itemize{
#' \item{"type": a character string specifying the type of phylogeny to be drawn; it must be one of "phylogram" (the default), "cladogram", "fan", "unrooted", "radial" or any unambiguous abbreviation of these}
#' \item{"direction": a character string specifying the direction of the tree. Four values are possible: "rightwards" (the default), "leftwards", "upwards", and "downwards"}
#' \item{"lab4ut": (= labels for unrooted trees) a character string specifying the display of tip labels for unrooted trees: either "horizontal" where all labels are horizontal (the default), or "axial" where the labels are displayed in the axis of the corresponding terminal branches. This option has an effect only if type = "unrooted"}
#' \item{"edge.color": a vector of mode character giving the colours used to draw the branches of the plotted phylogeny. These are taken to be in the same order than the component edge of phy. If fewer colours are given than the length of edge, then the colours are recycled}
#' \item{"edge.width": a numeric vector giving the width of the branches of the plotted phylogeny. These are taken to be in the same order than the component edge of phy. If fewer widths are given than the length of edge, then these are recycled}
#' \item{"edge.lty": same than the previous argument but for line types; 1: plain, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, 6: twodash}
#' \item{"font": an integer specifying the type of font for the labels: 1 (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic)}
#' \item{"cex": a numeric value giving the factor scaling of the tip and node labels (Character EXpansion). The default is to take the current value from the graphical parameters}
#' \item{"adj": a numeric specifying the justification of the text strings of the labels: 0 (left-justification), 0.5 (centering), or 1 (right-justification). This option has no effect if type="unrooted". If NULL (the default) the value is set with respect of direction (see details)}
#' \item{"srt": a numeric giving how much the labels are rotated in degrees (negative values are allowed resulting in clock-like rotation); the value has an effect respectively to the value of direction (see Examples). This option has no effect if type="unrooted"}
#' \item{"no.margin": a logical. If TRUE, the margins are set to zero and the plot uses all the space of the device}
#' \item{"label.offset": a numeric giving the space between the nodes and the tips of the phylogeny and their corresponding labels. This option has no effect if type="unrooted"}
#' \item{"rotate.tree": for "fan", "unrooted", or "radial" trees: the rotation of the whole tree in degrees (negative values are accepted}
#' }
#' A list of main parameters used in the function "ape::nodelabels":
#' \itemize{
#' \item{"text": a vector of mode character giving the text to be printed. By default, the labels for internal nodes (see "node.label"), that is, the bootstrap values associated with internal nodes}
#' \item{"node": a vector of mode numeric giving the numbers of the nodes where the text or the symbols are to be printed. By default, indexes for internal nodes, that is, (\eqn{Ntip}+1):(\eqn{Ntip}+\eqn{Nnode}), where \eqn{Ntip} is the number of tip nodes and \eqn{Nnode} for the number of internal nodes}
#' \item{"adj": one or two numeric values specifying the horizontal and vertical, respectively, justification of the text or symbols. By default, the text is centered horizontally and vertically. If a single value is given, this alters only the horizontal position of the text}
#' \item{"frame": a character string specifying the kind of frame to be printed around the text. This must be one of "rect" (the default), "circle", "none", or any unambiguous abbreviation of these}
#' \item{"cex": a numeric value giving the factor scaling of the tip and node labels (Character EXpansion). The default is to take the current value from the graphical parameters}
#' \item{"font": an integer specifying the type of font for the labels: 1 (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic)}
#' \item{"col": a character string giving the color to be used for the text or the plotting symbols; this is eventually recycled}
#' \item{"bg": a character string giving the color to be used for the background of the text frames or of the plotting symbols if it applies; this is eventually recycled. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}}
#' }
#' @export
#' @importFrom ape nj fastme.ols fastme.bal boot.phylo consensus is.rooted unroot root mrca write.tree read.tree plot.phylo nodelabels
#' @seealso \code{\link{visTreeBootstrap}}
#' @include visTreeBootstrap.r
#' @examples
#' # 1) generate an iid normal random matrix of 100x10 
#' data <- matrix( rnorm(100*10,mean=0,sd=1), nrow=100, ncol=10)
#' colnames(data) <- paste(rep('S',10), seq(1:10), sep="")
#' data <- t(data)
#'
#' # 2) build neighbor-joining tree with bootstrap values and visualise it by default
#' visTreeBootstrap(data)
#'
#' # 3) only display those internal nodes with bootstrap values > 30
#' # 3a) generate the bootstrapped tree (without visualisation)
#' tree_bs <- visTreeBootstrap(data, visTree=FALSE)
#' # 3b) look at the bootstrap values and ordered row names of input matrix
#' # the bootstrap values
#' tree_bs$node.label
#' # ordered row names of input matrix
#' tree_bs$tip.label
#' # 3c) determine internal nodes that should be displayed
#' Ntip <- length(tree_bs$tip.label) # number of tip nodes
#' Nnode <- length(tree_bs$node.label) # number of internal nodes
#' flag <- as.numeric(tree_bs$node.label) > 30
#' text <- tree_bs$node.label[flag]
#' node <- Ntip + (1:Nnode)[flag]
#' visTreeBootstrap(data, nodelabels.arg=list(text=text,node=node))
#'
#' # 4) obtain the consensus tree
#' tree_cons <- visTreeBootstrap(data, consensus=TRUE, num.bootstrap=10)

visTreeBootstrap <- function(data, algorithm=c("nj","fastme.ols","fastme.bal"), metric=c("euclidean","pearson","spearman","cos","manhattan","kendall","mi"), num.bootstrap=100, consensus=FALSE, consensus.majority=0.5, reroot="min.bootstrap", plot.phylo.arg=NULL, nodelabels.arg=NULL, visTree=T, verbose=T, ...)
{


    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    algorithm <- match.arg(algorithm)
    metric <- match.arg(metric)

    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
    }else if(is.null(data)){
        stop("The input data must be matrix.\n")
    }

    if(is.null(rownames(data))) {
        warning("There are no row names of the input data, which are assigned by a series of integers.")
        rownames(data) <- 1:nrow(data)
    }
    tip.label <- rownames(data)
    
    ## temporally make sure the row names are unique
    rownames(data) <- paste(rownames(data), 1:nrow(data), sep=".")
    
    if(verbose){
        message(sprintf("First, build the tree (using %s algorithm and %s distance) from input matrix (%d by %d)...", algorithm, metric, nrow(data), ncol(data)), appendLF=T)
    }
    ## build the tree
    d <- as.dist(sDistance(data, metric=metric))
    tr <- do.call(algorithm, list(d))
    
    if(verbose){
        message(sprintf("Second, perform bootstrap analysis with %d replicates...", num.bootstrap), appendLF=T)
    }
    ## perform bootstrap analysis
    f <- function(x) {
        d <- as.dist(sDistance(x, metric=metric))
        do.call(algorithm, list(d))
    }
    if(consensus==TRUE){
        res <- ape::boot.phylo(tr, data, f, B=num.bootstrap, block=1, quiet=TRUE, trees=TRUE)
        bp <- res$BP
        bp <- ceiling(100*bp/num.bootstrap)
        
        ## calculate the consensus tree
        if(consensus.majority>=50 & consensus.majority<=100){
            consensus.majority <- consensus.majority/100
        }
        if(consensus.majority<0.5 | consensus.majority>1){
            consensus.majority <- 0.5
        }
        
        if(verbose){
            message(sprintf("Finally, obtain consensus tree based on %1.2f majority...", consensus.majority), appendLF=T)
        }
        tree_cons <- ape::consensus(res$trees, p=consensus.majority, check.labels=T)
        tree_cons$tip.label <- sub("\\.\\d+$", "", tree_cons$tip.label)
        
        #return(tree_cons)
    }else{
        bp <- ape::boot.phylo(tr, data, f, B=num.bootstrap, block=1, quiet=TRUE, trees=FALSE)
        bp <- ceiling(100*bp/num.bootstrap)
    
        ## restore the original names
        tr$tip.label <- sub("\\.\\d+$", "", tr$tip.label)
        ## assign the bootstrap values onto node.label
        tr$node.label <- bp
    
        if(verbose){
            message(sprintf("Finally, visualise the bootstrapped tree..."), appendLF=T)
        }
     
        ########################################################################
        ## unroot the tree
        if(ape::is.rooted(tr)) tr <- ape::unroot(tr)
    
        ## reroot the tree either according to miminum bootstrap value or the bp vector index (if given)
        tmp_bs <- as.numeric(tr$node.label)
        if(is.integer(reroot) & reroot>=1 & reroot<=length(tr$node.label)){
            ## miminum bootstrap value and the node index
            reroot_index_mrca <- reroot+length(tr$tip.label)
            tree_bs <- ape::root(tr, node=reroot_index_mrca, resolve.root=F, interactive=F)
        }else if(reroot=="min.bootstrap"){
            ## miminum bootstrap value and the node index
            min_bs <- min(tmp_bs[!is.na(tmp_bs)])
            min_bs_index <- which(tmp_bs==min_bs)[1]
            reroot_index_mrca <- min_bs_index+length(tr$tip.label)
            tree_bs <- ape::root(tr, node=reroot_index_mrca, resolve.root=F, interactive=F)
        }else{
            tree_bs <- tr
        }
    
        ## write and re-read the reroot tree
        ## this ensures the tree structure indexing is "as-is"
        newick_tmp <- "MyNewickTreefile.tre"
        ape::write.tree(tree_bs, file="MyNewickTreefile.tre")
        tree_bs <- ape::read.tree(file=newick_tmp)
        unlink(newick_tmp) # delete the file
    }
    
    if(consensus==FALSE & visTree==TRUE){
        ########################################################################
        ## Define the tree struct dimensions and indexing
        Ntip <- length(tree_bs$tip.label)
        Nnode <- tree_bs$Nnode
        tip_index <- c(1:Ntip)
        node_index <- c((Ntip+1):(Ntip+Nnode))
        
        if(0){
        ## most recent common ancestor (MRCA) for each pair of tips and nodes
        mrca_node <- ape::mrca(tree_bs, full=T)
        
        visHeatmapAdv(mrca_node, Rowv=F,Colv=F, zlim=c(Ntip+1,Ntip+Nnode), colormap="gray-black", add.expr=abline(v=c(1,Ntip+1,(ncol(mrca_node)+1))-0.5, h=c(1,Nnode+1,(ncol(mrca_node)+1))-0.5, col="white"), KeyValueName="MRCA index", lmat=rbind(c(4,3), c(2,1)), lhei=c(1,5), lwid=c(1,3))
        
        ## connectivity linking each ancestor to its all children
        ## matrix of Nnode by (Ntip+Nnode)
        connectivity <- array(0, c(Nnode,Ntip+Nnode))
        for (i in 1:Nnode) {
            node_tmp <- i+Ntip
            child <- which(mrca_node[node_tmp,]==node_tmp, arr.ind=T)
            connectivity[i,child] <- 1
            # exclude self
            connectivity[i,node_tmp] <- 0
        }
        rownames(connectivity) <- node_index
        visHeatmapAdv(connectivity, Rowv=F,Colv=F, zlim=c(0,1), colormap="gray-black", add.expr=abline(v=c(1,Ntip+1,(ncol(mrca_node)+1))-0.5, col="white"), key=F, cexRow=1,cexCol=1)
        }
        
        ########################################################################

        ## a function to update the parameters
        update_parameters <- function(default, update){
            for(i in 1:length(names(default))){
                item <- names(default)[i]
                if(item %in% names(update)){
                    default[i] <- update[which(names(update)==item)]
                }
            }
            return(default)
        }
    
        ################
        #plot.phylo.arg <- list(direction="upwards",lab4ut="axial")
    
        ## the default parameters for "plot.phylo"
        plot.phylo.arg.default <- list(
                    type = "phylogram",
                    direction = "rightwards",
                    lab4ut = "horizontal",
                    edge.color = "black",
                    edge.width = 1,
                    edge.lty = 1,
                    font = 3,
                    cex = 1,
                    adj = 0,
                    srt = 0,
                    no.margin = T,
                    label.offset = 0,
                    rotate.tree = 0
                    )
    
        ## update parameters for "plot.phylo"
        plot.phylo.arg.default <- update_parameters(plot.phylo.arg.default, plot.phylo.arg)
    
        ################
        #nodelabels.arg <- list(frame="rect")
    
        ## the default parameters for "nodelabels"
        nodelabels.arg.default <- list(
                    text = tree_bs$node.label,
                    node = Ntip+ (1:Nnode),
                    adj = c(0.5, 0.5),
                    frame = "circle",
                    cex = 1, 
                    font = 4,
                    col = "black",
                    bg = "white-lightyellow-darkorange"
                    )
    
        ## update parameters for "nodelabels"
        nodelabels.arg.default <- update_parameters(nodelabels.arg.default, nodelabels.arg)
    
        ##################
        ape::plot.phylo(tree_bs, 
                        type = plot.phylo.arg.default$type,
                        direction = plot.phylo.arg.default$direction,
                        lab4ut = plot.phylo.arg.default$lab4ut,
                        edge.color = plot.phylo.arg.default$edge.color,
                        edge.width = plot.phylo.arg.default$edge.width,
                        edge.lty = plot.phylo.arg.default$edge.lty,
                        font = plot.phylo.arg.default$font,
                        cex = plot.phylo.arg.default$cex,
                        adj = plot.phylo.arg.default$adj,
                        srt = plot.phylo.arg.default$srt,
                        no.margin = plot.phylo.arg.default$no.margin,
                        label.offset = plot.phylo.arg.default$label.offset,
                        rotate.tree = plot.phylo.arg.default$rotate.tree,
                        show.tip.label = TRUE,
                        show.node.label = FALSE,
                        ...
                        )
        ape::nodelabels(
                    text = nodelabels.arg.default$text,
                    node = nodelabels.arg.default$node,
                    adj = nodelabels.arg.default$adj,
                    frame = nodelabels.arg.default$frame,
                    cex = nodelabels.arg.default$cex, 
                    font = nodelabels.arg.default$font,
                    col = nodelabels.arg.default$col,
                    bg = visColormap(colormap=nodelabels.arg.default$bg)(101)[1 + as.numeric(nodelabels.arg.default$text)]
                   )
   
        ########################################################################
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    if(consensus==TRUE){
        return(tree_cons)
    }else{
        return(tree_bs)
    }
    
}