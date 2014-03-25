#' Function to implement Random Walk with Restart (RWR) on the input graph. to pre-compute affinity matrix for the input graph
#'
#' \code{dRWR} is supposed to implement Random Walk with Restart (RWR) on the input graph. If the seeds (i.e. a set of starting nodes) are given, it intends to calculate the affinity score of all nodes in the graph to the given seeds. If the seeds are not give, it will pre-compute affinity matrix for nodes in the input graph with respect to each starting node (rather than a set of nodes) by looping over every node in the graph.
#'
#' @param g an object of class "igraph" or "graphNEL"
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param setSeeds an input matrix used to define sets of starting seeds. One column corresponds to one set of seeds that a walker starts with. The input matrix must have row names, coming from node names of input graph, i.e. V(g)$name, since there is a mapping operation. The non-zero entries mean that the corresonding rows (i.e. the gene/row names) are used as the seeds, and non-zero values can be viewed as how to weight the relative importance of seeds. By default, this option sets to "NULL" for the codebook matrix, suggesting each node in the graph will be used as a set of the seed to pre-compute affinity matrix for the input graph. This default does not scale for large input graphs since it will loop over every node in the graph; however, the pre-computed affinity matrix can be extensively reused for obtaining affinity scores between any combinations of nodes/seeds, allows for some flexibility in the downstream use, in particular when sampling a large number of random node combinations for statistical testing
#' @param restart the restart probability used for RWR. The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' 
#' @return 
#' When the seeds are NOT given, it returns:
#' \itemize{
#'  \item{\code{PTmatrix}: pre-computated affinity matrix with the dimension of n X n, where n is the number of nodes in the input graph. Columns stand for starting nodes walking from, and rows for ending nodes walking to. Therefore, a column for a starting node represents a steady-state affinity vector that the starting node will visit all the ending nodes in the graph}
#' }
#' When the seeds are given, it returns:
#' \itemize{
#'  \item{\code{PTmatrix}: affinity matrix with the dimension of n X nset, where n is the number of nodes in the input graph, and nset for the number of the sets of seeds (i.e. the number of columns in setSeeds). Each column stands for the steady probability vector, storing the affinity score of all nodes in the graph to the starting nodes/seeds. This steady probability vector can be viewed as the "influential impact" over the graph imposed by the starting nodes/seeds.}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{dNetInduce}}
#' @include dRWR.r
#' @examples
#' # 1) generate a random graph according to the ER model
#' g <- erdos.renyi.game(100, 1/100)
#'
#' # 2) produce the induced subgraph only based on the nodes in query
#' subg <- dNetInduce(g, V(g), knn=0)
#' V(subg)$name <- 1:vcount(subg)
#'
#' # 3) obtain the pre-computated affinity matrix
#' PTmatrix <- dRWR(subg, normalise="laplacian", restart=0.75)
#' # visualise affinity matrix
#' visHeatmapAdv(PTmatrix, Rowv=FALSE, Colv=FALSE, colormap="wyr", KeyValueName="Affinity")
#' 
#' # 4) obtain affinity matrix given sets of seeds
#' # define sets of seeds
#' # each seed with equal weight (i.e. all non-zero entries are '1')
#' aSeeds <- c(1,0,1,0,1)
#' bSeeds <- c(0,0,1,0,1)
#' setSeeds <- data.frame(aSeeds,bSeeds)
#' rownames(setSeeds) <- 1:5
#' # calcualte affinity matrix
#' PTmatrix <- dRWR(subg, normalise="laplacian", setSeeds=setSeeds, restart=0.75)
#' PTmatrix


dRWR <- function(g, normalise=c("laplacian","row","column","none"), setSeeds=NULL, restart=0.75)
{
    
    normalise <- match.arg(normalise)
    
    if(class(g)=="graphNEL"){
        ig <- igraph.from.graphNEL(g)
    }else{
        ig <- g
    }
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
    
    if ("weight" %in% list.vertex.attributes(ig)){
        adjM <- get.adjacency(ig, type="both", attr="weight", edges=F, names=T, sparse=F)
    }else{
        adjM <- get.adjacency(ig, type="both", attr=NULL, edges=F, names=T, sparse=F)
    }
    
    A <- adjM!=0
    if(normalise == "row"){
        D <- diag(apply(A,1,sum)^(-1))
        nadjM <- adjM %*% D
    }else if(normalise == "column"){
        D <- diag(apply(A,1,sum)^(-1))
        nadjM <- D %*% adjM
    }else if(normalise == "laplacian"){
        D <- diag(apply(A,1,sum)^(-0.5))
        nadjM <- D %*% adjM %*% D
    }else{
        nadjM <- adjM
    }
    
    #nig <- graph.adjacency(nadjM, mode=c("undirected"), weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)
    ## update the vertex attributes
    nattr <- list.vertex.attributes(ig)
    for(attr in nattr){
        #nig <- set.vertex.attribute(nig, name=attr, index=V(nig), get.vertex.attribute(ig,attr))
    }
    ## update the edge attributes    
    eattr <- list.edge.attributes(ig)
    for(attr in eattr){
        if (!("weight" %in% attr)){
            #nig <- set.edge.attribute(nig, name=attr, index=E(nig), get.edge.attribute(ig,attr))
        }
    }
    
    ################## RWR
    ## restarting prob
    if(is.null(restart) || is.na(restart) || restart<0 || restart>100){
        r <- 0.75
    }else if(restart>1 && restart<100){
        r <- restart/100
    }else{
        r <- restart
    }
    ## stopping critera
    stop_delta <- 1e-10   # L1 norm of successive estimates 'PT' below the threshold 'stop_delta'
    stop_step <- 100      # maximum steps of iterations
    
    if(is.null(setSeeds)){
        P0matrix <- matrix(as.numeric(nadjM>0),nrow=nrow(nadjM),ncol=ncol(nadjM))
        rownames(P0matrix) <- V(ig)$name
        colnames(P0matrix) <- V(ig)$name
    }else{
        ## check input data
        if(is.matrix(setSeeds) | is.data.frame(setSeeds)){
            data <- as.matrix(setSeeds)
        }else if(is.vector(setSeeds)){
            data <- as.matrix(setSeeds, ncol=1)
        }

        if(is.null(rownames(data))) {
            stop("The function must require the row names of the input setSeeds.\n")
        }else if(any(is.na(rownames(data)))){
            warning("setSeeds with NA as row names will be removed")
            data <- data[!is.na(rownames(data)),]
        }
        cnames <- colnames(data)
        if(is.null(cnames)){
            cnames <- seq(1,ncol(data))
        }
    
        ## check mapping between input data and graph
        ind <- match(rownames(data), V(ig)$name)
        nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
        if(length(nodes_mapped)!=vcount(ig)){
            warning("The row names of input setSeeds do not contain all those in the input graph.\n")
        }
        P0matrix <- matrix(0,nrow=nrow(nadjM),ncol=ncol(data))
        P0matrix[ind[!is.na(ind)],] <- as.matrix(data[!is.na(ind),])
        
        ## make sure the sum of elements in each starting probability vector is one
        P0matrix <- sapply(1:ncol(P0matrix), function(i){
            P0matrix[,i]/sum(P0matrix[,i])
        })
        rownames(P0matrix) <- V(ig)$name
        colnames(P0matrix) <- cnames
    }
    
    PTmatrix <- matrix(0, nrow=nrow(P0matrix), ncol=ncol(P0matrix))
    for(j in 1:ncol(P0matrix)){
        P0 <- as.matrix(P0matrix[,j],ncol=1)
        
        ## Initializing variables
        step <- 0
        delta <- 1
        
        PT <- P0
        ## Iterative update till convergence (delta<=1e-10)
        while (delta>stop_delta && step<=stop_step){
            PX <- (1-r) * nadjM %*% PT + r * P0
    
            # p-norm of v: sum((abs(v).p)^(1/p))
            delta <- sum(abs(PX-PT))
    
            PT <- PX
            step <- step+1
        }
        PTmatrix[,j] <- PT
    }
    
    ## make sure the sum of elements in each steady probability vector is one
    PTmatrix <- sapply(1:ncol(PTmatrix), function(i){
        PTmatrix[,i]/sum(PTmatrix[,i])
    })
    rownames(PTmatrix) <- rownames(P0matrix)
    colnames(PTmatrix) <- colnames(P0matrix)
        
    invisible(PTmatrix)
}