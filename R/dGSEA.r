#' Function to conduct gene set enrichment analysis given the input data and the ontology in query
#'
#' \code{dGSEA} is supposed to conduct gene set enrichment analysis given the input data and the ontology in query. It returns an object of class "eTerm". 
#'
#' @param data a data frame or matrix of input data. It must have row names, either Entrez Gene ID or Symbol
#' @param identity the type of gene identity (i.e. row names of input data), either "symbol" for gene symbols (by default) or "entrez" for Entrez Gene ID. The option "symbol" is preferred as it is relatively stable from one update to another; when gene symbols cannot be matched, synonyms will be searched against
#' @param genome the genome identity. It can be either "mm" for mouse genome or "hs" for human genome
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "MP" for Mammalian Phenotype, "HP" for Human Phenotype, "DO" for Disease Ontology, and "PS" for phylostratific age 
#' @param sizeRange the minimum and maximum size of members of each gene set in consideration. By default, it sets to a minimum of 10 but no more than 1000
#' @param which_distance which distance of terms in the ontology is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param weight type of score weigth. It can be "0" for unweighted (an equivalent to Kolmogorov-Smirnov, only considering the rank), "1" for weighted by input gene score (by default), and "2" for over-weighted, and so on
#' @param nperm the number of random permutations. For each permutation, gene-score associations will be permutated so that permutation of gene-term associations is realised
#' @param fast logical to indicate whether to fast calculate expected results from permutated data. By default, it sets to true
#' @param sigTail the tail used to calculate the statistical significance. It can be either "two-tails" for the significance based on two-tails  or "one-tail" for the significance based on one tail
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @return 
#' an object of class "eTerm", a list with following components:
#' \itemize{
#'  \item{\code{set_info}: a matrix of nSet X 4 containing gene set information, where nSet is the number of gene set in consideration, and the 4 columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
#'  \item{\code{gs}: a list of gene sets, each storing gene members. Always, gene sets are identified by "setID" and gene members identified by "Entrez ID"}
#'  \item{\code{data}: a matrix of nGene X nSample containing input data in consideration. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{es}: a matrix of nSet X nSample containing enrichment score, where nSample is the number of samples (i.e. the number of columns in input data}
#'  \item{\code{nes}: a matrix of nSet X nSample containing normalised enrichment score. It is the version of enrichment score but after being normalised by gene set size}
#'  \item{\code{pvalue}: a matrix of nSet X nSample containing nominal p value}
#'  \item{\code{adjp}: a matrix of nSet X nSample containing adjusted p value. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{gadjp}: a matrix of nSet X nSample containing globally adjusted p value in terms of all samples}
#'  \item{\code{fdr}: a matrix of nSet X nSample containing false discovery rate (FDR). It is the estimated probability that the normalised enrichment score represents a false positive finding}
#'  \item{\code{qvalue}: a matrix of nSet X nSample containing q value. It is the monotunically increasing FDR}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The interpretation of returned components:
#' \itemize{
#' \item{"es": enrichment score for the gene set is the degree to which this gene set is overrepresented at the top or bottom of the ranked list of genes in each column of input data;}
#' \item{"nes": normalised enrichment score for the gene set is enrichment score that has already normalised by gene set size. It is comparable across analysed gene sets;}
#' \item{"pvalue": nominal p value is the statistical significance of the enrichment score. It is not adjusted for multiple hypothesis testing, and thus is of limited use in comparing gene sets;}
#' \item{"adjp": adjusted p value by Benjamini & Hochberg method. It is comparable across gene sets;}
#' \item{"gadjp": globally adjusted p value by Benjamini & Hochberg method. Unlike "adjp", it is adjusted in terms of all samples;}
#' \item{"fdr": false discovery rate is the estimated probability that the normalised enrichment score represents a false positive finding. Unlike "adjp" or "gadjp" (also aliased as "fdr") that is derived from a list of p values, this version of fdr is directly calculate from the statistic (i.e. normalised enrichment score);}
#' \item{"qvalue": q value is the monotunically increasing FDR so that the higher "nes", the lower "qvalue".}
#' }
#' @export
#' @seealso \code{\link{dGSEAview}}, \code{\link{dGSEAwrite}}, \code{\link{visGSEA}}
#' @include dGSEA.r
#' @examples
#' \dontrun{
#' load(url("http://dnet.r-forge.r-project.org/data/Hiratani_TableS1.RData"))
#' data <- RT[,1:2]
#' eTerm <- dGSEA(data, identity="symbol", genome="mm", ontology="MP", which_distance=c(1,2))
#' }

dGSEA <- function(data, identity=c("symbol","entrez"), genome=c("mm", "hs"), ontology=c("GOBP","GOMF","GOCC","MP","DO","PS"), sizeRange=c(10,1000), which_distance=NULL, weight=1, nperm=100, fast=T, sigTail=c("two-tails","one-tail"), verbose=T)
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    identity <- match.arg(identity)
    genome <- match.arg(genome)
    ontology <- match.arg(ontology)
    sigTail <- match.arg(sigTail)
    
    if (is.vector(data)){
        data <- matrix(data, nrow=length(data), ncol=1)
    }else if(is.matrix(data) | is.data.frame(data)){
        data <- as.matrix(data)
    }else if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }
    
    if(is.null(rownames(data))){
        stop("The input data must have row names with attached gene id or symbols.\n")
    }
    if(is.null(colnames(data))){
        colnames(data) <- seq(1,ncol(data))
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("First, load the ontology %s and its gene associations in the genome %s (%s) ...", ontology, genome, as.character(now)), appendLF=T)
    }
    
    EG <- list()
    GS <- list()
    if(genome == "mm"){
        
        #data("org.Mm.eg", package="dnet")
        load(url("http://dnet.r-forge.r-project.org/data/org.Mm.eg.RData"))
        eval(parse(text="EG <- org.Mm.eg"))
        #do.call(assign, list("EG", org.Mm.eg))
        
        if(ontology == "GOBP"){
            #data("org.Mm.egGOBP", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Mm.egGOBP.RData"))
            eval(parse(text="GS <- org.Mm.egGOBP"))
        }else if(ontology == "GOMF"){
            #data("org.Mm.egGOMF", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Mm.egGOMF.RData"))
            eval(parse(text="GS <- org.Mm.egGOMF"))
        }else if(ontology == "GOCC"){
            #data("org.Mm.egGOCC", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Mm.egGOCC.RData"))
            eval(parse(text="GS <- org.Mm.egGOCC"))
        }else if(ontology == "MP"){
            #data("org.Mm.egMP", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Mm.egMP.RData"))
            eval(parse(text="GS <- org.Mm.egMP"))
        }else if(ontology == "HP"){
            #data("org.Mm.egHP", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Mm.egHP.RData"))
            eval(parse(text="GS <- org.Mm.egHP"))
        }else if(ontology == "DO"){
            #data("org.Mm.egDO", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Mm.egDO.RData"))
            eval(parse(text="GS <- org.Mm.egDO"))
        }else if(ontology == "PS"){
            #data("org.Mm.egPS", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Mm.egPS.RData"))
            eval(parse(text="GS <- org.Mm.egPS"))
        }
    }else if(genome == "hs"){
        
        #data("org.Hs.eg", package="dnet")
        load(url("http://dnet.r-forge.r-project.org/data/org.Hs.eg.RData"))
        eval(parse(text="EG <- org.Hs.eg"))
        #do.call(assign, list("EG", org.Hs.eg))
        
        if(ontology == "GOBP"){
            #data("org.Hs.egGOBP", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Hs.egGOBP.RData"))
            eval(parse(text="GS <- org.Hs.egGOBP"))
        }else if(ontology == "GOMF"){
            #data("org.Hs.egGOMF", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Hs.egGOMF.RData"))
            eval(parse(text="GS <- org.Hs.egGOMF"))
        }else if(ontology == "GOCC"){
            #data("org.Hs.egGOCC", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Hs.egGOCC.RData"))
            eval(parse(text="GS <- org.Hs.egGOCC"))
        }else if(ontology == "MP"){
            #data("org.Hs.egMP", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Hs.egMP.RData"))
            eval(parse(text="GS <- org.Hs.egMP"))
        }else if(ontology == "HP"){
            #data("org.Hs.egMP", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Hs.egHP.RData"))
            eval(parse(text="GS <- org.Hs.egHP"))
        }else if(ontology == "DO"){
            #data("org.Hs.egDO", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Hs.egDO.RData"))
            eval(parse(text="GS <- org.Hs.egDO"))
        }else if(ontology == "PS"){
            #data("org.Hs.egPS", package="dnet")
            load(url("http://dnet.r-forge.r-project.org/data/org.Hs.egPS.RData"))
            eval(parse(text="GS <- org.Hs.egPS"))
        }
    }
    
    allGeneID <- EG$gene_info$GeneID
    allSymbol <- as.vector(EG$gene_info$Symbol)
    allSynonyms <- as.vector(EG$gene_info$Synonyms)
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Then, do mapping based on %s (%s) ...", identity, as.character(now)), appendLF=T)
    }
    
    if(identity == "symbol"){
    
        Symbol <- rownames(data)
        
        if(1){
            ## for those starting with 'Mar' in a excel-input date format
            a <- Symbol
            flag <- grep("-Mar$", a, ignore.case=T, perl=T, value=F)
            if(length(flag)>=1){
                b <- a[flag]
                c <- sub("-Mar$", "", b, ignore.case=T, perl=T)
                d <- sub("^0", "", c, ignore.case=T, perl=T)
                e <- sapply(d, function(x) paste(c("March",x), collapse=""))
                a[flag] <- e
                Symbol <- a
            }

            ## for those starting with 'Sep' in a excel-input date format
            a <- Symbol
            flag <- grep("-Sep$", a, ignore.case=T, perl=T, value=F)
            if(length(flag)>=1){
                b <- a[flag]
                c <- sub("-Sep$", "", b, ignore.case=T, perl=T)
                d <- sub("^0", "", c, ignore.case=T, perl=T)
                e <- sapply(d, function(x) paste(c("Sept",x), collapse=""))
                a[flag] <- e
                Symbol <- a
            }
        }
        
        ## case-insensitive
        match_flag <- match(tolower(Symbol),tolower(allSymbol))
        
        ## match Synonyms (if not found via Symbol)
        na_flag <- is.na(match_flag)
        a <- Symbol[na_flag]

        ###
        tmp_flag <- is.na(match(tolower(allSymbol), tolower(Symbol)))
        tmp_Synonyms <- allSynonyms[tmp_flag]
        Orig.index <- seq(1,length(allSynonyms))
        Orig.index <- Orig.index[tmp_flag]
        ###

        b <- sapply(1:length(a), function(x){
            tmp_pattern1 <- paste("^",a[x],"\\|", sep="")
            tmp_pattern2 <- paste("\\|",a[x],"\\|", sep="")
            tmp_pattern3 <- paste("\\|",a[x],"$", sep="")
            tmp_pattern <- paste(tmp_pattern1,"|",tmp_pattern2,"|",tmp_pattern3, sep="")
            tmp_result <- grep(tmp_pattern, tmp_Synonyms, ignore.case=T, perl=T, value=F)
            ifelse(length(tmp_result)==1, Orig.index[tmp_result[1]], NA)
        })
        match_flag[na_flag] <- b
        GeneID <- allGeneID[match_flag]
        
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tAmong %d symbols of input data, there are %d mappable via official gene symbols, %d mappable via alias and %d unmappable", length(Symbol), (length(Symbol)-length(a)), sum(!is.na(b)), sum(is.na(b))), appendLF=T)
        }
        
    }else{
        GeneID <- rownames(data)
        match_flag <- match(GeneID,allGeneID)
        GeneID <- allGeneID[match_flag]
    }
       
    flag <- !is.na(GeneID)
    data <- as.matrix(data[flag,])
    GeneID <- GeneID[flag]
    
    ## Take average for those identical rows
    GeneID_uni <- sort(unique(GeneID))
    tmp_matrix <- matrix(0, nrow=length(GeneID_uni),ncol=ncol(data))
    for(i in 1:length(GeneID_uni)){
        flag <- which(GeneID_uni[i]==GeneID)
        if(length(flag)==1){
            tmp_matrix[i,] <- data[flag,]
        }else{
            tmp_matrix[i,] <- colMeans(data[flag,])
        }
    }
    colnames(tmp_matrix) <- colnames(data)
    rownames(tmp_matrix) <- GeneID_uni
    data <- tmp_matrix

    nGene <- nrow(data)
    nSample <- ncol(data)
    geneid <- rownames(data) ## only those genes in question are considered
    
    
    ## filter based on "which_distance"
    if(!is.null(which_distance) & sum(is.na(GS$set_info$distance))==0){
        set_filtered <- sapply(which_distance, function(x) {
            GS$set_info$setID[(GS$set_info$distance==as.integer(x))]
        })
        set_filtered <- unlist(set_filtered)
    }else{
        set_filtered <- GS$set_info$setID
    }
    
    ## derive the "gs" of interest
    gs <- list()
    set_info <- data.frame()
    gs_names <- vector()
    k <- 1
    for(i in 1:length(GS$gs)){
        ## make sure only those meeting the distance criteria
        if(!is.na(match(names(GS$gs[i]),set_filtered))){
            gs_origin <- GS$gs[[i]]
            flag <- match(gs_origin, geneid)
            gs_filtered <- geneid[flag[!is.na(flag)]]
            if(length(gs_filtered)>=sizeRange[1] & length(gs_filtered)<=sizeRange[2]){
                gs[k] <- list(gs_filtered)
                gs_names[k] <- names(GS$gs[i])
                
                flag <- match(names(GS$gs[i]),GS$set_info$setID)
                set_info <- rbind(set_info,GS$set_info[flag,])
                
                k <- k+1
            }
        }
    }
    names(gs) <- gs_names
    nSet <- length(gs)
    
    ## Enrichment score for the gene set; that is, the degree to which this gene set is overrepresented at the top or bottom of the ranked list of genes in the expression dataset
    SS.es <- matrix(0, nrow=nSet, ncol=nSample)
    colnames(SS.es) <- colnames(data)
    rownames(SS.es) <- names(gs)
    ## Normalized enrichment score; that is, the enrichment score for the gene set after it has been normalized across analyzed gene sets
    SS.nes <- matrix(0, nrow=nSet, ncol=nSample)
    colnames(SS.nes) <- colnames(data)
    rownames(SS.nes) <- names(gs)
    ## Nominal p value; that is, the statistical significance of the enrichment score. The nominal p value is not adjusted for gene set size or multiple hypothesis testing; therefore, it is of limited use in comparing gene sets
    SS.pvalue <- matrix(1, nrow=nSet, ncol=nSample)
    colnames(SS.pvalue) <- colnames(data)
    rownames(SS.pvalue) <- names(gs)
    ## Familywise-error rate; that is, a more conservatively estimated probability that the normalized enrichment score represents a false positive finding. For the goal of GSEA is to generate hypotheses, recommend focusing on the FDR statistic
    SS.fwer <- matrix(1, nrow=nSet, ncol=nSample)
    colnames(SS.fwer) <- colnames(data)
    rownames(SS.fwer) <- names(gs)
    ## False discovery rate; that is, the estimated probability that the normalized enrichment score represents a false positive finding
    SS.fdr <- matrix(1, nrow=nSet, ncol=nSample)
    colnames(SS.fdr) <- colnames(data)
    rownames(SS.fdr) <- names(gs)
    ## FDR q-values; that is, monotunically increasing FDR 
    SS.qvalue <- matrix(1, nrow=nSet, ncol=nSample)
    colnames(SS.qvalue) <- colnames(data)
    rownames(SS.qvalue) <- names(gs)
    
    ## adjusted p value: Adjust P-values for Multiple Comparisons
    SS.adjp <- matrix(1, nrow=nSet, ncol=nSample)
    colnames(SS.adjp) <- colnames(data)
    rownames(SS.adjp) <- names(gs)
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Third, perform GSEA analysis (%s) ...", as.character(now)), appendLF=T)
    }
    
    for(j in 1:nSample){
        
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tSample %d is being processed at (%s) ...", j, as.character(now)), appendLF=T)
        }
        
        rank.score <- data[,j]
        
        ind <- order(rank.score, decreasing=T)
        rank.score.sorted <- rank.score[ind]
        geneid.sorted <- geneid[ind]
        geneid2ind <- ind
        
        ## shuffle the gene labels
        geneid.shuffled <- matrix(0, nrow=length(geneid.sorted), ncol=nperm)
        for(r in 1:nperm){
            geneid.shuffled[,r] <- sample(geneid2ind)
        }
        
        es.position <- matrix(0, nrow=nSet, ncol=1)
        es.observed <- matrix(0, nrow=nSet, ncol=1)
        es.expected <- matrix(0, nrow=nSet, ncol=nperm)
        for(k in 1:nSet){
            
            if(verbose){
                if(k %% 100==0 | k==nSet){
                    now <- Sys.time()
                    message(sprintf("\t\t %d of %d gene sets have been processed (%s) ...", k, nSet, as.character(now)), appendLF=T)
                }
            }
            
            nHit <- length(gs[[k]])
            nMiss <- nGene - nHit
            
            ########################################
            ## observed
            observed.point <- rep(-1/nMiss, nGene)
            flag <- match(gs[[k]], geneid.sorted)
            flag_obs <- flag
            if(weight==0) {
                observed.point[flag] <- 1/nHit
            }else if(weight==1){
                hit_tmp <- abs(rank.score.sorted[flag])
                observed.point[flag] <- hit_tmp/sum(hit_tmp)
            }else{
                hit_tmp <- abs(rank.score.sorted[flag] ** weight)
                observed.point[flag] <- hit_tmp/sum(hit_tmp)      
            }
            observed.cumsum <- cumsum(observed.point)
            max.ES <- max(observed.cumsum)
            min.ES <- min(observed.cumsum)
            es.observed[k] <- signif(ifelse(max.ES>abs(min.ES), max.ES, min.ES), digits=5)
            es.position[k] <- ifelse(max.ES>abs(min.ES), which.max(observed.cumsum), which.min(observed.cumsum))
            
            ########################################
            ## expected
            if(fast == F){
                
                for(r in 1:nperm){
                    expected.point <- rep(-1/nMiss, nGene)
                    flag <- geneid.shuffled[flag_obs,r]
                    if(weight==0) {
                        expected.point[flag] <- 1/nHit
                    }else if(weight==1){
                        hit_tmp <- abs(rank.score.sorted[flag])
                        expected.point[flag] <- hit_tmp/sum(hit_tmp)
                    }else{
                        hit_tmp <- abs(rank.score.sorted[flag] ** weight)
                        expected.point[flag] <- hit_tmp/sum(hit_tmp)      
                    }
                    expected.cumsum <- cumsum(expected.point)
                    
                    max.ES <- max(expected.cumsum)
                    min.ES <- min(expected.cumsum)
                    es.expected[k,r] <- signif(ifelse(max.ES>abs(min.ES), max.ES, min.ES), digits=5)
                }  
                
            }else{
                
                ## fast calculation
                ## faster (x5) by cumsum only those in hits
                ## supposed to assess the enrichment of random permutations
                
                N <- nGene
                Nh <- nHit
                Nm <- nMiss
                
                for(r in 1:nperm){
                    flag <- geneid.shuffled[flag_obs,r]
            
                    loc.vector <- vector(length=N, mode="numeric")
                    peak.res.vector <- vector(length=Nh, mode="numeric")
                    valley.res.vector <- vector(length=Nh, mode="numeric")
                    tag.correl.vector <- vector(length=Nh, mode="numeric")
                    tag.diff.vector <- vector(length=Nh, mode="numeric")
                    tag.loc.vector <- vector(length=Nh, mode="numeric")
            
                    loc.vector <- seq(1, N)
                    tag.loc.vector <- loc.vector[flag]
                    tag.loc.vector <- sort(tag.loc.vector, decreasing=F)
    
                    if(weight==0) {
                        tag.correl.vector <- rep(1, Nh)
                    }else if(weight==1){
                        tag.correl.vector <- rank.score.sorted[tag.loc.vector]
                        tag.correl.vector <- abs(tag.correl.vector)
                    }else{
                        tag.correl.vector <- rank.score.sorted[tag.loc.vector] ** weight
                        tag.correl.vector <- abs(tag.correl.vector)    
                    }

                    norm.tag <- 1.0/sum(tag.correl.vector)
                    tag.correl.vector <- tag.correl.vector * norm.tag
                    norm.no.tag <- 1.0/Nm
                    tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
                    tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
                    tag.diff.vector <- tag.diff.vector * norm.no.tag
                    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
                    valley.res.vector <- peak.res.vector - tag.correl.vector
                    max.ES <- max(peak.res.vector)
                    min.ES <- min(valley.res.vector)
                    es.expected[k,r] <- signif(ifelse(max.ES>abs(min.ES), max.ES, min.ES), digits=5)
                }

            }
        }
        
        ##########
        ## P-value
        if(sigTail=="one-tail"){
            # one-tail
            pES <- sapply(1:length(es.observed), function(x){
                #sum(es.observed[x]<=es.expected[x,])/nperm
                signif(ifelse(es.observed[x]>=0, sum(es.observed[x]<=es.expected[x,])/nperm, sum(es.observed[x]>=es.expected[x,])/nperm), digits=5)
            })
        }else{
            # two-tails
            pES <- sapply(1:length(es.observed), function(x){     
                if(sum(es.expected[x,]>=0)==nperm){
                    # all positives
                    signif(sum(es.observed[x]<=es.expected[x,])/nperm, digits=5)
                }else if(sum(es.expected[x,]>=0)==0){
                    # all negatives
                    signif(sum(es.observed[x]>=es.expected[x,])/nperm, digits=5)
                }else{
                    pos.phi <- es.expected[x,es.expected[x,]>=0]
                    neg.phi <- es.expected[x,es.expected[x,]<0]
                    signif(ifelse(es.observed[x]>=0, sum(es.observed[x]<=pos.phi)/length(pos.phi), sum(es.observed[x]>=neg.phi)/length(neg.phi)), digits=5)
                }
            })
        }

        ## adjusted p-value
        adjP <- p.adjust(pES, method="BH")
        
        ##########        
        ## normalised ES score
        if(sigTail=="one-tail"){
            nES.observed <- es.observed / rowMeans(es.expected)
            nES.expected <- es.expected / rowMeans(es.expected)
        }else{
            # two-tails
            nES.observed <- matrix(0, nrow=nSet, ncol=1)
            nES.expected <- matrix(0, nrow=nSet, ncol=nperm)
            for(i in 1:length(es.observed)){
                pos.flag <- es.expected[i,]>=0
                neg.flag <- es.expected[i,]<0

                if(sum(pos.flag)==nperm){
                    # all positives
                    # observed
                    nES.observed[i] <- signif(es.observed[i]/mean(es.expected[i,]), digits=5)
                    # expected
                    nES.expected[i,] <- signif(es.expected[i,]/mean(es.expected[i,]), digits=5)
                    
                }else if(sum(pos.flag)==0){
                    # all negatives
                    # observed
                    nES.observed[i] <- signif(es.observed[i]/abs(mean(es.expected[i,])), digits=5)
                    # expected
                    nES.expected[i,] <- signif(es.expected[i,]/abs(mean(es.expected[i,])), digits=5)
                }else{
                    pos.m <- mean(es.expected[i,pos.flag])
                    neg.m <- abs(mean(es.expected[i,neg.flag]))
                    # observed
                    nES.observed[i] <- signif(ifelse(es.observed[i]>=0, es.observed[i]/pos.m, es.observed[i]/neg.m), digits=5)
                    # expected
                    nES.expected[i,pos.flag] <- signif(es.expected[i,pos.flag]/pos.m, digits=5)
                    nES.expected[i,neg.flag] <- signif(es.expected[i,neg.flag]/neg.m, digits=5)
                }
            }
        }
        
        ##########
        ## FWER
        if(sigTail=="one-tail"){
            # one-tail
            max.expected.flag <- apply(abs(nES.expected), 2, which.max)
            max.expected <- sapply(1:length(max.expected.flag), function(x) nES.expected[max.expected.flag[x],x])
            FWER <- sapply(1:length(nES.observed), function(x){
                ifelse(nES.observed[x]>=0, sum(max.expected>=nES.observed[x])/nperm, sum(max.expected<=nES.observed[x])/nperm)
            })
        }else{
            # two-tails
            max.expected.flag <- apply(abs(nES.expected), 2, which.max)
            max.expected <- sapply(1:length(max.expected.flag), function(x) nES.expected[max.expected.flag[x],x])
            max.expected.pos <- max.expected[max.expected>=0]
            max.expected.neg <- max.expected[max.expected<0]
            FWER <- sapply(1:length(nES.observed), function(x){
                ifelse(nES.observed[x]>=0, sum(max.expected.pos>=nES.observed[x])/length(max.expected.pos), sum(max.expected.neg<=nES.observed[x])/length(max.expected.neg))
            })
        }
        
        ##########        
        ## FDR
        if(sigTail=="one-tail"){
            # one-tail
            FDR <- sapply(1:length(nES.observed), function(x){
                median(apply(nES.expected>=nES.observed[x],2,sum)) / sum(nES.observed>=nES.observed[x])
            })
            FDR <- ifelse(FDR>1,1,FDR)
        }else{
            # two-tails
            nES.observed.pos <- nES.observed[nES.observed>=0]
            nES.observed.neg <- nES.observed[nES.observed<0]
            FDR <- sapply(1:length(nES.observed), function(x){
                if(nES.observed[x] >= 0){
                    called <- sum(nES.observed.pos>=nES.observed[x])
                    ifelse(called>0, median(apply(nES.expected>=nES.observed[x],2,sum)) / called, 1)
                }else{
                    called <- sum(nES.observed.neg<=nES.observed[x])
                    ifelse(called>0, median(apply(nES.expected<=nES.observed[x],2,sum)) / called, 1)
                }
            })
            FDR <- ifelse(FDR>1, 1, FDR)
        }
        
        ##########        
        ## FDR q-value
        flag <- order(nES.observed, decreasing=T)
        # for origin index
        Orig.index <- seq(1,nSet)
        Orig.index <- order(Orig.index[flag], decreasing=F)
        # FDR being sorted according to nES.observed
        qES.sorted <- FDR[flag]
        # for positive part
        pos.part <- length(nES.observed[nES.observed>=0])
        if(pos.part > 0){
            tmp.min <- qES.sorted[pos.part]
            # (FDR) from higest to the lowest
            for (k in seq(pos.part, 1, -1)) {
                if (qES.sorted[k] < tmp.min) {
                    tmp.min <- qES.sorted[k]
                }else{
                    qES.sorted[k] <- tmp.min
                }
            }
        }
        # for negative part
        if(pos.part < nSet){
            neg.part <- pos.part + 1
            tmp.min <- qES.sorted[pos.part]
            for (k in seq(neg.part, nSet)) {
                if (qES.sorted[k] < tmp.min) {
                    tmp.min <- qES.sorted[k]
                }else{
                    qES.sorted[k] <- tmp.min
                }
            }
        }
        # put back to the origin
        qES <- qES.sorted[Orig.index]

        ############################
        SS.es[,j] <- es.observed
        SS.nes[,j] <- nES.observed
        SS.pvalue[,j] <- pES
        SS.adjp[,j] <- adjP
        SS.fwer[,j] <- FWER
        SS.fdr[,j] <- FDR
        SS.qvalue[,j] <- qES
    }
    
    ## globally adjusted p value for Multiple Comparisons
    vec <- c(SS.pvalue)
    gadjp <- p.adjust(vec, method="BH")
    SS.gadjp <- matrix(gadjp, nrow=nSet, ncol=nSample)
    colnames(SS.gadjp) <- colnames(data)
    rownames(SS.gadjp) <- names(gs)

    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("End at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    eTerm <- list(set_info = set_info,
                  gs       = gs,
                  data     = data,
                  es       = SS.es,
                  nes      = SS.nes,
                  pvalue   = SS.pvalue,
                  adjp     = SS.adjp,
                  gadjp    = SS.gadjp,
                  fwer     = SS.fwer,
                  fdr      = SS.fdr,
                  qvalue   = SS.qvalue,
                  call     = match.call()
                 )
    class(eTerm) <- "eTerm"
    
    invisible(eTerm)
}