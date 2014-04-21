# This is a demo for TCGA mutational profile dataset from Kandoth et al
# 
# This dataset is available from TCGA (see <a href="http://www.ncbi.nlm.nih.gov/pubmed/24132290" target="24132290">http://www.ncbi.nlm.nih.gov/pubmed/24132290</a>), containing somatic mutational profiles for 3096 cancer samples with survival data. These cancer samples belong to one of 12 major cancer types, including breast adenocarcinoma (BRCA), lung adenocarcinoma (LUAD), lung squamous cell carcinoma (LUSC), uterine corpus endometrial carcinoma (UCEC), glioblastoma multiforme (GBM), head and neck squamous cell carcinoma (HNSC), colon and rectal carcinoma (COAD/READ), bladder urothelial carcinoma (BLCA), kidney renal clear cell carcinoma (KIRC), ovarian serous carcinoma (OV) and acute myeloid leukaemia (LAML). For each patient sample, somatic mutations are represented as a profile of  states on genes, where non-zero entry indicates a gene for which how many mutations have occurred in the tumor relative to germ line. The dataset is provided as an 'ExpressionSet' object.
## assayData: exprs(TCGA_mutations), a matrix of 19171 genes X 3096 samples;
## phenoData: pData(TCGA_mutations), variables describing sample phenotypes (i.e. columns in assayData), including clinical/survival information about samples: "time" (i.e. survival time in days), "status" (i.e., survival status: 0=alive; 1=dead), "Age" (the patient age in years), "Gender" (the patient gender: male/female), "TCGA_tumor_type", "Tumor_stage", "Tumor_grade"
## featureData: fData(TCGA_mutations), variables describing features (i.e. rows in assayData), including information about features/genes: "EntrezID" for gene EntrezID, "Symbol" for gene symbol, "Desc" for gene description, "Synonyms" for gene symbol alias
###############################################################################
library(dnet)

# Load or install packages specifically used in this demo
list.pkg <- c("affy", "survival")
source("http://bioconductor.org/biocLite.R")
for(pkg in list.pkg){
    if(!require(pkg, character.only=T)){
        biocLite(pkg)
        lapply(pkg, library, character.only=T)
    }
}

# load an "ExpressionSet" object
load(url("http://dnet.r-forge.r-project.org/data/Datasets/TCGA_mutations.RData"))
#load("RData_Rd/data/Datasets/TCGA_mutations.RData")
eset <- TCGA_mutations
# extract information about phenotype data
pd <- pData(eset)
pd[1:5,]
# extract information about feature/gene data
fd <- fData(eset)
fd[1:5,]
# extract information about mutational data
md <- exprs(eset)
md[1:5,1:5]
# number of samples for each cancer type
table(pData(eset)$TCGA_tumor_type)
tumor_type <- sort(unique(pData(eset)$TCGA_tumor_type))

# Kaplan-Meier survival curves for individual tumor types
fit <- survfit(Surv(time, status) ~ TCGA_tumor_type, data=pd, type=c("kaplan-meier","fh")[1])
plot(fit, xscale=365.25, xlab = "Years", ylab="Survival", lty=1:2, col=rainbow(length(tumor_type))) 
legend("topright", tumor_type, lty=1:2, col=rainbow(length(tumor_type)))

############################################################################
# Survival analysis across tumor types using Cox proportional hazards model
# Cox regression yields an equation for the hazard/risk as a function of several explanatory variables
# Except for the gene mutational data in subject,  other explanatory variables (or called covariates) include: age, gender, and tumor type
## only those genes with mutations at least 1% of samples will be analysed
flag <- sapply(1:nrow(md), function(i) ifelse(sum(md[i,]!=0)>=0.01*ncol(md), T, F))
esetGene <- eset[flag, ]
md_selected <- md[flag,]
## survival analysis to obtain hazard ratio (LR) and pvaules
gene_signif <- matrix(1, nrow=nrow(md_selected), ncol=2)
rownames(gene_signif) <- rownames(md_selected)
colnames(gene_signif) <- c("LR", "pvalue")
for(i in 1:nrow(md_selected)){
    if(i %% 100==0 || i==1 || i==nrow(md_selected)) message(sprintf("%d genes processed (%s)", i, as.character(Sys.time())), appendLF=T)
    ## fit a Cox proportional hazards model
    #gene_mut <- 1*(md_selected[i,]!=0)
    gene_mut <- md_selected[i,]
    data <- cbind(pd,gene=gene_mut)
    fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + gene, data=data)
    res <- as.matrix(anova(fit))
    ## 2nd: likelyhood ratio; 4th: pvalue
    gene_signif[i,] <- res[5,c(2,4)]
}
LR <- gene_signif[,1]
pvals <- gene_signif[,2]

############################################################################
## survival analysis of each tumor type
type_gene_signif <- list()
for (k in 1:length(tumor_type)){
    type <- tumor_type[k]

    ind <- match(pData(eset)$TCGA_tumor_type, type)
    type_md_selected <- md_selected[, !is.na(ind)]
    
    tmp <- matrix(1, nrow=nrow(type_md_selected), ncol=2)
    rownames(tmp) <- rownames(type_md_selected)
    colnames(tmp) <- c("LR", "pvalue")
    for(i in 1:nrow(type_md_selected)){
        if(i %% 100==0 || i==1 || i==nrow(type_md_selected)) message(sprintf("%s: %d genes processed (%s)", type, i, as.character(Sys.time())), appendLF=T)
        ## fit a Cox proportional hazards model
        gene_mut <- type_md_selected[i,]
        if(sum(gene_mut!=0)>=0.01*length(gene_mut)){
            data <- cbind(pd[!is.na(ind),],gene=gene_mut)
            fit <- coxph(formula=Surv(time,status) ~ Age + Gender + gene, data=data)
            res <- as.matrix(anova(fit))
            ## 2nd: likelyhood ratio; 4th: pvalue
            tmp[i,] <- res[4,c(2,4)]
        }else{
            tmp[i,] <- c(1,1)
        }
    }
    type_gene_signif[[k]] <- tmp
}
names(type_gene_signif) <- tumor_type

type_LR <- matrix(1, nrow=nrow(md_selected), ncol=length(tumor_type))
rownames(type_LR) <- rownames(md_selected)
colnames(type_LR) <- tumor_type
type_pvals <- matrix(1, nrow=nrow(md_selected), ncol=length(tumor_type))
rownames(type_pvals) <- rownames(md_selected)
colnames(type_pvals) <- tumor_type
for(k in 1:length(type_gene_signif)){
    tmp <- type_gene_signif[[k]]
    type_LR[,k] <- tmp[,1]
    type_pvals[,k] <- tmp[,2]
}



############################################################################
# An igraph object that contains a functional protein association network in human. The network is extracted from the STRING database (version 9.1). Only those associations with medium confidence (score>=400) are retained.
load(url("http://dnet.r-forge.r-project.org/data/Hs/org.Hs.string.RData"))
#load("RData_Rd/data/Hs/org.Hs.string.RData")
org.Hs.string
# restrict to those edges with high confidence (score>=700)
g <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=700])
g


# extract network that only contains genes in eset
ind <- match(V(g)$symbol, names(pvals))
## for extracted graph
nodes_mapped <- V(g)$name[!is.na(ind)]
network <- dNetInduce(g, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=T)
V(network)$name <- V(network)$symbol
#E(network)$weight <- E(network)$combined_score
network

# Identification of gene-active subnetwork
# 2) identification of gene-active subnetwork
#pvals <- pvals[V(network)$symbol]
## restrict the identified subnetwork to have the node size of 40 or so
#g <- dNetPipeline(g=network, pval=pvals, method="fdr", nsize=40)
g <- dNetPipeline(g=network, pval=pvals, method="fdr", fdr=3e-02)
g


############################################################################
# Estimate the robustness of the identified subnetwork based on jackknife sampling strategy
B <- 1000
gene_pvals_random <- matrix(1, nrow=nrow(md_selected), ncol=B)
rownames(gene_pvals_random) <- rownames(md_selected)
colnames(gene_pvals_random) <- 1:B
for(j in 1:B){
    for(i in 1:nrow(md_selected)){
        if(i %% 1000==0 || i==1 || i==nrow(md_selected)) message(sprintf("Round %d: %d genes processed (%s)", j, i, as.character(Sys.time())), appendLF=T)
        # at least 1% of patients
        gene_mut <- md_selected[i,]
        ## fit a Cox proportional hazards model
        data <- cbind(pd,gene=gene_mut)
        ## do sampling (80% samples without replacement)
        flag <- sample(1:nrow(data), size=0.8*nrow(data), replace=F)
        data <- data[flag,]
        ## fit
        fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + gene, data=data)
        res <- as.matrix(anova(fit))
        ## 2nd: likelyhood ratio; 4th: pvalue
        gene_pvals_random[i,j] <- res[5,4]
    }
}
save(list=c("gene_pvals_random","modules"), file="gene_pvals_random.RData")

modules <- list()
for(j in 1:B){
    message(sprintf("###Jacknife resampling: %d replicate", j), appendLF=T)
    # For each jackknife replicate, recalculate survival p-values, which is used to calculate node scores. The same fdr as before is used
    #modules[[j]] <- dNetPipeline(g=network, pval=gene_pvals_random[,j], method="fdr", fdr=3e-02)
    modules[[j]] <- dNetPipeline(g=network, pval=gene_pvals_random[,j], method="fdr", nsize=vcount(g))
}

# Estimate the robustness of the identified subnetwork by removal of one tumor type (defined by tissue of origin) and re-identification of the subnetwork based on the remaining tumor types
gene_pvals_loo <- matrix(1, nrow=nrow(md_selected), ncol=length(tumor_type))
rownames(gene_pvals_loo) <- rownames(md_selected)
colnames(gene_pvals_loo) <- tumor_type
for(k in 1:length(tumor_type)){
    type <- tumor_type[k]
    ind <- match(pData(eset)$TCGA_tumor_type, type)
    type_md_selected <- md_selected[, is.na(ind)]
    
    for(i in 1:nrow(type_md_selected)){
        if(i %% 1000==0 || i==1 || i==nrow(type_md_selected)) message(sprintf("%s: %d genes processed (%s)", type, i, as.character(Sys.time())), appendLF=T)
        gene_mut <- type_md_selected[i,]
        ## fit a Cox proportional hazards model
        data <- cbind(pd[is.na(ind),],gene=gene_mut)
        ## fit
        fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + gene, data=data)
        res <- as.matrix(anova(fit))
        ## 2nd: likelyhood ratio; 4th: pvalue
        gene_pvals_loo[i,k] <- res[5,4]
    }
}

modules <- list()
for(j in 1:length(tumor_type)){
    message(sprintf("###Jacknife resampling: %d replicate", j), appendLF=T)
    # For each leave-one-out, recalculate survival p-values, which is used to calculate node scores
    modules[[j]] <- dNetPipeline(g=network, pval=gene_pvals_loo[,j], method="fdr", nsize=vcount(g))
}

############################################################################


# append the confidence information from the source graphs into the target graph
cmodule <- dNetConfidence(target=g, sources=modules, plot=F)
E(cmodule)$edgeConfidence <- ceiling(E(cmodule)$edgeConfidence)
# visualise the confidence target graph
visNet(cmodule, layout=glayout, edge.width=E(cmodule)$edgeConfidence/10, edge.label=E(cmodule)$edgeConfidence, edge.label.cex=0.6)

# Cumulative distribution of the edges as a function of the confidence score (the times an edge repeatedly identified)
# half of edges can be repeatedly identified at the confidence level of 60%, after removing one tumor type each time
# the network identified is not dominated by one single tumor type
plot(ecdf(E(cmodule)$edgeConfidence), verticals=TRUE, do.points=F, xlim=c(0,100), ylab="Fraction of edges", xlab="Confidence scores", lwd=3)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)



# 3) visualisation of the gene-active subnetwork itself
## the layout of the network visualisation (fixed in different visuals) 
#glayout <- layout.fruchterman.reingold(g)
## color nodes according to communities (identified via a spin-glass model and simulated annealing)
com <- spinglass.community(g, spins=25)
com$csize <- sapply(1:length(com),function(x) sum(com$membership==x))
vgroups <- com$membership
colormap <- "yellow-darkorange"
palette.name <- visColormap(colormap=colormap)
mcolors <- palette.name(length(com))
vcolors <- mcolors[vgroups]
com$significance <- dCommSignif(g, com)
## node sizes according to degrees
vdegrees <- igraph::degree(g)
## highlight different communities
mark.groups <- communities(com)
mark.col <- visColoralpha(mcolors, alpha=0.2)
mark.border <- visColoralpha(mcolors, alpha=0.2)
edge.color <- c("#C0C0C0", "#000000")[crossing(com,g)+1]
edge.color <- visColoralpha(edge.color, alpha=0.5)
## visualise the subnetwrok
#visNet(g, glayout=glayout, vertex.label=V(g)$geneSymbol, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color, vertex.label.color="darkblue", vertex.label.dist=0, vertex.label.font=2)

wth <- 3600
png("a.png", width=wth, height=wth, res=wth*72/480)
visNet(g, glayout=glayout, vertex.label=V(g)$geneSymbol, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color, newpage=F, vertex.label.color="blue", vertex.label.dist=0.4, vertex.label.font=2, edge.width=E(cmodule)$edgeConfidence/10)
legend_name <- paste("C",1:length(mcolors)," (n=",com$csize,", pval=",signif(com$significance,digits=2),")",sep='')
legend("topleft", legend=legend_name, fill=mcolors, bty="n", cex=0.6)
dev.off()


esetGene <- eset
# For each gene, calculate number of samples (within a tumor type) having the mutated gene 
flag_net_mutated <- sapply(tumor_type, function(x) {
    e <- esetGene[, which(pData(esetGene)$TCGA_tumor_type==x)]
    tmp <- apply(exprs(e)[V(g)$name,]!=0, 2, sum)
})
for(i in 1:length(flag_net_mutated)){
    tmp <- flag_net_mutated[[i]]
    data <- cbind(pd[names(tmp),], net=tmp!=0)
    # Kaplan-Meier survival curves for individual tumor types
    fit <- survfit(Surv(time, status) ~ net, data=data, type=c("kaplan-meier","fh")[1])
    dev.new()
    plot(fit, xscale=365.25, xlab = "Years", ylab="Survival", lty=1:2, col=rainbow(2)) 
    text_legend <- paste(names(flag_net_mutated[i]), c(" without net"," with net"), sep="")
    legend("topright", text_legend, lty=1:2, col=rainbow(2))
}


# fit a Cox proportional hazards model using the subnetwork
## for the whole network
data_g <- t(md[V(g)$name,])
data_g <- apply(data_g!=0, 1, sum)
data <- cbind(pd, net=data_g)
fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + net, data=data)
res <- as.matrix(anova(fit))
## 2nd: likelyhood ratio; 4th: pvalue
LR_g <- res[5,2]
pvals_g <- res[5,4]
## for the cumulative nodes from the network
cg_names <- names(sort(LR[V(g)$name], decreasing=T))
cg_signif <- matrix(1, nrow=length(cg_names), ncol=2)
rownames(cg_signif) <- cg_names
colnames(cg_signif) <- c("LR", "pvalue")
for(i in 1:length(cg_names)){
    data_g <- t(md[cg_names[1:i],])
    if(i!=1){
        data_g <- apply(data_g!=0, 1, sum)
    }else{
        data_g <- as.vector(data_g!=0)
    }
    data <- cbind(pd, cnet=data_g)
    fit <- coxph(formula=Surv(time,status) ~ Age + Gender + TCGA_tumor_type + cnet, data=data)
    res <- as.matrix(anova(fit))
    cg_signif[i,] <- res[5,c(2,4)]
}
cg_signif[cg_signif[,2]==0,2] <- min(cg_signif[cg_signif[,2]!=0,2])

# median: 0.61, 6.2 66.1
bp.LR.list <- list(All=LR, Neti=LR[cg_names], Netc=cg_signif[2:nrow(cg_signif),1])
par(las=1, mar=c(5,8,4,2)) # all axis labels horizontal
boxplot(bp.LR.list, outline=F, horizontal=F, names=c("All genes", "Individual genes\n in the network", "Combined genes\n in the network"), col=c("red","green","blue"), border=par("fg"), ylab="Hazard ratio", log="y", ylim=c(0.1,100), yaxt="n",xaxt="n")
axis(2, at=c(0.1,1,10,100), labels=c(0.1,1,10,100), las=2)
axis(1, at=1:3, labels=c("All genes", "Individual genes\n in the network", "Combined genes\n in the network"), las=2)

## Two-sample Kolmogorov-Smirnov test
# p-value = p-value < 2.2e-16
stats::ks.test(x=LR, y=LR[cg_names], alternative=c("two.sided","less", "greater")[1], exact=NULL)
# p-value = p-value < 2.2e-16
stats::ks.test(x=LR[cg_names], y=cg_signif[2:nrow(cg_signif),1], alternative=c("two.sided","less", "greater")[1], exact=NULL)
# p-value = p-value < 2.2e-16
stats::ks.test(x=LR, y=cg_signif[2:nrow(cg_signif),1], alternative=c("two.sided","less", "greater")[1], exact=NULL)


plot(cg_signif[,1])
lines(LR[cg_names])

x <- 1:length(cg_names)
y <- LR[cg_names]
y <- cg_signif[,1]
df <- data.frame(x,y)
# model
mod <- lm(y ~ x, data = df)
# predicts + interval
newx <- df$x
preds <- predict(mod, newdata=data.frame(x=newx), interval='confidence')
# plot
plot(y ~ x, data=df, type='l')
# add fill
polygon(c(rev(newx),newx), c(rev(preds[,3]), preds[,2]), col='grey80', border=NA)
# model
abline(mod)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')



# 4) visualisation of the gene-active subnetwork overlaid by the node/gene score
max_colorbar <- ceiling(quantile(abs(V(g)$score),0.5))
visNet(g, glayout=glayout, pattern=V(g)$score, zlim=c(-1*max_colorbar,max_colorbar), vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 5) visualisation of the gene-active subnetwork overlaid by p-values
colormap <- "darkgreen-lightgreen-lightpink-darkred"
logFC <- -1*log10(pvals[V(g)$name])
visNet(g, glayout=glayout, pattern=logFC, colormap=colormap, vertex.shape="sphere", zlim=c(0,2), mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 6) visualisation of the gene-active subnetwork overlaid by log-LR
colormap <- "darkgreen-lightgreen-lightpink-darkred"
logFC <- log2(LR[V(g)$name])
visNet(g, glayout=glayout, pattern=logFC, colormap=colormap, vertex.shape="sphere", zlim=c(-2,2), mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 7) Network-based sample classifications and visualisations on 2D sample landscape
# it uses the gene-active subnetwork overlaid by all replication timing data
frac_mutated <- sapply(tumor_type, function(x) {
    e <- eset[, which(pData(eset)$TCGA_tumor_type==x)]
    apply(exprs(e)!=0,1,sum)/ncol(e)
})
rownames(frac_mutated) <- fData(esetGene)$Symbol
frac_mutated[1:10,]

data <- frac_mutated[V(g)$name,]
sReorder <- dNetReorder(g, data, feature="edge", node.normalise="degree", amplifier=3, metric="none")
visNetReorder(g=g, data=data, sReorder=sReorder, height=ceiling(sqrt(ncol(data)))*3, newpage=T, glayout=glayout, colormap=colormap, vertex.label=NA,vertex.shape="sphere", vertex.size=16,mtext.cex=0.8,border.color="888888", zlim=c(0,0.10), mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# 8) heatmap of replication timing data in the subnetwork
visHeatmapAdv(data, colormap=colormap, zlim=c(0,0.15), KeyValueName="Fraction of mutations")

# 9) output the subnetwork and their replication timing data
## Write the subnetwork into a SIF-formatted file (Simple Interaction File)
sif <- data.frame(source=get.edgelist(g)[,1], type="interaction", target=get.edgelist(g)[,2])
write.table(sif, file=paste("Survival_TCGA.sif", sep=""), quote=F, row.names=F,col.names=F,sep="\t")
## Output the corresponding replication timing data
hmap <- data.frame(Symbol=rownames(data), data)
write.table(hmap, file=paste("Survival_TCGA.txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")

# 9) enrichment analysis for genes in the subnetwork
## get a list of genes in the subnetwork
data <- V(g)$name
data

## 9a) GOBP enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="GOBP")
## visualise the top significant terms in the GOBP heirarchy
## first, load the GOBP ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOBP.RData"))
g <- ig.GOBP
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$adjp)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[,2:3], nSet=sapply(eTerm$gs,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,]
## induce the shortest paths (one for each significant term) to the ontology root
subg <- dDAGinduce(g, nodes_query, path.mode="shortest_paths")
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", node.attrs=list(color=nodes.highlight))

## 9b) GOMF enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="GOMF")
## visualise the top significant terms in the GOMF heirarchy
## first, load the GOMF ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOMF.RData"))
g <- ig.GOMF
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$adjp)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[,2:3], nSet=sapply(eTerm$gs,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,]
## induce the shortest paths (one for each significant term) to the ontology root
subg <- dDAGinduce(g, nodes_query, path.mode="shortest_paths")
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", node.attrs=list(color=nodes.highlight))

## 9c) MP enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="MP", ontology.algorithm="elim")
## visualise the top significant terms in the MP heirarchy
## first, load the MP ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.MP.RData"))
g <- ig.MP
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$pvalue)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[nodes_query,2:3], cbind(zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,])
## induce all possible paths to the ontology root
subg <- dDAGinduce(g, nodes_query)
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info=c("none","term_id","term_name","both","full_term_name")[5], layout.orientation=c("left_right","top_bottom","bottom_top","right_left")[1], node.attrs=list(color=nodes.highlight))

## 9d) DO enrichment analysis
eTerm <- dEnricher(data, identity="symbol", genome="Hs", ontology="DO", ontology.algorithm="pc")
## visualise the top significant terms in the DO heirarchy
## first, load the DO ontology
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.DO.RData"))
g <- ig.DO
## select the top most significant 10 terms
nodes_query <- names(sort(eTerm$adjp)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
cbind(eTerm$set_info[nodes_query,2:3], cbind(zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,])
## induce all possible shortest paths to the ontology root
subg <- dDAGinduce(g, nodes_query)
## color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", zlim=c(0,4), node.attrs=list(color=nodes.highlight))


##########################################################################################
## 9e) PS enrichment analysis
data <- V(g)$name
eTerm <- dEnricher(data, test="HypergeoTest", identity="symbol", genome="Hs", ontology="PS2", RData.location="RData_Rd/data")
## Look at the evolution relevance along the path to the eukaryotic common ancestor
cbind(eTerm$set_info[,2:3], nSet=sapply(eTerm$gs,length), nOverlap=sapply(eTerm$overlap,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)

out <- as.data.frame(cbind(eTerm$set_info[,2:3], nSet=sapply(eTerm$gs,length), nOverlap=sapply(eTerm$overlap,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp))
write.table(out, file=paste("PS.txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")

#load(url("http://dnet.r-forge.r-project.org/data/Hs/org.Hs.eg.RData"))
load("RData_Rd/data/Hs/org.Hs.eg.RData")
names(org.Hs.eg)
gene_info <- org.Hs.eg$gene_info
entrez <- unlist(eTerm$overlap[6], use.names=F)
gene_info[match(entrez,rownames(gene_info)),]

data <- frac_mutated[V(g)$name,]

tree_bs <- visTreeBootstrap(t(data), nodelabels.arg=list(cex=0.4,bg="white-pink-violet"), metric=c("euclidean","pearson","spearman","cos","manhattan","kendall","mi")[3], num.bootstrap=2000, plot.phylo.arg=list(type="p", direction="downwards", cex=0.4, edge.width=1.2), show.tip.label=F)
flag <- match(tree_bs$tip.label, colnames(data))

base <- sapply(eTerm$overlap, function(x){
    as.character(gene_info[match(x,rownames(gene_info)),2])
})
if(1){
    ## reordering via hierarchical clustering
    cluster_order <- matrix(1, nrow=length(base))
    base_order <- matrix(1, nrow=length(base))
    for(i in 1:length(base)){
        tmp <- base[[i]]
        ind <- match(tmp, rownames(data))
        if(length(ind)>0){
            base_order[ind] <- i
            tmpD <- data[ind,]
            if(length(tmp) != 1){
                distance <- as.dist(sDistance(tmpD, metric="pearson"))
                cluster <- hclust(distance, method="average")
                cluster_order[ind] <- cluster$order
            }else if(length(tmp) == 1){
                cluster_order[ind] <- 1
            }
        }
    }
    ## contruct data frame including 1st column for temporary index, 2nd for cluster order, 3rd for base/cluster ID
    df <- data.frame(ind=1:nrow(data), cluster_order, base_order)
    # order by: first base, then hexagon
    ordering <- df[order(base_order,cluster_order),]$ind
}

RowSideColors <- sapply(1:length(base), function(x) base_order==x)
RowSideColors <- t(RowSideColors)
rslab <- ifelse(eTerm$adjp<0.05," (FDR<0.05)","")
rslab <- paste(gsub(".*:","",eTerm$set_info$name), rslab, sep="")
rownames(RowSideColors) <- rslab
colnames(RowSideColors) <- rownames(data)
RowSideColors <- ifelse(RowSideColors==T, "gray","white")
RowSideColors <- RowSideColors[, ordering]

base_order1 <- base_order[ordering]
basesep_index <- sapply(unique(base_order1), function(x) which(base_order1[length(base_order1):1]==x)[1])
basesep_index <- basesep_index[1:length(basesep_index)-1]

labRow <- sapply(pvals[match(V(g)$name, names(pvals))], function(x){
    if(x < 0.005){
        " ***"
    }else if(x < 0.01){
        " **"
    }else if(x<0.05){
        " *"
    }else{
        ""
    }
})
labRow <- paste(rownames(data), labRow, sep="")

visHeatmapAdv(data=data[ordering,flag], Rowv=F, Colv=F, colormap="lightyellow-orange", zlim=c(0,0.12), keysize=1.5, RowSideColors=RowSideColors, RowSideWidth=2, RowSideLabelLocation="top", add.expr=abline(h=(basesep_index-0.5), lty=2,lwd=1,col="black"), offsetRow=-0.5, labRow=labRow[ordering], KeyValueName="Mutation frequency", margins = c(6,6))


# Cross-tumor mutational ubiquity
ind <- match(V(g)$name, rownames(denseness))
net_denseness <- denseness[ind]
net_denseness <- net_denseness[ordering]
names(net_denseness) <- rownames(data)[ordering]
library(beeswarm)
par(las=2, mar=c(12,8,4,2)) # all axis labels horizontal
beeswarm(net_denseness ~ base_order1, col=4, pch=16, horizontal=F, "ylab"="Cross-tumor mutational ubiquity", "xlab"="", labels="", ylim=c(0,1))
lbls <- eTerm$set_info$name[unique(base_order1)]
lbls <- gsub(".*:","",lbls)
boxplot(net_denseness ~ base_order1, add=T, horizontal=F, names=lbls)
# p-value = 0.05472
stats::ks.test(x=net_denseness[base_order1==6], y=net_denseness, alternative=c("two.sided","less", "greater")[1], exact=NULL)
# p-value = 0.006488
stats::ks.test(x=net_denseness[base_order1==6], y=net_denseness[base_order1<6], alternative=c("two.sided","less", "greater")[1], exact=NULL)
# p-value = 0.2121
stats::ks.test(x=net_denseness[base_order1==6], y=net_denseness[base_order1>6], alternative=c("two.sided","less", "greater")[1], exact=NULL)


###############
# data output
rnames <- rownames(data[ordering,flag])
output <- list()
ind <- match(rnames, gene_info$Symbol)
output$gene_info <- gene_info[ind,c(1,2,3,4,5)]
output$ancestor <- as.character(eTerm$set_info$name[base_order][ordering])
ind <- match(rnames, names(pvals))
output$pval_overall <- pvals[ind]
output$HR_overall <- LR[ind]
output$freq <- data[ordering,flag]
ind <- match(rnames, rownames(denseness))
output$denseness <- denseness[ind]
ind <- match(rnames, rownames(type_pvals))
output$pval_type <- type_pvals[ind, flag]
output$HR_type <- type_LR[ind, flag]
output$community <- vgroups[ordering]

out <- as.data.frame(output)
write.table(out, file=paste("out.txt", sep=""), quote=F, row.names=F,col.names=T,sep="\t")
##########################################################################################



data <- V(g)$name
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="SF", RData.location="RData_Rd/data")

eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC2KEGG", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC2CP", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC2REACTOME", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC3TFT", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC3MIR", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC5BP", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC5MF", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC6", RData.location="RData_Rd/data")
eTerm <- dEnricher(data, test=c("FisherTest","HypergeoTest","BinomialTest")[2], identity="symbol", genome="Hs", ontology="MsigdbC7", RData.location="RData_Rd/data")
nodes_query <- names(sort(eTerm$pvalue)[1:20])
cbind(eTerm$set_info[nodes_query,2:3], cbind(nSet=sapply(eTerm$gs,length), nOverlap=sapply(eTerm$overlap,length), zscore=eTerm$zscore, pvalue=eTerm$pvalue, adjp=eTerm$adjp)[nodes_query,])



esetGene <- eset
# For each gene, calculate number of samples (within a tumor type) having the mutated gene 
num_mutated <- sapply(tumor_type, function(x) {
    e <- esetGene[, which(pData(esetGene)$TCGA_tumor_type==x)]
    apply(exprs(e)!=0,1,sum)
})
rownames(num_mutated) <- fData(esetGene)$Symbol
num_mutated[1:10,]

# For each gene, calculate fraction of samples (within a tumor type) having the mutated gene 
frac_mutated <- sapply(tumor_type, function(x) {
    e <- esetGene[, which(pData(esetGene)$TCGA_tumor_type==x)]
    apply(exprs(e)!=0,1,sum)/ncol(e)
})
rownames(frac_mutated) <- fData(esetGene)$Symbol
frac_mutated[1:10,]

# define the "mutation sparseness" of genes in terms of a vector which stores the fraction of samples (within a tumor type) having the mutated gene 
# sparseness for a vector is: 1) one if the vector contains only a single non-zero value; 2) zero if and only if all elements are equal; 3) otherwise, the value interpolates smoothly between the two extremes
sparseness <- sapply(1:nrow(frac_mutated), function(i){
    v <- frac_mutated[i,]
    n <- length(v)
    norm1 <- sum(abs(v))
    norm2 <- sqrt(sum(v^2))
    (sqrt(n)-norm1/norm2) / (sqrt(n)-1)
})
sparseness <- matrix(sparseness, ncol=1)
rownames(sparseness) <- rownames(frac_mutated)
# derive the "mutation denseness/ubiquity" of genes: mutational fraction with the same type, and fraction consistent across different types
denseness <- 1- sparseness
hist(denseness,100, xlab="Cross-tumor mutation ubiquity", xlim=c(0,1))


eTerm <- dGSEA(data=cbind(frac_mutated,denseness), identity="symbol", genome="Hs", ontology="Customised", customised.genesets=V(g)$name, weight=0, nperm=5000, RData.location="RData_Rd/data")

frac_pvalue <- as.vector(eTerm$pvalue)
frac_fdr <- stats::p.adjust(frac_pvalue, method="BH")
frac_nes <- as.vector(eTerm$nes)
frac_es <- as.vector(eTerm$es)
df <- cbind(frac_es, frac_nes, frac_pvalue, frac_fdr)
rownames(df) <- colnames(eTerm$es)
rownames(df)[nrow(df)] <- "Mutation\nubiquity"
ind <- sort.int(frac_es, index.return=T)$ix
data <- df[ind,]

par(las=1) # make label text perpendicular to axis
par(mar=c(5,8,4,2)) # increase y-axis margin.
z <- data[,2]
barY <- barplot(z, xlab="Normalised enrichment score (NES)", horiz=TRUE, names.arg=rownames(data), cex.names=0.7, cex.lab=0.7, cex.axis=0.7, col="transparent")
#text(cex=0.8, y=barY, x=z-0.4, sprintf("adjusted p-value = %.1e",data[,3]), xpd=TRUE) 

for(i in 1:length(eTerm$nes)){
    which_sample <- i
    res <- dGSEAview(eTerm, which_sample=which_sample, top_num=5, sortBy="adjp", decreasing=FALSE, details=TRUE)
    wth <- 3600
    png(paste(colnames(eTerm$es)[which_sample], ".png", sep=""), width=wth, height=wth, res=wth*72/480)
    visGSEA(eTerm, which_sample=which_sample, which_term=rownames(res)[1], newpage=F)
    dev.off()
}

output <- dGSEAwrite(eTerm, which_content="gadjp", which_score="gadjp", filename="eTerm.txt")


####### GSEA for each sample
data <- exprs(eset)+0.001
eTerm <- dGSEA(data=data, identity="symbol", genome="Hs", ontology="Customised", customised.genesets=V(g)$name, weight=1, nperm=5000, RData.location="RData_Rd/data")
save(list=c("eTerm"), file="GSEA_each_sample.RData")

####
load(url("http://dnet.r-forge.r-project.org/data/Hs/org.Hs.string.RData"))
#load("RData_Rd/data/Hs/org.Hs.string.RData")
org.Hs.string
# restrict to those edges with high confidence (score>=700)
g_temp <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=700])
ind <- match(V(g_temp)$symbol, rownames(eset))
## for extracted graph
nodes_mapped <- V(g_temp)$name[!is.na(ind)]
g_temp <- dNetInduce(g_temp, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=T)
V(g_temp)$name <- V(g_temp)$symbol
setSeeds <- exprs(eset)[ind[!is.na(ind)],]
PTmatrix <- dRWR(g=g_temp, normalise=c("laplacian","row","column","none")[1], setSeeds=setSeeds, restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)
eTerm <- dGSEA(data=PTmatrix, identity="symbol", genome="Hs", ontology="Customised", customised.genesets=V(g)$name, weight=0, nperm=5000, RData.location="RData_Rd/data")
save(list=c("eTerm"), file="GSEA_each_sample_weight.RData")
####
sample_pvalue <- as.vector(eTerm$pvalue)
sample_fdr <- stats::p.adjust(sample_pvalue, method="BH")
sample_nes <- as.vector(eTerm$nes)
sample_es <- as.vector(eTerm$es)
names(sample_pvalue) <- names(sample_fdr) <- names(sample_nes) <- names(sample_es) <- colnames(eTerm$pvalue)
df <- data.frame(pData(eset), pvalue=sample_pvalue, fdr=sample_fdr, nes=sample_nes, es=sample_es)
save(list=c("eTerm","df"), file="GSEA_each_sample_weight.RData")
load("GSEA_each_sample_weight.RData")


#scp gg-pc4:./GSEA_each_sample.RData ./
load("GSEA_each_sample.RData")
sample_pvalue <- as.vector(eTerm$pvalue)
sample_fdr <- stats::p.adjust(sample_pvalue, method="BH")
sample_nes <- as.vector(eTerm$nes)
sample_es <- as.vector(eTerm$es)
names(sample_pvalue) <- names(sample_fdr) <- names(sample_nes) <- names(sample_es) <- colnames(eTerm$pvalue)
df <- data.frame(pData(eset[ind[!is.na(ind)],]), pvalue=sample_pvalue, fdr=sample_fdr, nes=sample_nes, es=sample_es)

## for beeswarm
library(beeswarm)
wth <- 3600
png("GSEA_samples_nes.png", width=wth*2.2, height=wth, res=wth*72/480)
par(las=2, mar=c(9,4,1,1)) # all axis labels horizontal
lbls <- paste(names(table(df$TCGA_tumor_type)), " (n=", table(df$TCGA_tumor_type), ")",  sep="")
beeswarm(nes ~ TCGA_tumor_type, data=df, corral= c("none", "gutter", "wrap", "random", "omit")[1], pch=16, pwcol=2+as.numeric(df$pvalue>5e-2), cex=0.6, horizontal=F, ylab="Normalized enrichment score (NES)", xlab="", labels=lbls, ylim=c(0.6,1.6))
legend("topright", legend = c("pvalue<0.05", "pvalue>0.05"), pch=16, col=2:3, box.col="transparent", border="transparent", cex=0.8)
boxplot(nes ~ TCGA_tumor_type, data=df, add=T, outline=F, horizontal=F, names=lbls, col=c("transparent"),border=par("fg"))
dev.off()

## for pie
tmp <- table(df$TCGA_tumor_type)
sig <- sapply(names(tmp), function(x) sum(df$pvalue[df$TCGA_tumor_type==x]<0.05))
nonsig <- sapply(names(tmp), function(x) sum(df$pvalue[df$TCGA_tumor_type==x]>=0.05))
frac_sig <- sig/(sig+nonsig)
tmp_df <- cbind(sig, nonsig, frac_sig)
#ind <- sort.int(frac_sig, index.return=T)$ix
#tmp_df <- tmp_df[ind,]
wth <- 3600
png("GSEA_samples_nes_pie.png", width=wth*3, height=wth, res=wth*72/480)
par(mfrow=c(1,11), mar=c(2,2,2,2))
for(i in 1:length(tmp)){
    pie(tmp_df[i,1:2], col=c("red","green"), labels=paste(c("",""), tmp_df[i,1:2], sep=""), main=rownames(tmp_df)[i], cex=1.5)
}
dev.off()




require(grDevices)
     pie(rep(1, 24), col = rainbow(24), radius = 0.9)
     
     pie.sales <- c(0.12, 0.3, 0.26, 0.16, 0.04, 0.12)
     names(pie.sales) <- c("Blueberry", "Cherry",
         "Apple", "Boston Cream", "Other", "Vanilla Cream")
     pie(pie.sales) # default colours
     pie(pie.sales, col = c("purple", "violetred1", "green3",
                            "cornsilk", "cyan", "white"))
     pie(pie.sales, col = gray(seq(0.4, 1.0, length = 6)))
     pie(pie.sales, density = 10, angle = 15 + 10 * 1:6)
     pie(pie.sales, clockwise = TRUE, main = "pie(*, clockwise = TRUE)")
     segments(0, 0, 0, 1, col = "red", lwd = 2)
     text(0, 1, "init.angle = 90", col = "red")
     
     n <- 200
     pie(rep(1, n), labels = "", col = rainbow(n), border = NA,
         main = "pie(*, labels=\"\", col=rainbow(n), border=NA,..")




# random walk with restart using denseness as seeds
PTmatrix <- dRWR(g=network, normalise=c("laplacian","row","column","none")[1], setSeeds=data.frame(frac_mutated, denseness), restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)

PTmatrix <- dRWR(g=network, normalise=c("laplacian","row","column","none")[1], setSeeds=frac_mutated, restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)


data <- t(PTmatrix)
n <- nrow(data)
res <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
    res[i, ] <- (data[i, ] %*% t(data))
}
rownames(res) <- rownames(data)
colnames(res) <- rownames(data)
visHeatmapAdv(res, Rowv=F, Colv=F, zlim=c(0.0001,0.0004))

adjmatrix <- res > quantile(res[lower.tri(res)], 0.75)

sg <- graph.adjacency(adjmatrix, mode="undirected", weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)
visNet(sg, edge.width=E(sg)$weight)



ind <- sample(1:nrow(frac_mutated))
seeds_random <- frac_mutated[ind,]
rownames(seeds_random) <- rownames(frac_mutated)
PT_random <- dRWR(g=network, normalise=c("laplacian","row","column","none")[1], setSeeds=seeds_random, restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)
data <- t(PT_random)
n <- nrow(data)
res_random <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
    res_random[i, ] <- (data[i, ] %*% t(data))
}
rownames(res_random) <- rownames(data)
colnames(res_random) <- rownames(data)
res_random[lower.tri(res_random)]

a <- cbind(res[lower.tri(res)], res_random[lower.tri(res_random)])


# GSEA using PS
eTerm <- dGSEA(data=PTmatrix, identity="symbol", genome="Hs", ontology="PS", sizeRange=c(10,100000), which_distance=NULL, sigTail=c("two-tails","one-tail")[2])
which_sample=12
res <- dGSEAview(eTerm, which_sample=1, top_num=5, sortBy="nES", decreasing=T, details=TRUE)
visGSEA(eTerm, which_sample=which_sample, which_term=rownames(res)[1])
output <- dGSEAwrite(eTerm, which_content="nES", which_score="nES", cutoff=0, filename="eTerm.txt",keep.significance=F)

# visualise using advanced heatmap
data <- output[,6:ncol(output)]
rownames(data) <- paste(output$setID,output$setSize,output$name, sep="_")
visHeatmapAdv(data, Rowv=F, Colv=T, colormap="darkgreen-lightgreen-lightpink-darkred", zlim=c(0,2), margins = c(7,14),cexRow=1,cexCol=1)


