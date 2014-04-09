# This is a demo for TCGA mutational profile dataset from Kandoth et al
# 
# This dataset is available from TCGA (see <a href="http://www.ncbi.nlm.nih.gov/pubmed/24132290" target="24132290">http://www.ncbi.nlm.nih.gov/pubmed/24132290</a>), containing somatic mutational profiles for 3130 cancer samples (i.e. only those having clinical data). These cancer samples belong to one of 12 major cancer types, including breast adenocarcinoma (BRCA), lung adenocarcinoma (LUAD), lung squamous cell carcinoma (LUSC), uterine corpus endometrial carcinoma (UCEC), glioblastoma multiforme (GBM), head and neck squamous cell carcinoma (HNSC), colon and rectal carcinoma (COAD/READ), bladder urothelial carcinoma (BLCA), kidney renal clear cell carcinoma (KIRC), ovarian serous carcinoma (OV) and acute myeloid leukaemia (LAML). For each patient sample, somatic mutations are represented as a profile of binary (1, 0) states on genes, where '1' indicates a gene for which mutation has occurred in the tumor relative to germ line. The dataset is stored in "TCGA_mutations" (as an 'ExpressionSet' object), including:
## assayData: exprs(TCGA_mutations), a matrix of 19477 genes X 3130 samples;
## phenoData: pData(TCGA_mutations), variables describing sample phenotypes (i.e. columns in assayData), including information about samples: "Age", "Days_to_death", "Days_to_last_followup", "Vita_status", "Gender", "Date_diagnosis", "TCGA_tumor_type", "Tumor_stage", "Tumor_grade"
## featureData: fData(TCGA_mutations), variables describing features (i.e. rows in assayData), including information about features/genes: "EntrezID" for gene EntrezID, "Symbol" for gene symbol, "Desc" for gene description, "Synonyms" for gene symbol alias
###############################################################################
library(dnet)

# Load or install packages specifically used in this demo
list.pkg <- c("affy")
for(pkg in list.pkg){
    if(!require(pkg, character.only=T)){
        install.packages(pkg,repos="http://www.stats.bris.ac.uk/R",dependencies=TRUE)
        lapply(pkg, library, character.only=T)
    }
}

# load an "ExpressionSet" object
load(url("http://dnet.r-forge.r-project.org/data/Datasets/TCGA_mutations.RData"))
## remove those genes without mutations in any samples
flag <- apply(exprs(TCGA_mutations)!=0, 1, sum)
eset <- TCGA_mutations[which(flag!=0),]
# extract information about the first 5 samples
pData(eset)[1:5,]
# extract information about the first 5 features
fData(eset)[1:5,]
# number of samples for each cancer type
table(pData(eset)$TCGA_tumor_type)
tumor_type <- sort(unique(pData(eset)$TCGA_tumor_type))

# Entrez Gene information for the human
load(url("http://dnet.r-forge.r-project.org/data/Hs/org.Hs.eg.RData"))
gene_info <- org.Hs.eg$gene_info
gene_info[1:2,]

# An igraph object that contains a functional protein association network in human. The network is extracted from the STRING database (version 9.1). Only those associations with medium confidence (score>=400) are retained.
load(url("http://dnet.r-forge.r-project.org/data/Hs/org.Hs.string.RData"))
org.Hs.string
# restrict to those edges with high confidence (score>=700)
g <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=700])
g

# extract network that only contains genes in eset
ind <- match(V(g)$symbol, rownames(eset))
## for extracted expression
esetGene <- eset[ind[!is.na(ind)],]
esetGene
## for extracted graph
nodes_mapped <- V(g)$name[!is.na(ind)]
network <- dNetInduce(g, nodes_query=nodes_mapped, knn=0, remove.loops=F, largest.comp=F)
V(network)$name <- V(network)$symbol
#E(network)$weight <- E(network)$combined_score
network

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
# derive the "mutation denseness" of genes
denseness <- 1- sparseness
hist(denseness,100)

# random walk with restart using denseness as seeds
PTmatrix <- dRWR(g=network, normalise=c("laplacian","row","column","none")[1], setSeeds=data.frame(frac_mutated, denseness), restart=0.5, normalise.affinity.matrix=c("none","quantile")[1], verbose=T)

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
