\name{dGSEA}
\alias{dGSEA}
\title{Function to conduct gene set enrichment analysis given the input
data and the ontology in query}
\usage{
dGSEA(data, identity = c("symbol", "entrez"),
genome = c("mm", "hs"),
ontology = c("GOBP", "GOMF", "GOCC", "MP", "DO", "PS"),
sizeRange = c(10, 1000), which_distance = NULL,
weight = 1, nperm = 100, fast = T,
sigTail = c("two-tails", "one-tail"), verbose = T)
}
\arguments{
\item{data}{a data frame or matrix of input data. It must
have row names, either Entrez Gene ID or Symbol}

\item{identity}{the type of gene identity (i.e. row names
of input data), either "symbol" for gene symbols (by
default) or "entrez" for Entrez Gene ID. The option
"symbol" is preferred as it is relatively stable from one
update to another; when gene symbols cannot be matched,
synonyms will be searched against}

\item{genome}{the genome identity. It can be either "mm"
for mouse genome or "hs" for human genome}

\item{ontology}{the ontology supported currently. It can
be "GOBP" for Gene Ontology Biological Process, "GOMF"
for Gene Ontology Molecular Function, "GOCC" for Gene
Ontology Cellular Component, "MP" for Mammalian
Phenotype, "HP" for Human Phenotype, "DO" for Disease
Ontology, and "PS" for phylostratific age}

\item{sizeRange}{the minimum and maximum size of members
of each gene set in consideration. By default, it sets to
a minimum of 10 but no more than 1000}

\item{which_distance}{which distance of terms in the
ontology is used to restrict terms in consideration. By
default, it sets to 'NULL' to consider all distances}

\item{weight}{type of score weigth. It can be "0" for
unweighted (an equivalent to Kolmogorov-Smirnov, only
considering the rank), "1" for weighted by input gene
score (by default), and "2" for over-weighted, and so on}

\item{nperm}{the number of random permutations. For each
permutation, gene-score associations will be permutated
so that permutation of gene-term associations is
realised}

\item{fast}{logical to indicate whether to fast calculate
expected results from permutated data. By default, it
sets to true}

\item{sigTail}{the tail used to calculate the statistical
significance. It can be either "two-tails" for the
significance based on two-tails or "one-tail" for the
significance based on one tail}

\item{verbose}{logical to indicate whether the messages
will be displayed in the screen. By default, it sets to
false for no display}
}
\value{
an object of class "eTerm", a list with following
components: \itemize{ \item{\code{set_info}: a matrix of
nSet X 4 containing gene set information, where nSet is
the number of gene set in consideration, and the 4
columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term
Name"), "namespace" and "distance"} \item{\code{gs}: a
list of gene sets, each storing gene members. Always,
gene sets are identified by "setID" and gene members
identified by "Entrez ID"} \item{\code{data}: a matrix of
nGene X nSample containing input data in consideration.
It is not always the same as the input data as only those
mappable are retained} \item{\code{es}: a matrix of nSet
X nSample containing enrichment score, where nSample is
the number of samples (i.e. the number of columns in
input data} \item{\code{nes}: a matrix of nSet X nSample
containing normalised enrichment score. It is the version
of enrichment score but after being normalised by gene
set size} \item{\code{pvalue}: a matrix of nSet X nSample
containing nominal p value} \item{\code{adjp}: a matrix
of nSet X nSample containing adjusted p value. It is the
p value but after being adjusted for multiple
comparisons} \item{\code{gadjp}: a matrix of nSet X
nSample containing globally adjusted p value in terms of
all samples} \item{\code{fdr}: a matrix of nSet X nSample
containing false discovery rate (FDR). It is the
estimated probability that the normalised enrichment
score represents a false positive finding}
\item{\code{qvalue}: a matrix of nSet X nSample
containing q value. It is the monotunically increasing
FDR} \item{\code{call}: the call that produced this
result} }
}
\description{
\code{dGSEA} is supposed to conduct gene set enrichment
analysis given the input data and the ontology in query.
It returns an object of class "eTerm".
}
\note{
The interpretation of returned components: \itemize{
\item{"es": enrichment score for the gene set is the
degree to which this gene set is overrepresented at the
top or bottom of the ranked list of genes in each column
of input data;} \item{"nes": normalised enrichment score
for the gene set is enrichment score that has already
normalised by gene set size. It is comparable across
analysed gene sets;} \item{"pvalue": nominal p value is
the statistical significance of the enrichment score. It
is not adjusted for multiple hypothesis testing, and thus
is of limited use in comparing gene sets;} \item{"adjp":
adjusted p value by Benjamini & Hochberg method. It is
comparable across gene sets;} \item{"gadjp": globally
adjusted p value by Benjamini & Hochberg method. Unlike
"adjp", it is adjusted in terms of all samples;}
\item{"fdr": false discovery rate is the estimated
probability that the normalised enrichment score
represents a false positive finding. Unlike "adjp" or
"gadjp" (also aliased as "fdr") that is derived from a
list of p values, this version of fdr is directly
calculate from the statistic (i.e. normalised enrichment
score);} \item{"qvalue": q value is the monotunically
increasing FDR so that the higher "nes", the lower
"qvalue".} }
}
\examples{
\dontrun{
load(url("http://dnet.r-forge.r-project.org/data/Hiratani_TableS1.RData"))
data <- RT[,1:2]
eTerm <- dGSEA(data, identity="symbol", genome="mm", ontology="MP",
which_distance=c(1,2))
}
}
\seealso{
\code{\link{dGSEAview}}, \code{\link{dGSEAwrite}},
\code{\link{visGSEA}}
}

