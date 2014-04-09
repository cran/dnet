\name{dDAGroot}
\alias{dDAGroot}
\title{Function to find the root node of a direct acyclic graph (DAG)}
\usage{
dDAGroot(g)
}
\arguments{
\item{g}{an object of class "igraph" or "graphNEL"}
}
\value{
\itemize{ \item{\code{root}: the root name (i.e. Term ID)}
}
}
\description{
\code{dDAGroot} is supposed to find the root node of a
direct acyclic graph (DAG; an ontology). It return the name
(i.e Term ID) of the root node.
}
\note{
none
}
\examples{
# 1) load GO Molelular Function as igraph object
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOMF.RData"))
g <- ig.GOMF

# 2) find the root
root <- dDAGroot(g)
}
\seealso{
\code{\link{dDAGroot}}
}

