\name{dDAGreverse}
\alias{dDAGreverse}
\title{Function to reverse the edge direction of a direct acyclic graph (DAG)}
\usage{
dDAGreverse(g)
}
\arguments{
\item{g}{an object of class "igraph" or "graphNEL"}
}
\value{
\itemize{ \item{\code{gr}: a graph being reversed, an
object of class "igraph" or "graphNEL"} }
}
\description{
\code{dDAGreverse} is supposed to reverse the edge
direction of a direct acyclic graph (DAG; an ontology). The
return graph remains all attributes associated on nodes and
edges.
}
\note{
none
}
\examples{
# 1) load GO Molelular Function as igraph object
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOMF.RData"))
g <- ig.GOMF

# 2) the graph with reverse edge direction
gr <- dDAGreverse(g)
}
\seealso{
\code{\link{dDAGreverse}}
}

