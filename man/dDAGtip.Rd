\name{dDAGtip}
\alias{dDAGtip}
\title{Function to find the tip node(s) of a direct acyclic graph (DAG)}
\usage{
dDAGtip(g)
}
\arguments{
\item{g}{an object of class "igraph" or "graphNEL"}
}
\value{
\itemize{ \item{\code{tip}: the tip name (i.e. Term ID)} }
}
\description{
\code{dDAGtip} is supposed to find the tip node(s) of a
direct acyclic graph (DAG; an ontology). It return the name
(i.e Term ID) of the tip node(s).
}
\note{
none
}
\examples{
# 1) load GO Molelular Function as igraph object
load(url("http://dnet.r-forge.r-project.org/data/Obo/ig.GOMF.RData"))
g <- ig.GOMF

# 2) find tips
tips <- dDAGtip(g)
}
\seealso{
\code{\link{dDAGtip}}
}

