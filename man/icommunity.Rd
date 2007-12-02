\name{icommunity}
\alias{icommunity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identification of Communities or Modules in a Network }
\description{
Given a network represented by an adjacency matrix, the function returns the communities or modules in the network. Community detection is an optimization problem which is achieved by the leading eigenvector of the so-called modularity matrix.
}
\usage{
icommunity(A, nodes = NULL, partite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ a symmetric adjacency matrix with elements 0's or 1's representing the network }
  \item{nodes}{ a vector of character strings for the labels of the vertices }
  \item{partite}{ a logical variable indicating whether within- or between-community edge densities are maximized in network partitioning for community identification }
}
\details{
Networks or graphs are represented by the so-called adjacency matrices whose elements are either zeros or ones. A one (zero) at element \eqn{ij}{ij} indicates that vertices \eqn{i}{i} and \eqn{j}{j} in the network are (not) connected by an edge. The edges are undirected. The input adjacency matrix should therefore be symmetric. \code{rownames(A)} or \code{colnames(A)}, if not \code{NULL}, will be used as vertex labels. Or provide your labels via \code{nodes}. Otherwise "1", "2", "3", ..., etc will be used as labels for the vertices. Highly interacting vertices form a community or module. The algorithm in the implementation seeks to partition the network such that the densities of the edges within the sought communities are larger than the average density. This mode of community detection is the default \code{partite=FALSE}. When \code{partite} is set to \code{TRUE}, the returned communities have larger between- than within-community edge densities.   
}
\value{
  \item{A }{a copy of the input adjacency matrix. It is stored here for subsequent plotting of the network by \code{\link{inetplot}} in the package}
  \item{Q }{a measure of modularity. In the network partitioning algorithm, the original network is divided into two subnetworks. The modularity of the network gains as a result of the division. Each of the subnetworks is recursively subject to division until the modularity of the network gains no more. The gains in the modularity in each division are stored here and used in the hierarchical plotting of the network by \code{\link{ihierarchy}} in the package }
  \item{sizes }{an integer array holding the sizes of the identified communities }
  \item{indivisibles }{vertex labels of the identified communities }
  \item{ringleaders }{the ringleadership of the vertices in their communities; the larger the ringleadership, the larger drop in the modularity if the vertex is moved out of her designated community to any other community } 
  \item{arrows }{the horizontal and vertical bars used for the hierarchical plotting by \code{\link{ihierarchy}} }
  \item{xcoordinates}{the x-coordinates of the vertex labels for \code{\link{ihierarchy}} }
  \item{ycoordinates}{the y-coordinates of the vertex labels for \code{\link{ihierarchy}} }
  \item{indices }{matrix indices of the singlets and identified community members }
  \item{singlets }{the vertices that are not connected to others in the network }
}
\references{ 
\item{}{
M.E.J. Newman,
Modularity and community structure in networks,
\emph{Proc. Natl. Acad. Sci. U.S.A.} 103 (2006) 8577--8582}

\item{}{M.E.J. Newman,
Finding community structure in networks using the eigenvectors of matrices,
\emph{Phys. Rev. E} 74 (2006) 036104-1--19}

}
\author{ Sun-Chong Wang }
\note{ 
Computing time/resource depends on the size and complexity of the network,
and scales as \eqn{O(n^2\log(n))}{O(n^2log(n))}.
}
\seealso{ \code{\link{ihierarchy}}, \code{\link{inetplot}} }
\examples{
## load example networks 
data(icashflow)

## community identification
cluster3 <- icommunity(cf3,labelcf3,partite=FALSE)
## network plotting
inetplot(cluster3,shaft=2,circle=FALSE,theta=33,points=FALSE)

## poly-partite partition
partite3 <- icommunity(cf3,labelcf3,partite=TRUE)
## plot the multi-partite graph
inetplot(partite3,shaft=10,circle=FALSE,theta=33,points=FALSE)

## to extract the community members:
partite3_community <- numeric()
for (i in 1:length(partite3$sizes)) { # loop over no. of communities
    partite3_community <- c(partite3_community,rep(i,partite3$sizes[i]))
}
for (i in 1:length(partite3$sizes)) { 
    print(c(paste("community",as.character(i),"="),partite3$indivisibles[which(partite3_community == i)]))
}
## the ringleadership scores of the nodes are extracted in the same fashion
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ misc }% 
