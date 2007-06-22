\name{icommunity}
\alias{icommunity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identification of Communities or Modules in a Network }
\description{
Given a network represented by an adjacency matrix, the function returns the communities in the network. Community detection is an optimization problem which is achieved by the leading eigenvector of the so-called modularity matrix.
}
\usage{
icommunity(A, nodes = NULL, partite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ an adjacency matrix with elements 0's or 1's representing the network }
  \item{nodes}{ a vector of strings for the labels of the vertices }
  \item{partite}{ a logical variable indicating whether within- or between-community edge densities are maximized in network partitioning for community identification }
}
\details{
Networks or graphs are represented by the so-called adjacency matrices whose elements are either zeros or ones. A one (zero) at element ij indicates that vertices i and j in the network are (not) connected by an edge. The edges are undirected. The input adjacency matrix should therefore be square and symmetric. If the character strings for the vertex labels are not provided, i.e. \code{nodes=NULL}, "1", "2", "3", ..., etc will be used as the labels for the vertices. Highly interacting vertices form a community or module. The algorithm in the implementation seeks to partition the network such that the densities of the edges within the sought communities are larger than the average density. This mode of community detection is the default \code{partite=FALSE}. When \code{partite} is set to \code{TRUE}, the returned communities have larger between- than within-community edge densities.   
}
\value{
  \item{A }{A copy of the input adjacency matrix. It is stored here for subsequent plotting of the network by \code{\link{inetplot}} in the package}
  \item{Q }{A measure of modularity. In the network partitioning algorithm, the original network is divided into two subnetworks. The modularity of the network gains as a result of the division. Each of the subnetworks is recursively subject to division until the modularity of the network gains no more. The gains in the modularity in each division are used in the hierarchical plotting of the network by \code{\link{ihierarchy}} in the package}
  \item{sizes }{An integer array holding the sizes of the identified communities }
  \item{indivisibles }{Vertex labels of the identified communities }
  \item{ringleaders }{The ringleadership of the vertices in their communities; the larger the ringleadership, the larger drop in the modularity if the vertex is moved out of her designated community to any other community } 
  \item{arrows }{The horizontal and vertical bars used for the hierarchical plotting by \code{\link{ihierarchy}} }
  \item{xcoordinates}{The x-coordinates of the vertex labels for \code{\link{ihierarchy}} }
  \item{ycoordinates}{The y-coordinates of the vertex labels for \code{\link{ihierarchy}} }
  \item{indices }{Matrix indices of the identified community members }
  \item{singlets }{The vertices that are not connected to others in the network }
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
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{ihierarchy}}, \code{\link{inetplot}} }
\examples{
## load example networks 
data(cashflow)

## community identification
cluster3 <- icommunity(cf3,labelcf3,partite=FALSE)
## network plotting
inetplot(cluster3,shaft=2,circle=FALSE,theta=33,points=FALSE)

## anti-community identification
partite3 <- icommunity(cf3,labelcf3,partite=TRUE)
## plot the identified anti-communities
inetplot(partite3,shaft=10,circle=FALSE,theta=33,points=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ misc }% 
