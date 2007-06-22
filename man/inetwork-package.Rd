\name{inetwork-package}
\alias{inetwork-package}
\alias{inetwork}
\docType{package}
\title{
Network Analysis and Plotting
}
\description{
\code{inetwork} implements a network partitioning algorithm to identify communities (or modules) in a network. The network plotting function then utilizes the identified community structure to position the vertices for plotting. The package also contains functions to calculate the assortativity, i.e. degree correlation, and transitivity, i.e. clustering coefficient, of a vertex.
}
\details{
\tabular{ll}{
Package: \tab inetwork\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-06-22\cr
License: \tab LGPL\cr
}
Identification of communities (or modules) in an undirected network (or graph) by \code{\link{icommunity}}, followed by organizing the communities within communities by \code{\link{ihierarchy}} and plotting of the network by \code{\link{inetplot}}. \code{\link{assortativity}} and \code{\link{transitivity}} in the package calculate a vertex's degree-correlation and clustering coefficient.
}
\author{
Sun-Chong Wang\cr 
Institute of Systems Biology and Bioinformatics, National Central University, Taiwan\cr
Epigenetics Laboratory, Centre for Addiction and Mental Health, Toronto Canada\cr
Maintainer: Sun-Chong Wang <sun-chong\_wang@camh.net>
}
\references{
M.E.J. Newman,
Modularity and community structure in networks,
\emph{Proc. Natl. Acad. Sci. U.S.A.} 103 (2006) 8577--8582.
}
\keyword{graphs}
\keyword{misc}
\seealso{
\code{\link{icommunity}}, \code{\link{ihierarchy}}, \code{\link{inetplot}}, \code{\link{assortativity}}, \code{\link{transitivity}}
}
\examples{
## load data
data(cashflow)

## network partitioning to find communities
ihierarchy(icommunity(cf5))

## network plotting showing the communities
inetplot(icommunity(cf5,partite=TRUE))
}
