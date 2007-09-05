\name{itransitivity}
\alias{itransitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Clustering Coefficient }
\description{
The clustering coefficient of each vertex of an undirected and unweighted network is calculated given the network's adjacency matrix
}
\usage{
itransitivity(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ a symmetric adjacency matrix with elements zeros or ones representing an undirected and unweighted network }
}
\details{
The clustering coefficient or transitivity \eqn{c_i}{c_i} of vertex \eqn{i}{i} is defined as
\deqn{c_i={1\over k_i(k_i-1)}\sum_{j,h}A_{ij}A_{ih}A_{jh}}{c_i = (1/(k_i(k_i-1))) sum_jh (A_ij*A_ih*A_jh)}
where \eqn{k_i}{k_i} is the degree of vertex \eqn{i}{i} and \eqn{A_{ij}}{A_ij} the adjacency matrix of the network. The clustering coefficient measures the connectivity of the neighboring vertices of a vertex. A decreasing \eqn{c}{c} with \eqn{k}{k} indicates that the network is hierarchical as observed in many nonrandom networks.
}
\value{
  \item{c }{an array holding the clustering coefficients of the vertices }
  \item{Ck }{the averaged c's for the vertices having the same degree k}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Sun-Chong Wang }
\note{ 
For large networks consisting of 700 vertices, the calculation can take ten
minutes. The computation time/resource scales as \eqn{O(n^3)}{O(n^3)}.
}
\seealso{ \code{\link{iassortativity}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ misc }% 
