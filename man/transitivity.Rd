\name{transitivity}
\alias{transitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Clustering Coefficient }
\description{
The clustering coefficient of each vertex of an undirected and unweighted network is calculated given the network's adjacency matrix
}
\usage{
transitivity(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ A symmetric adjacency matrix with elements zeros or ones representing an undirected and unweighted network }
}
\details{
The clustering coefficient or transitivity \eqn{c_i}{c_i} of vertex \eqn{i}{i} is defined as
\deqn{c_i={1\over k_i(k_i-1)}\sum_{j,h}a_{ij}a_{ih}a_{jh}}{c_i = (1/(k_i(k_i-1))) sum_jh (a_ij*a_ih*a_jh)}
where \eqn{k_i}{k_i} is the degree of vertex \eqn{i}{i} and \eqn{a_{ij}}{a_ij} the adjacency matrix of the network. The clustering coefficient measures the connectivity of the neighboring vertices of a vertex. A decreasing \eqn{c}{c} with \eqn{k}{k} indicates that the network is hierarchical as observed in many nonrandom networks.
}
\value{
  \item{c }{An array holding the clustering coefficients of the vertices }
  \item{Ck }{The averaged c's for the vertices having the same degree k}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Sun-Chong Wang }
\note{ 
For large networks consisting of 700 vertices, the calculation can take ten
minutes.
}
\seealso{ \code{\link{assortativity}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ misc }% 
