\name{iassortativity}
\alias{iassortativity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Assortative Mixings }
\description{
The function calculates the assortative mixing of each vertex in an undirected, unweighted graph (or network)
}
\usage{
iassortativity(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ a symmetric adjacency matrix with elements zeros or ones representing an undirected and unweighted network }
}
\details{
Assortative mixing \eqn{knn_i}{knn_i} of vertex \eqn{i}{i} sums the degrees of its neighboring vertices. The summation is then normalized to the vertex's 
degree:
\deqn{knn_i={1\over k_i}\sum_j A_{ij}k_j}{knn_i=(1/k_i)(sum_j A_ij*k_j)}
where \eqn{k_j}{k_j} is the degree of vertex \eqn{j}{j} and \eqn{A_{ij}}{A_ij} the adjacency matrix of the network. A scatter plot of \eqn{knn}{knn} versus \eqn{k}{k} shows if hub (i.e. high-degree) vertices tend to attach to other hub vertices. 
}
\value{
  \item{knn }{ an array holding the assortative mixing of each vertex }
  \item{Knnk }{ the average knn of the vertices having the same degree k. An increasing Knnk with k indicates that the network is assortative characteristic of social networks. Most nonsocial networks such as technological and biological networks show disassortativity. }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Sun-Chong Wang }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{itransitivity}} }
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ misc }% 
