\name{ihierarchy}
\alias{ihierarchy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot the Communal or Modular Hierarchy of a Network }
\description{
Given the communities identified by \code{icommunity}, the function displays the hierarchical structure of the communities within communities in the network.
}
\usage{
ihierarchy(clusters, fan = FALSE, spread = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{ the object returned by \code{\link{icommunity}} }
  \item{fan}{ a logical variable indicating whether the plot is arched or not }
  \item{spread}{ the distance between the communities along the arc }
}
\details{
In the hierarchy, the length of a vertical bar indicates the gain in modularity resulting from the division underneath the shoulder of the bar. (Specifically, length=\eqn{(\Delta Q)^(-1/3)}{(Delta Q)^(-1/3)} with \eqn{\Delta Q}{Delta Q} in \code{clusters$Q}.) One of the applications is to best identify the (number of) communities together with information other than network connectivity. When the size of the network and/or the number of identified communities become large, display of the communities/vertices can be along concentric arcs or circles by setting \code{fan=TRUE}. The distances between the communities on the arcs can be further tuned by \code{spread}.
}
%\value{
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
%\references{ ~put references to the literature/web site here ~ }
\author{ Sun-Chong Wang }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{inetplot}}, \code{\link{icommunity}} }
\examples{
## load example networks
data(icashflow)

## hierarchical community structure in the network
ihierarchy(icommunity(cf9,labelcf9))

## display the communities along arcs
ihierarchy(icommunity(cf9,labelcf9),fan=TRUE,spread=0.5)

## the separations are increased
ihierarchy(icommunity(cf9,labelcf9),fan=TRUE,spread=1.7)
}
\keyword{ graphs }
\keyword{ misc }% 
