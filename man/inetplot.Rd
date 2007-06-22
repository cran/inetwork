\name{inetplot}
\alias{inetplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Network or Graph Plotting }
\description{
Given the communities identified by \code{icommunity}, the plotting algorithm places the communities with decreasing community sizes along a spiral, in an effort to show the hierarchical structure of the network. The vertices belonging to a community are arranged on a circle at the designated position on the spiral. The radii of the circles are proportional to the community sizes.
}
\usage{
inetplot(clusters, theta = 30, shaft = 1, circle = TRUE, singlets = FALSE, 
         labels = TRUE, edges = TRUE, points = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{ An object returned from the network partitioning function \code{\link{icommunity}} }
  \item{theta}{ Angle in degree between neighboring communities on the spiral }
  \item{shaft}{ A scaling factor for the arm of the spiral }
  \item{circle}{ A logical to indicate whether to position the vertices in a community on a circle (default) or spiral }
  \item{singlets}{ A logical to turn off (default) displaying isolated vertices }
  \item{labels}{ A logical to show vertex labels (default) or not }
  \item{edges}{ A logical to draw an edge between connected vertices (by default) or not }
  \item{points}{ A logical to draw dots on vertices (by default) or not }
}
\details{
The Cartesian coordinates of community \eqn{i}{i} on the spiral are given by
\deqn{x_i=shaft\times (i-1)\times \cos(i\theta)}{xi = shaft * (i-1) * cos(i*theta)}
\deqn{y_i=shaft\times (i-1)\times \sin(i\theta)}{yi = shaft * (i-1) * sin(i*theta)} 
where \eqn{i}{i} is the \eqn{i}{i}th largest community in the network. The pitch angle \eqn{\theta}{theta} is set by the argument \code{theta}. If \code{circle=FALSE}, the vertices of a community are arranged along a spiral originating from the designated position on the backbone spiral. The \code{theta} and \code{shaft} apply to the backbone as well as to the component spirals. 
}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
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
\seealso{ \code{\link{icommunity}}, \code{\link{ihierarchy}} }
\examples{
## load example networks
data(cashflow)

## identify the communities within the network
partite7 <- icommunity(cf7,labelcf7,partite=TRUE)

## arrange the communities on a spiral
inetplot(partite7,shaft=10,circle=TRUE,labels=FALSE)

## arrange the community members along component spirals
inetplot(partite7,shaft=50,circle=FALSE,labels=TRUE,points=FALSE)

## show the isolated vertices
inetplot(partite7,shaft=10,points=FALSE,singlets=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ misc }% 
