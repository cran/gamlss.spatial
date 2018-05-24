\name{gamlss.spatial-package}
\alias{gamlss.spatial-package}
\alias{gamlss.spatial}
\docType{package}
\title{
\packageTitle{gamlss.spatial}
}
\description{
\packageDescription{gamlss.spatial}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{gamlss.spatial}
\packageIndices{gamlss.spatial}
}
\author{
\packageAuthor{gamlss.spatial}

Maintainer: \packageMaintainer{gamlss.spatial}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.
Rue and Held (2005) \emph{Gaussian markov random fields: theory and applications}, Chapman & Hall, USA. 

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\keyword{ package }
\examples{
library(mgcv)
data(columb)
data(columb.polys)
m1 <- MRFA(columb$crime, columb$district, polys=columb.polys)
draw.polys(columb.polys, m1) 
}
