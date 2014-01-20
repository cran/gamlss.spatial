\name{gamlss.mrfa}
\alias{gamlss.mrfa}
\alias{mrfa}
\alias{mrf}
\alias{gamlss.mrf}


\title{Markov Random Field fitting within GAMLSS 
}
\description{
The functions \code{mrf()} and \code{mrfa()} can be used to fit Markov Random Field additive terms within GAMLSS. The two functions should produce identical results.
}
\usage{
mrf(x, precision = NULL, neighbour = NULL, polys = NULL, area = NULL, 
     adj.weight = 1000, control = mrf.control(...), ...)
mrfa(x, precision = NULL, neighbour = NULL, polys = NULL, area = NULL, 
     start = 10, df = NULL, adj.weight = 1000, ...)

gamlss.mrf(x, y, w, xeval = NULL, ...)
gamlss.mrfa(x, y, w, xeval = NULL, ...)
}

\arguments{
  \item{x}{a factor containing the areas}
  \item{precision}{the precision matrix if set}
 \item{neighbour}{an object containing the neighbour information for the area if set}
  \item{polys}{the polygon information if set}
  \item{area}{this argument is here to allow more areas than the levels of the factor \code{x}, see example below} 
  \item{start}{starting value for the smoothing parameter \code{lambda}
}
  \item{df}{degrees of freedom for fitting if required
}
  \item{adj.weight}{a value to adjust the iterative weight if necessary}
  \item{y}{working response variable}
  \item{w}{iterative weights}
  \item{xeval}{whether to predict or not}
   \item{control}{to be use for some of the srgument of \code{MRF()}.}
  \item{\dots}{for extra arguments}
}
\details{
The function \code{mrfa()} is to support  the function \code{MRFA()} within GAMLSS.
It is intended to be called within a GAMLSS formula. The function \code{gamlss.mrf()} is not intended to be used directly. It is calling the  function \code{MRFA()} within the GAMLSS fitting algorithm.
}
\value{
a fitted gamlss object
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rue and Held (2005) \emph{Gaussian markov random fields: theory and applications}, Chapman  & Hall, USA. 

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Fernanda De Bastiani, Mikis Stasinopoulos, Robert Rigby and Vlasios Voudouris. 

Maintainer: Fernanda <fernandadebastiani@gmail.com>
}


\seealso{
 \code{\link{MRF}}, \code{\link{MRFA}}
}
\examples{
library(gamlss)
library(mgcv)
data(columb)
data(columb.polys)
vizinhos=polys2nb(columb.polys)
precisionC <- nb2prec(vizinhos,x=columb$district)
# MRFA
 m1<- gamlss(crime~ mrf(district, polys=columb.polys), data=columb)
 m2<- gamlss(crime~ mrfa(district, polys=columb.polys), data=columb)
AIC(m1,m2, k=0)
draw.polys(columb.polys, getSmo(m2), scheme="topo")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{spatial}% __ONLY ONE__ keyword per line
