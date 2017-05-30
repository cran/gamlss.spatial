#This function was create to plot the interaction between factors

int.term<-function(object, 
                   xvar, 
                   position, 
                   fac,
                   which.lev=NULL,
                   xlabel=NULL, 
                   ylabel=NULL,  
                   scheme = c( "shaded", "lines"),
                   col.term = "darkred",  
                   lwd.term = 1.5,   
                   lty.se = 2 ,
                   lwd.se = 1,
                   col.se = "orange" ,  
                   col.shaded = "gray",
                   factor.plots = FALSE,
                   main="factor",...)
{
  xlabel<- if(is.null(xlabel)) "x" else xlabel
  ylabel<- if(is.null(ylabel)) "f(x)" else ylabel
      t1<-lpred(object, type="terms", se=TRUE)
     fv <- as.vector(t1$fit[,position])
    var <- t1$se.fit[,position]
  upper <- fv+2*var
  lower <- fv -2*var
    ran <- range(upper, lower)*c(.95, 1.05)
  if (!factor.plots)
  {
    plot(fv[order(xvar)]~xvar[order(xvar)], type="n", ylim=ran,
         ylab=ylabel, xlab=xlabel , main="interaction",  col = col.term, 
         lwd = lwd.term)   
  }
  nlevs <- nlevels(fac)
  if(!is.null(which.lev))
   {    
    if(which.lev%in%levels(fac)) lev.fac<-which.lev
    else stop("The level provide is not in the level of the factor")
   }
  else
   {
    lev.fac<- levels(fac) 
   }
  for (i in lev.fac)
  {
    if (factor.plots)
    {
      plot(fv[order(xvar)]~xvar[order(xvar)], type="n", ylim=ran,
           ylab=ylabel, xlab=xlabel , main=main,  col = col.term, 
           lwd = lwd.term)   
    }
       fvi <- fv[fac==i]
        xi <-  xvar[fac==i]
        ox <- order(xi)
    upperi <- upper[fac==i]
    loweri <- lower[fac==i]
    if (scheme=="lines")
    {
      lines(xi[ox],  upperi[ox] , lty = lty.se, lwd = lwd.se, col = col.se)
      lines(xi[ox],  loweri[ox],  lty = lty.se, lwd = lwd.se, col = col.se)
      lines(xi[ox],  fvi[ox], col = col.term, lwd = lwd.term)  
    } else
    {
      x1 <- xi[ox] 
      xx <- c(x1,rev(x1))
      yy <- c(loweri[ox] , upperi[rev(ox)])
      polygon(xx, yy, col = col.shaded, border = col.shaded)
      lines(xi[ox], fvi[ox], col = col.term, lwd = lwd.term)  
    } 
  }
}   


