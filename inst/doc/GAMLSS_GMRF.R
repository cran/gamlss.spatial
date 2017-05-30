## ----message=FALSE-------------------------------------------------------
library(gamlss.spatial) 
library(mgcv)
# bring the data
data(columb)
names(columb)
# getting the polygons file 
data(columb.polys)
head(columb.polys,1)

## ------------------------------------------------------------------------
# getting the neighbours object from the polygons object
vizinhos <- polys2nb(columb.polys)
vizinhos[[1]]["0"]

## ------------------------------------------------------------------------
# getting the precision matrix from the neighbours object
precisionC <- nb2prec(vizinhos,x=columb$district)
precisionC[1:10, 1:20]

## ------------------------------------------------------------------------
# fit using the  polygone information
# MRFA alternaing 
system.time(m11<-MRFA(columb$crime, columb$district, polys=columb.polys))
# MRF Q-function
system.time(m21<-MRF(columb$crime, columb$district, polys=columb.polys))
# fit using the neighbour information
# MRFA alternaing
system.time(m12<-MRFA(columb$crime, columb$district,  neighbour=vizinhos))
# MRF Q-function
system.time(m22<-MRF(columb$crime, columb$district, neighbour=vizinhos))
# fit using the percision matrix
# MRFA alternaing
system.time(m13<-MRFA(columb$crime, columb$district, precision=precisionC))
# MRF Q-function
system.time(m23<-MRF(columb$crime, columb$district, precision=precisionC))
AIC(m11, m21, m12, m22, m13, m23, k=0)

## ------------------------------------------------------------------------
summary(m11)
summary(m21)

## ----eval=FALSE, echo=FALSE----------------------------------------------
#  plot(m11)
#  wp(m11)
#  # both look terrible
#  # model for sigma?

## ----ColumbData, fig.show='hide', fig.asp=1------------------------------
cr <- columb$crime
names(cr) <- as.character(columb$district)
draw.polys(columb.polys, cr,   scheme="topo",swapcolors=TRUE)
title("(a) crime")
draw.polys(columb.polys, m11, scheme="topo",swapcolors=TRUE)
title("(b) smooth values")

## ------------------------------------------------------------------------
# the dimension of the original data
dim(columb)
# removing one area (district ''4'') from  the data 
columb2 <- columb[-5,]
# drop unused level from a factor
columb2$district <-droplevels(columb2$district)
dim(columb2)
nlevels(columb2$district)
# fitting the reduced data 
# using  polys
r1<-MRF(columb2$crime, columb2$district, polys=columb.polys, 
        area=columb$district)
# using neighbours 
r2<-MRF(columb2$crime, columb2$district,  neighbour=vizinhos, 
        area=columb$district)
# using the old precision
r3<-MRF(columb2$crime, columb2$district, precision=precisionC,
        area=columb$district)
# creating new precision matrix
precisionC2 <- nb2prec(vizinhos, x=columb2$district, 
                       area=columb$district)
dim(precisionC2)
# fitting  using the new 49 x 49 precision
r4<-MRF(columb2$crime, columb2$district, precision=precisionC2,
         area=columb$district)
# checking the results
AIC(r1,r2,r3,r4, k=0)

## ----ColumbDataMis, fig.show='hide', fig.asp=1---------------------------
draw.polys(columb.polys, fitted(r1), scheme="heat", swapcolors=TRUE )
title("(a)")
draw.polys(columb.polys, r1, scheme="heat",swapcolors=TRUE  )
title("(b)")

## ----warning=FALSE-------------------------------------------------------
# fit the model
g1 <- gamlss(crime~gmrf(district,precision=precisionC), data=columb)
# comparing the fitted values
head(cbind(fitted(m21),fitted(g1)))
tail(cbind(fitted(m21),fitted(g1)))
#  the log-sigma coefficients
coef(m21)
coef(getSmo(g1))
# get sigma_b
m21$sigb
getSmo(g1)$sigb
# get sigma_e
m21$sige
fitted(g1,"sigma")[1]
# comparing the deviances
deviance(g1)
deviance(m21)
# get degrees of freedom for mu
g1$mu.df

## ----message=FALSE-------------------------------------------------------
library(gamlss.add)
g2 <- gamlss(crime~ga(~s(x,y)), data=columb)
AIC(g1,g2)

## ----ColumbDataTwo, fig.show='hide', fig.asp=1---------------------------
names(g1$mu.fv) <- names(g2$mu.fv) <- as.character(columb$district)
draw.polys(columb.polys, fitted(g2), scheme="terrain",swapcolors=TRUE  )
title("(a) 2-d smoothing")
draw.polys(columb.polys, fitted(g1), scheme="terrain", swapcolors=TRUE )
title("(b) IAR")

## ----eval=FALSE, echo=FALSE----------------------------------------------
#  plot(fitted(g1)~fitted(g2))

## ----message=FALSE-------------------------------------------------------
columb <- transform(columb, logos=log(open.space+1))
e0 <- gamlss(crime~1, data=columb)
# linear
e1 <- stepGAIC(e0, scope=list(lower=~1, upper=~area+income+home.value+logos))
formula(e1)

## ----ColumbDataLinear, fig.show='hide', fig.asp=1------------------------
term.plot(e1)

## ----eval=FALSE----------------------------------------------------------
#  ge1<- gamlss(crime~income+home.value+gmrf(district,precision=precisionC),
#               data=columb)

## ----message=FALSE, warning=FALSE, cache=TRUE----------------------------
e2 <- stepGAIC(e0, scope=list(lower=~1, upper=~pb(area)+pb(income)+
                      pb(home.value)+pb(logos)))
formula(e2)

## ----ColumbDataAdditive, fig.show='hide', fig.asp=1----------------------
term.plot(e2)

## ----ColumbDataNeural, fig.show='hide', fig.asp=1, warning=FALSE---------
e3 <- gamlss(crime~nn(~income+home.value+logos+area, decay=0.1), data=columb)
term.plot(e3)

## ----ColumbDataTree, fig.show='hide', fig.asp=1, warning=FALSE-----------
e4 <- gamlss(crime~tr(~income+home.value+logos+area), data=columb)
term.plot(e4)

## ------------------------------------------------------------------------
AIC(g1,g2, e1, e2, e3, e4)
AIC(g1,g2, e1, e2, e3, e4, k=log(49))

## ----ColumbDataFour, fig.show='hide', fig.asp=1--------------------------
names(e4$mu.fv) <- names(e3$mu.fv) <- names(e2$mu.fv) <-
  as.character(columb$district)
draw.polys(columb.polys, fitted(g1), scheme="terrain", swapcolors=TRUE )
title("(a) IAR")
draw.polys(columb.polys, fitted(e2), scheme="terrain",swapcolors=TRUE  )
title("(b) Additive Smooth")
draw.polys(columb.polys, fitted(e3), scheme="terrain",swapcolors=TRUE  )
title("(c) neural network")
draw.polys(columb.polys, fitted(e2), scheme="terrain",swapcolors=TRUE  )
title("(d) decision tree")

## ----eval=FALSE, echo=FALSE----------------------------------------------
#  # additive
#  e4 <- gamlss(crime~pb(income)+pb(home.value)+pb(log(open.space+1)), data=columb)
#  # nn surface
#  
#  k4 <- gamlss(crime~pbz(income)+pbz(home.value)+pbz(log(open.space+1))+pbz(area), data=columb)
#  k6 <- gamlss(crime~nn(~income+home.value+open.space+area, decay=0.1), data=columb)
#  AIC(g1,k2, k3, k4)
#  AIC(g1,k2, k3, k4, k=log(49))
#  
#  # MRF fit using the  polygone information
#  #tem.time(m1<-MRFA(columb$crime, columb$district, polys=columb.polys))
#  
#  # library(gamlss.add)
#  # k3 <- gamlss(crime~ga(~s(x,y)), data=columb)
#  # fv <- fitted(k3)
#  # names(fv) <- names(columb.polys)
#  # draw.polys(columb.polys, fv)
#  # k4 <- gamlss(crime~ga(~s(x,y)), data=columb, family=)
#  # plot(k4)
#  # fv <- fitted(k3)
#  # names(fv) <- names(columb.polys)
#  # draw.polys(columb.polys, fv)
#  #
#  # k1 <- gamlss(crime~gmrf(district,precision=precisionC), data=columb, family=GU, bf.cyc=1, n.cyc=5, glm.t)
#  # #plot(k0)
#  # #draw.polys(columb.polys,fitted(k0), scheme="topo",swapcolors=TRUE)
#  # #draw.polys(columb.polys,resid(k0))
#  #
#  # k1 <- gamlss(crime~income+home.value+open.space+gmrf(district, precision=precisionC, method="A"), data=columb)
#  # #k1 <- gamlss(crime~income+home.value+open.space+gmrf(district, precision=precisionC, method="Q"), data=columb)
#  # summary(k1)
#  # draw.polys(columb.polys,getSmo(k1))

## ----warning=FALSE, message=FALSE----------------------------------------
library(gamlss.spatial)
data(rent99)
data(rent99.polys)
rent99$cheating<-relevel(rent99$cheating,"1")
# creating new variables for interactions
# heating and years interaction
cy<-(as.numeric(rent99$cheating)-1)*rent99$yearc
# kitchen and years interation
ky<-(as.numeric(rent99$kitchen)-1)*rent99$yearc
# kitchen and area interation
ka<-(as.numeric(rent99$kitchen)-1)*rent99$area
# heating has its relevant level changed from 0 to 1 
heating<-relevel(rent99$cheating,"1")
rent99 <- transform(rent99,heating=heating, cy=cy, ky=ky, ka=ka) 

## ----Rent99_data_resp, fig.show='hide', fig.asp=1------------------------
hist(rent99$rent,ylab="f(y)",main="Histogram of rent", xlab="rent")
boxplot(rent99$rent)

## ----Rent99_data_xvar, fig.show='hide', fig.asp=1------------------------
plot(rent99$rent~rent99$area, data=rent, col=gray(0.7), 
     pch = 15, cex = 0.5, xlab="area", ylab="rent")
plot(rent99$rent~rent99$yearc, data=rent, col=gray(0.7), 
     pch = 15, cex = 0.5, xlab="yearc", ylab="rent")
plot(rent99$rent~rent99$location, data=rent, col=gray(0.7), 
     pch = 15, cex = 0.5, xlab="location", ylab="rent")
plot(rent99$rent~rent99$bath, data=rent, col=gray(0.7), 
     pch = 15, cex = 0.5, xlab="bath", ylab="rent")
plot(rent99$rent~rent99$kitchen, data=rent, col=gray(0.7), 
     pch = 15, cex = 0.5, xlab="kitchen", ylab="rent")
plot(rent99$rent~rent99$cheating, data=rent, col=gray(0.7), 
     pch = 15, cex = 0.5, xlab="cheating", ylab="rent")

## ----eval=FALSE----------------------------------------------------------
#  m0<-gamlss(rent~location+bath+kitchen+cheating+area+yearc+pb(area)+pb(yearc),
#             sigma.fo=~area + yearc + pb(area) + pb(yearc),
#             nu.fo=~area + yearc + pb(area) + pb(yearc), family=BCCGo,
#             data=rent99)
#  m1 <- stepGAICAll.A(m0, scope=list(lower=~location+bath+kitchen+cheating+area+
#          yearc+pb(area)+pb(yearc), upper=~(location+bath+kitchen+cheating+
#          area+yearc)^2+pb(area)+pb(yearc)), sigma.scope=list(lower=~area+yearc,
#          upper=~location+bath+kitchen+cheating+area+yearc+pb(area)+pb(yearc)),
#          nu.scope=list(lower=~area+yearc, upper=~location+bath+kitchen+
#          cheating+area+yearc+pb(area)+pb(yearc)),k=4)

## ------------------------------------------------------------------------
fd<-as.factor(rent99$district)
farea<-as.factor(names(rent99.polys))
vizinhos <- polys2nb(rent99.polys)
#creating the precision matrix
precision <- nb2prec(vizinhos,fd,area=farea)
#adding spatial effect for mu

## ----eval=FALSE----------------------------------------------------------
#  m2<- gamlss(formula = rent ~ location + bath + kitchen + cheating +
#              pb(area) + pb(yearc) + cy + ky + ka +
#              gmrf(fd, area = farea, precision = precision, method="A"),
#              sigma.formula = ~area +  pb(yearc) + cheating,
#              nu.formula = ~pb(area) + pb(yearc) +
#              kitchen, family = BCCGo, data = rent99, start.from=m1)

## ----eval=TRUE, echo=FALSE, results='hide'-------------------------------
#load("/Users/MikisStasinopoulos/Dropbox/gamlss/R-code/Spatial/modelm2.Rdata")
load("mrf-rent-data-analysis-final2.Rdata")
#save(m2, file="/Users/MikisStasinopoulos/Dropbox/gamlss/R-code/Spatial/modelm2.Rdata")

## ------------------------------------------------------------------------
AIC(m1,m2, k=2)

## ----eval=FALSE----------------------------------------------------------
#  rent99$nyearc<-rent99$yearc-mean(rent99$yearc)
#  rent99$narea<-rent99$area-mean(rent99$area)
#  m2final<- gamlss(formula = rent ~ location + bath +
#            cheating*nyearc + kitchen*nyearc +
#            kitchen*narea + pb(area) + pb(yearc) +
#            gmrf(fd, area = farea, precision = precision),
#            sigma.formula = ~area  + cheating + pb(yearc),
#            nu.formula = ~  kitchen + pb(area) + pb(yearc),
#            family = BCCGo, data = rent99, start.from=m2)

## ----termplotmu, fig.show='hide', fig.asp=1, warning=FALSE---------------
#to get the termplot for the factors (without interaction)
term.plot(m2final, what="mu", ylim="free")

## ----termplotmu2, fig.show='hide', fig.asp=1, warning=FALSE--------------
#to get the termplot for the continuous variables 
term.plot(m2, what="mu", ylim="free")

## ----termplotmuinteractions, fig.show='hide', fig.asp=1, warning=FALSE----
source("int-term-plot.R")
#to find out the position of the interaction terms
head(lpred(m2final, type="terms")) 
int.term(object=m2final, xvar=rent99$yearc, position=10, 
         fac=rent99$cheating, factor.plots=TRUE, xlabel="yearc", 
         ylabel="cheating", which.lev="0")
int.term(object=m2final, xvar=rent99$yearc, position=11, 
         fac=rent99$kitchen, factor.plots=TRUE, 
         xlabel="yearc", ylabel="kitchen", which.lev="1")
int.term(object=m2final, xvar=rent99$area, position=12, 
         fac=rent99$kitchen, factor.plots=TRUE, 
         xlabel="area", ylabel="kitchen", which.lev="1")

## ----fittedspatial, fig.show='hide', fig.asp=1---------------------------
draw.polys(rent99.polys,getSmo(m2final, what="mu", which=3), 
           scheme="heat")

## ----termplotsigma, fig.show='hide', fig.asp=1---------------------------
term.plot(m2final, what="sigma", ylim="free")

## ----termplotnu, fig.show='hide', fig.asp=1------------------------------
term.plot(m2final, what="nu", ylim="free")

## ------------------------------------------------------------------------
summary(m2final)

## ----plotwp, fig.show='hide', fig.asp=1, warning=FALSE, message=FALSE----
wp(m2final, ylim.all=0.7)

## ----plotvwp, fig.show='hide', fig.asp=1, warning=FALSE, message=FALSE----
wp(m2final, xvar=~yearc+area, n.inter=4, ylim.worm=1)

