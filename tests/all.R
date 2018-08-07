library(chebpol)

cat("Have FFTW:",havefftw(),'\n')
set.seed(42)

a <- array(rnorm(24),c(x=2,y=3,z=4))
chebcoef(a)

if(havefftw()) {
# A long one-dimensional
f <- function(x) ifelse(x==0,0,sin(1/x))
ch <- ipol(f,dims=50000,method='chebyshev')
cat('Long test:',ch(0.03),f(0.03),ch(0.031),f(0.031),'\n')
}
f <- function(x) exp(-sum(x^2))

dims <- c(8,7,6)
ch <- ipol(f,dims=dims,method='cheb')
s <- runif(3,-1,1)
ch(s)-f(s)
iv <- list(c(1,2),c(1,4),c(1,3))
ch <- ipol(f,dims=dims,intervals=iv)
s <- c(1.4,2.3,1.9)
ch(s) - f(s)
a <- c(1.01,1.01,1.01) ; ch(a)- f(a)

sum(evalongrid(f,dims))
# vector valued function
g <- function(x) c(sum(x),prod(x),exp(-sum(x^2)))
gv <- evalongrid(g,dims)
sum(gv)
dim(gv)


chebknots(17)
chebknots(c(x=3,y=4,z=5))
# test chebappxg
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
ch <- ipol(exp(s),grid=list(s), method='general')
## test it:
r <- runif(1); cat('true:',exp(r),'appx:',ch(r),'\n')

#multivariate chebappxg
su <- seq(0,1,length.out=11)
grid <- list(su,su^2,su^3)
dims <- lapply(grid,length)

fv <- structure(apply(expand.grid(grid),1,f),dim=lapply(grid,length))
ch <- ipol(fv,grid=grid,method='general')
s <- runif(3)
cat('true:',f(s),'appx:',ch(s),'\n')

# multi linear
s <- runif(3)
lip <- ipol(fv,grid=grid,method='multi')
cat('true',f(s), 'appx:', ch(s), 'lip:',lip(s),'\n')

# test dct transform
a <- array(rnorm(24),c(2,3,4))
chebcoef(a,TRUE)

# uniform grid stuff
# Runge function
f <- function(x) 1/(1+25*x^2)
grid <- seq(-1,1,length.out=15)
val <- f(grid)
uc <- ipol(val,method='uniform')
# and the Chebyshev
ch <- ipol(val,dims=15,method='cheb')
# test it at 10 random points
a <- runif(10)
cbind(uc(a),ch(a),sapply(a,f))

uc <- ipol(f,dims=15,intervals=c(-1,1),method='uniform')
a <- runif(1,-1,1)
uc(a); ch(a); f(a)

#polyharmonic splines
f <- function(x) 10/(10+sum(sqrt(x)))
knots <- matrix(runif(6000,0,20), 6)
phs <- ipol(f,knots=knots,k=3,method='poly')
round(fivenum(phs(knots)-apply(knots,2,f)),6)
# test it in a random point
a <- runif(6,0,20)
f(a); phs(a)
phs <- ipol(f,knots=knots,k=5,method='poly',normalize=FALSE)
nphs <- ipol(f,knots=knots,k=5,method='poly',normalize=TRUE)
phs(a); nphs(a)
phs <- ipol(f,knots=knots,k=-1,method='poly',normalize=FALSE)
nphs <- ipol(f,knots=knots,k=-1,method='poly',normalize=TRUE)
phs(a); nphs(a)


#Floater-Hormann

f <- function(x) sum(log(1+x^2))
grid <- list(seq(-1,1,len=50), seq(-1,1,len=30)^3L, chebknots(20)[[1]])
fh <- ipol(f,grid=grid,method='fh',k=c(3,4,4))
x <- matrix(runif(3*1e2),3)
fh.v <- fh(x,2)
f.v <- apply(x,2,f)
round(fivenum(fh.v-f.v),6)
# check some points:
round(f(sapply(grid,function(x) x[4])) - fh(sapply(grid,function(x) x[4]+1e-9)),7)
