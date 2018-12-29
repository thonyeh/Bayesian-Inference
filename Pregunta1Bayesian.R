
#### Pregunta 1######
### Distribucion Rayleigh
  #-----------------------
### Cargando librerias
library(VGAM)

# distribución Rayleigh
f <- function(x, sigma) {
  #if (any(x < 0)) return (0) 
  #stopifnot(sigma > 0)      
  (x/sigma^2)*exp(-x^2/(2*sigma^2))
}

#set.seed(1989) 

N <-10000       
x <- numeric(N)
sigma <- 2

## Valor inicial
x[1] <- 0.5

## 1.a

for(i in 2:N){
  #Valor anterior
  x.start <- x[i-1]
  
  # Proponer
  x.prop <-rchisq(1,df=x.start)
  
  # Proceso de aceptacion
  px.start <- dchisq(x.start,df=x.prop)
  r <- min(1,f(x.prop,sigma)*dchisq(x.prop,df=x.start)/(f(x.start,sigma)*dchisq(x.prop,df=x.start)))
  u=runif(1)
  if(r > u){x[i]<-x.prop}
  if(r <= u){x[i]<-x[i-1]}
}

acf(x)
x <- x[-c(1,1000)]


par(mfrow=c(1,2))
hist(x,prob=T,main=paste("Histograma of x"))

x_r = seq(0,6,0.1)
curve(drayleigh(x_r, scale=2), add=TRUE, col="blue")
ts.plot(z,col=2,main=paste("Series de tiempo de z"))

#### Quemando Cadenas(1000==10lag) o burning 
burn <- 1001      
z <- z[burn: N] 
par(mfrow=c(1,2))
hist(z,prob=T,main=paste("Histograma of z"))
curve(drayleigh(x, scale=2), add=TRUE, col="blue")
ts.plot(z,col=2,main=paste("Series de tiempo de z"))


### item b)
### Usando Normal truncada en [a,b]
#-----------------------
### Cargando paquete
library(msm)

set.seed(1990)  
N <-10000       
sigma <- 2
z<- numeric(N)
z[1] <- rtnorm(1,mean=1,sd=1,lower=0,upper=Inf) #Valor inicial

#### Corriendo algoritmo Metropolis-Hastings
for(h in 2:N){
  xt<-z[h-1]      # valor anterior
  y<-rtnorm(1,mean=xt, sd=1,lower=0, upper=Inf) # valor candidato
  u<-runif(1)
  alpha<-min(1,(f(y,sigma))/(f(xt,sigma)))
  # prob. de aceptar al candidato
  if(u<=alpha){z[h]<-y} #acepta al candidato
  if(u>alpha){z[h]<-xt} #rechaza al candidato
  
}

acf(z)

par(mfrow=c(1,2))
hist(z,prob=T,main=paste("Histograma de z"))
curve(drayleigh(x, scale=2), add=TRUE, col="blue")
ts.plot(z,col=2,main=paste("Series de tiempo de z"))

#------------------------------------------------------------------------------#
#c) Distribucion Gamma como función generadora 

set.seed(19910) 
N <-10000       
### Valor inicial 
z<- numeric(N)
z[1]  <-1       
sigma <-2
eta1  <-1
eta2  <-2

#### Corriendo Algoritmo Metropolis-Hastings
for(h in 2:N){
  xt<-z[h-1]      # valor anterior
  y<-rgamma(1,shape = eta1,
            rate = eta1/xt) # valor candidato
  u<-runif(1)
  r<-min(1,(f(y,sigma)*dgamma(1,shape = eta1,
                              rate = eta1/y))/
           (f(xt,sigma)*dgamma(1,shape = eta1,
                               rate = eta1/xt)))
  # prob. de aceptar al candidato
  if(u<=r){z[h]<-y} #acepta al candidato
  if(u>r){z[h]<-xt} #rechaza al candidato
  
}

acf(z)

par(mfrow=c(1,2))
hist(z,prob=T,main=paste("Histograma of z"))
curve(drayleigh(x, scale=2), add=TRUE, col="blue")
ts.plot(z,col=2,main=paste("Series de tiempo de z"))
mean(z)

#------------------------------------------------------------------------------

#### Corriendo Algoritmo Metropolis-Hastings
for(h in 2:N){
  xt<-z[h-1]      # valor anterior
  y<-rgamma(1,shape = eta2,
            scale = eta2/xt) # valor candidato
  u<-runif(1)
  r<-min(1,(f(y,sigma)*dgamma(1,shape = eta2,
                              scale = eta2/y))/
           (f(xt,sigma)*dgamma(1,shape = eta2,
                               scale = eta2/xt)))
  # prob. de aceptar al candidato
  if(u<=r){z[h]<-y} #acepta al candidato
  if(u>r){z[h]<-xt} #rechaza al candidato
  
}

acf(z)

par(mfrow=c(1,2))
hist(z,prob=T,main=paste("Histograma of z"))
curve(drayleigh(x, scale=2), add=TRUE, col="blue")
ts.plot(z,col=2,main=paste("Series de tiempo de z"))
mean(z)











