
install.packages("gtools")
install.packages("haven")
install.packages("mvtnorm")
install.packages("coda")
#### Cargando librerias
library(gtools)              
library(haven)
library(mvtnorm)
library(coda)



#####################################################################
#
# PREGUNTA 2- Regresión Logistica Simple
#
#####################################################################

credito        <- read.csv(file.choose())
credito        <- credito[,1:10]

names(credito) <- c("id","censor","tiempo","bancarizado","sexo","soltero","edad",     
                    "pcali","tarjeta","linea")  

# Data
model1 <- glm(censor ~ edad,family=binomial(),credito)
model1$coef

y <- credito$censor
n <- dim(credito)[1]
X <- cbind(rep(1,n),credito$edad)  # Matriz de diseño
p <- dim(X)[2]

#### item b)
I     <- 20000
BETA  <- matrix(0,I,p)
PI    <- rep(0,I)
beta  <- rep(0,p) ; var.beta <- mean(y)*(1-mean(y))*solve(t(X)%*%X)

lik   <- function(bet) sum(y*(bet[1] + bet[2]*credito$edad)) - sum(log(1+exp(bet[1] + bet[2]*credito$edad)))

for(i in 1:I)
{
  beta.start <- t(rmvnorm(1,beta,var.beta))
  log.r      <- lik(beta.start) - lik(beta) 
  if( log.r > log(runif(1))) { beta <- beta.start }    
  BETA[i,]  <- beta
  PI[i]     <- exp(beta[1] + beta[2]*30)/(1+exp(beta[1] + beta[2]*30))  
}

#
# Resultados
#

vBETA <- as.mcmc(BETA)
plot(vBETA)
summary(vBETA)

# item d)
quantile(PI,probs=c(0.025,0.975))
PI[i]

# item c)
pred.prob   <- function(x,beta) exp(beta[1] + beta[2]*x)/(1 + exp(beta[1] + beta[2]*x))

credito <- credito[order(credito$edad,decreasing=F),]
plot(credito$edad, pred.prob(credito$edad,beta=colMeans(vBETA)),
     ylab="Probabilidad de caer en morosidad",xlab="Edad",
     ylim=c(0,1),
     type="l")
lines(credito$edad,pred.prob(credito$edad,
                             beta=c(quantile(vBETA[,1],0.025),
                                    quantile(vBETA[,2],0.025))),col=2)
lines(credito$edad,pred.prob(credito$edad,
                             beta=c(quantile(vBETA[,1],0.975),
                                    quantile(vBETA[,2],0.975))),col=2)

