##########################################################
# Run the BART model on the simulation data in Kim (2020+)
##########################################################

# Load required libraries
library(MCMCpack)
library(rootSolve)
library(R.matlab)
library(truncnorm)
library(data.table)
library(Rcpp)
library(MASS)
library(plot.matrix)
library(mnormt)
library(reshape)

#################
# Simulation Data
#################

# Square lattice region
x.easting <- 1:20
x.northing <- 1:20
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid) # num. of locations (20 by 20 = 400)

# Distance and neighborhood matrix
Distance <- as.matrix(dist(Grid))
diag(Distance) <- 1
W <-array(0, c(K,K))
W[Distance==1] <- 0.1

# Define spatial connections (see Kim (2020+) for more details)
W.hole <- array(0.1, c(K,K))
for(i in 1:400){
	temp <- Grid[i,]
	W.hole[i,which((Grid[,1] > c(temp[1])) & (Grid[,2] > c(temp[2])))] <- 100
}
W.hole <- W.hole + t(W.hole)
W.hole <- ifelse(W.hole==0.2, 1000, 0.8)
diag(W.hole) <- 0


# Generate covariates and response data
x1 <- rnorm(K, -1+1*(Grid[,1]/10-1)^2, 2)
x2 <- rnorm(K, 1-1*(Grid[,2]/10-1)^2, 1)
x.multivariate <- rmnorm(K, c(-1, 1), 0.1*matrix(c(1,0.5,0.5,1), 2, 2))
x3 <- x.multivariate[,1]
x4 <- x.multivariate[,2]
x5 <- rpois(K, 2)
theta <- rnorm(K, mean=rep(0), sd=0.5)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=exp(-1 * Distance*W.hole)) # spatial effect
Y <- 5+1*x1 + 1*x2 + 1*x3 + 1.5*abs(x4) + 0.2*x4*x5 + theta + phi

# Randomly set missing responses (340 locations)
missing <- sample(1:400, 340)
Y[missing] <- NA

P <- 50      # Num. of predictors (including x1-x5)
n <- 400

cov <- list()
for(i in 1:20){
  cov[[i]] <- rnorm(n,0,1)  # draw x independently from a normal distribution
}
Xpred1 <- do.call(cbind, cov)

cov <- list()
for(i in 1:20){
  cov[[i]] <- rpois(n,i/10) # draw x independently from a poisson distribution
}
Xpred2 <- do.call(cbind, cov)

# Vector of X predictors
Xpred <- cbind(x1,x2,x3,x4,x5, 
               x1[sample(1:400)], x2[sample(1:400)], x3[sample(1:400)], x4[sample(1:400)], x5[sample(1:400)], 
               Xpred1, Xpred2)



####################
# Set the BART model
####################

# Import files for the spatial random effect
sourceCpp("CARBayes.cpp") # this C++ code from CARBayes
source("commonfunctions.R")

# Import files for the BART model
source("GROW.R")
source("PRUNE.R")
source("CHANGE.R")
source("Common.R")
source("Prediction.R")

P <- dim(Xpred)[2] # Num. of predictors
Xcut <- lapply(1:P, function(t) sort(unique(Xpred[,t]))) # unique values of the predictors
shift <- mean(Y, na.rm=TRUE)
Y <- Y - shift  # Shifting Y by the mean
n.full <- length(Y) # Num. of the locations
m <- 100 # Num. of the trees in the BART model (default:100)
n.iter <- 15000 # Num. of MCMC iterations

# Create empty objects
Tree <- matrix(0,nrow=n.full, ncol=m)
Effect <- matrix(nrow=n.full, ncol=n.iter)
TREE <- list()            # Predicted mean Values
TREE[[1]] <- list(T=matrix(nrow=n.full, ncol=m))
MU <- list()              # Mean parameters associated with terminal nodes
MU[[1]] <- list()
I.nodes <- list()         # Place holder for internal nodes (predcitors) + terminal nodes (0's)
I.nodes[[1]] <- list()
Rules <- list()           # Place holder for internal nodes (rules)
Rules[[1]] <- list()
Prob <- list()
pr <- matrix(nrow=P, ncol=1)
Prob[[1]] <- pr   
Sigma2 <- NULL            # Variance parameter
ind <- matrix(nrow=1, ncol=P)  # place holder for inclusion prob.

### New place holder for the next iteration
TREE[[2]] <- TREE[[1]]
MU[[2]] <- MU[[1]]
I.nodes[[2]] <- I.nodes[[1]]
Rules[[2]] <- Rules[[1]]
Prob[[2]] <-list(pr=Prob[[1]])


# Set (hyper-)parameters and initial values
prior.tau2 <- c(1, 0.01)
tau2.posterior.shape <- prior.tau2[1] + 0.5*n.full
proposal.sd.rho <- 0.02
spatial <- rnorm(n.full, 0, 0.01) 
rho <- 0.5
tau2 <- 0.1
R <- Y  # Initial values of R

# Missing index & Obj for imputed outcome
mis.ind <- which(is.na(Y))
Y.da <- matrix(nrow=length(mis.ind), ncol=15000)
Y.da[,1] <- rnorm(length(mis.ind), mean(Y, na.rm=TRUE), 0.01) # initial values
Y[mis.ind] <- Y.da[,1]

residuals <- lm(Y~Xpred)$residuals
sigma2 <- var(residuals)  # Initial value of SD^2
Sigma2[1] <- sigma2
sigma_mu <- max((min(Y)/(-2*sqrt(m)))^2, (max(Y)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y

# Set 'nu' and 'q' parameters (Chipman et al. 2010)
nu <- 3   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
lambda <- uniroot.all(f, c(0.1^5,10))

dir.alpha <- 5 # alpha parameter in the prior for selection probabilities
post.dir.alpha <- rep(dir.alpha, P)
alpha <- 0.95     # alpha (1+depth)^{-beta} where depth=0,1,2,...
beta <- 2         # default setting (alpha, beta) = (0.95, 2)
prop.prob <- rdirichlet(1, post.dir.alpha) # selection probabilities
Prob[[2]]$pr[,1] <- prop.prob

# Prob. of selecting an alteration (GROW, PRUNE and CHANGE)
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE
 

### Initial Setup for the "m" trees (depth = 0 and 1)
for(t in 1:m){
  i.nodes <- list()
  for(i in 1:2){
    i.nodes[[i]] <- rep(NA, 2^(i-1))
  }
  i.nodes[[1]] <- 999
  I.nodes[[1]][[t]] <- i.nodes
  
  rules <- list()
  for(i in 1:2){
    rules[[i]] <- rep(NA, 2^(i-1))
  }
  Rules[[1]][[t]] <- rules
}


# Weight matrix
a.weight <- 0.6 # initial value for a_w parameter in the weight
diag(Distance) <- NA 
W <- 1/Distance
W2 <- 1/Distance^2
diag(W) <- 0
diag(W2) <- 0
W <- ifelse(W == Inf, 0, W)
W2 <- ifelse(W2 == Inf, 0, W2)

wind_mat <- array(1, c(K,K))
for(i in 1:400){
  temp <- Grid[i,]
  wind_mat[i,which((Grid[,1] < c(temp[1])) | (Grid[,2] < c(temp[2])))] <- 0
}
wind_mat <- wind_mat + t(wind_mat)
wind_mat <- ifelse(wind_mat==0, 0, 1)

W2.wind <- wind_mat*a.weight*W2+(1-wind_mat)*(1-a.weight)*W2
W2.post <- common.Wcheckformat(W2.wind)
det.Q <- 0
Wstar <- 0
Wstar.eigen <- 0
Wstar.val <- 0


#####################################################################
############# Run main MCMC #########################################
#####################################################################

for(j in 2:n.iter){
  # update tree 't'
  for(t in 1:m){
    R <- Y - rowSums(Tree[,-t]) - spatial
    i.nodes <- I.nodes[[1]][[t]]
    rules <- Rules[[1]][[t]]
    
    ### Find the depth of the tree (0 or 1 or 2)
    tree.length <- length(which(!is.na(unlist(rules))))
    if(tree.length == 0){  # tree has no node yet. need to grow
        grow.step <- GROW.first(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes[[1]][[t]],rules=Rules[[1]][[t]], R=R, prop.prob=prop.prob)
        I.nodes[[2]][[t]] <- grow.step$i.nodes
        Rules[[2]][[t]] <- grow.step$rules
    } else {
      step <- sample(1:3, 1, prob=c(p.grow, p.prune, p.change))  # Pick a step
      if(step==1){  # GROW step
        grow.step <- GROW(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes[[1]][[t]],rules=Rules[[1]][[t]], R=R, prop.prob=prop.prob)
        I.nodes[[2]][[t]] <- grow.step$i.nodes
        Rules[[2]][[t]] <- grow.step$rules
      }
      if(step==2){   # PRUNE step
        prune.step <- PRUNE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes[[1]][[t]],rules=Rules[[1]][[t]], R=R, prop.prob=prop.prob)
        I.nodes[[2]][[t]] <- prune.step$i.nodes
        Rules[[2]][[t]] <- prune.step$rules
      }
      if(step==3){   # CHANGE step
        change.step <- CHANGE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes[[1]][[t]],rules=Rules[[1]][[t]], R=R, prop.prob=prop.prob)
        I.nodes[[2]][[t]] <- change.step$i.nodes
        Rules[[2]][[t]] <- change.step$rules
      }
    }
    
    Mean <- Mean.Parameter(Sigma2[j-1], sigma_mu, i.nodes=I.nodes[[2]][[t]], rules=Rules[[2]][[t]],  R)
    T <- Mean$T
    mu <- Mean$mu
    
    eval(parse(text=(paste("TREE[[",2,"]]$T[,",t,"]=T",sep=""))))
    eval(parse(text=(paste("Tree[,",t,"]=T",sep=""))))
    eval(parse(text=(paste("MU[[",2,"]]$mu",t,"=mu",sep=""))))
  }
  
  # Sample variance parameter
  Sigma2[j] <- rinvgamma(1, nu/2+n.full/2, scale = nu*lambda/2 + sum((Y-rowSums(Tree)-spatial)^2)/2)
  
  #######################################
  # Star update the spatial random effect
  #######################################
  offset <- Y-rowSums(Tree)
  spatial <- gaussiancarupdate(Wtriplet=W2.post$W.triplet, Wbegfin=W2.post$W.begfin, W2.post$W.triplet.sum, nsites=length(Y), phi=spatial, tau2=tau2, rho=rho, nu2=Sigma2[j], offset=offset)
  spatial <- spatial - mean(spatial)

  temp2 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y), spatial, spatial, rho)
  tau2.posterior.scale <- temp2 + prior.tau2[2] 
  tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale)) # update tau^2 parameter
  
  Wstar <- diag(apply(W2.wind,1,sum)) - W2.wind
  Wstar.eigen <- eigen(Wstar)
  Wstar.val <- Wstar.eigen$values
  det.Q <- 0.5 * sum(log((rho * Wstar.val + (1-rho))))    
  
  # update rho parameter
  proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)  
  temp.proposal <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y), spatial, spatial, proposal.rho)
  det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
  logprob.current <- det.Q - temp2 / tau2
  logprob.proposal <- det.Q.proposal - temp.proposal / tau2
  hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
  prob <- exp(logprob.proposal - logprob.current + hastings)
  # Accept or reject the proposal
  if(prob > runif(1)){
    rho <- proposal.rho
    det.Q <- det.Q.proposal
  }              

  
  # update a_w parameter in the weights
  proposal.a.weight <- rbeta(1, a.weight/(1-a.weight)*2, 2) 
  
  W2.wind.proposal <- wind_mat*proposal.a.weight*W2+(1-wind_mat)*(1-proposal.a.weight)*W2 # new weights 
  W2.post.proposal <- common.Wcheckformat(W2.wind.proposal)
  
  temp2 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y), spatial, spatial, rho)
  temp3 <- quadform(as.matrix(W2.post.proposal$W.triplet), W2.post.proposal$W.triplet.sum, W2.post.proposal$n.triplet, length(Y), spatial, spatial, rho)
  
  Wstar.proposal <- diag(apply(W2.wind.proposal,1,sum)) - W2.wind.proposal
  Wstar.eigen.proposal <- eigen(Wstar.proposal)
  Wstar.val.proposal <- Wstar.eigen.proposal$values
  
  det.Q.proposal <- 0.5 * sum(log((rho * Wstar.val.proposal + (1-rho))))              
  logprob.current <- det.Q - temp2 / tau2
  logprob.proposal <- det.Q.proposal - temp3 / tau2
  hastings <- log(dbeta(proposal.a.weight, 2,1))-log(dbeta(a.weight, 2, 1))+log(dbeta(a.weight, proposal.a.weight/(1-proposal.a.weight)*2, 2)) - log(dbeta(x=proposal.a.weight, a.weight/(1-a.weight)*2, 2)) 
  prob <- exp(logprob.proposal - logprob.current + hastings)
  
  #### Accept or reject the proposal
  if(prob > runif(1)){
    a.weight <- proposal.a.weight
    det.Q <- det.Q.proposal
    W2.wind <- W2.wind.proposal
    W2.post <- W2.post.proposal
    Wstar <- Wstar.proposal
    Wstar.eigen <- Wstar.eigen.proposal
    Wstar.val <- Wstar.val.proposal
  }              
  
  # Update alpha parameter in the Dirichlet prior for the selection probabilities
  add <- as.numeric(table(factor(unlist(I.nodes[[2]]), levels = 1:P)))
  if(j<n.iter/10){ # warm-up
    post.dir.alpha <- rep(1, P) + add 
  }else{
    p.dir.alpha <- max(rnorm(1, dir.alpha, 0.1), 0.1^10)
    SumS <-  log(ifelse(Prob[[2]]$pr<0.1^300, 0.1^300, Prob[[2]]$pr))
    dir_lik.p <- sum(SumS*(rep(p.dir.alpha/P, P)-1)) + lgamma(sum(rep(p.dir.alpha/P, P)))-sum(lgamma(rep(p.dir.alpha/P, P))) 
    dir_lik <- sum(SumS*(rep(dir.alpha/P, P)-1)) + lgamma(sum(rep(dir.alpha/P, P)))-sum(lgamma(rep(dir.alpha/P, P))) 
    ratio <- dir_lik.p + log((p.dir.alpha/(p.dir.alpha+P))^(0.5-1)*(P/(p.dir.alpha+P))^(1-1)*abs(1/(p.dir.alpha+P)-p.dir.alpha/(p.dir.alpha+P)^2 ) ) + dnorm(dir.alpha, p.dir.alpha,0.1, log=TRUE) - dir_lik - log((dir.alpha/(dir.alpha+P))^(0.5-1)*(P/(dir.alpha+P))^(1-1)*abs(1/(dir.alpha+P)-dir.alpha/(dir.alpha+P)^2 ) ) - dnorm(p.dir.alpha, dir.alpha,0.1, log=TRUE)
    if (log(runif(1))<ratio) {
      dir.alpha <- p.dir.alpha
    }
    post.dir.alpha <- rep(dir.alpha/P, P) + add 
  }
  
  # update selection prob.
  prop.prob <- rdirichlet(1, post.dir.alpha)
  Prob[[2]]$pr[,1] <- prop.prob
  
  # Store the updates
  TREE[[1]] <- TREE[[2]]
  MU[[1]] <- MU[[2]]
  I.nodes[[1]] <- I.nodes[[2]]
  Rules[[1]] <- Rules[[2]]
  Prob[[1]] <-Prob[[2]]
  
  # Predict the outcomes based on the current parameters
  Tree11 <- matrix(0,nrow=n.full, ncol=m)
  for(t in 1:m){
    i.nodes <- I.nodes[[2]][[t]]
    rules <- Rules[[2]][[t]]
    eval(parse(text=(paste("mu <- MU[[",2,"]]$mu",t,sep=""))))
    T <- Mean.Parameter_pred(sigma2=Simga2[j], sigma_mu=sigma_mu, i.nodes=i.nodes, rules=rules,  mu=mu, ind=3, xpred=Xpred)$T
    eval(parse(text=(paste("Tree11[,",t,"]=T",sep=""))))
  }
  
  # Data augmentation
  Y.da[,(j-1)] <- rnorm(length(mis.ind), c(rowSums(Tree11)+spatial)[mis.ind], sqrt(Sigma2[j]))
  Y[mis.ind] <- Y.da[,(j-1)]
  ind.temp<-ifelse(add > 0, 1, 0)
  ind <- rbind(ind, ind.temp) # to calculate the inclusion prob.
}

# Post-processing; burn-in the first 10000 iterations and shift the estimates by 'shift'.
mat <- matrix(nrow=n.full, ncol=5000)
for(i in 1:5000){
    Y[mis.ind] <- Y.da[,i+9999]
    mat[,i] <- Y+shift
}
mat <- rowMeans(mat)
ind <- ind[10001:15000,] 

