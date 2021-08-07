##########################################################
# Run the BART model on the simulation data in Kim (2020+)
##########################################################

# Load reauired libraries
library(MASS)
library(mnormt)
library(MCMCpack)
library(rootSolve)
library(R.matlab)
library(truncnorm)
library(data.table)
library(Rcpp)


#################
# Simulation Data
#################

# Square lattice region
x.easting <- 1:15
x.northing <- 1:15
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid) # num. of locations (20 by 20 = 400)
    
# Distance and neighborhood matrix    
Distance <- distance <- as.matrix(dist(Grid))
diag(Distance) <- 0.1
W <-array(0, c(K,K))
W[distance==1] <-0.1
    
# Define spatial connections (see Kim (2020+) for more details)
wind_mat <- array(0, c(K,K))
wind_seed <- rep(1,15^2)
wind_seed[which(lower.tri(matrix(NA,15,15), diag = FALSE)==TRUE)] <- -1
DIAG <- setdiff(which(lower.tri(matrix(NA,15,15), diag = TRUE)==TRUE),which(lower.tri(matrix(NA,15,15), diag = FALSE)==TRUE))
wind_seed[DIAG] <- 0
wind_seed.alt <- c(1,17,18,32:35, 48:52, 63:69, 79:86, 94:103, 110:120, 125:135, 141:150, 156:165, 172:180, 187:195, 203:210, 218:225)
wind_seed.alt1 <- c(10,11,12,13,14,15,26,27, 28, 29, 30,42,43, 44, 45,58,59, 60,74,75,90)
wind_seed[wind_seed.alt1] <- 2
    
for(i in 1:225){
  value <- wind_seed[i]
  if(value != 0){
    wind_mat[i,which(wind_seed == value)] <- 1
  }
  if(i %in% wind_seed.alt){
    wind_mat[i, wind_seed.alt] <- 1
  }
}
    
    
# Generate the covariates and response data
mean.vector <- rep(-2.5,15^2)
mean.vector[which(lower.tri(matrix(NA,15,15), diag = FALSE)==TRUE)] <- 0.5
DIAG <- setdiff(which(lower.tri(matrix(NA,15,15), diag = TRUE)==TRUE),which(lower.tri(matrix(NA,15,15), diag = FALSE)==TRUE))
mean.vector[DIAG] <- 0
    
x1 <- rnorm(K, 1.5+1*(Grid[,1]/10-1)^2, 1.5)
x2 <- rnorm(K, mean.vector, 1)
x3 <- rnorm(K, 0.5, 1.5)
x.com <- rmnorm(K, c(-1, 2.5), matrix(c(1.5,0.25,0.25,1.5), 2, 2))
x4 <- x.com[,1]
x5 <- x.com[,2]
    
mat <- matrix(-1, 15, 15)
for(i in 1:15){
  for(j in 1:15){
    if(i < j){mat[i,j] <- -2.5}
    if(j < i){mat[i,j] <- 0.5}
  }
}
mat[10:15,1] <- -3.5
mat[11:15,2] <- -3.5
mat[12:15,3] <- -3.5
mat[13:15,4] <- -3.5
mat[14:15,5] <- -3.5
mat[15,6] <- -3.5
    
low <- setdiff(which(lower.tri(matrix(NA,15,15), diag = FALSE)==TRUE), wind_seed.alt1)
up <- which(upper.tri(matrix(NA,15,15), diag = FALSE)==TRUE)
    
# generate spatial random effects
Sigma.alt <- wind_mat
diag(Sigma.alt) <- 2.7
Sigma <- cov2cor(Sigma.alt/Distance)
phi <- identity(wind_seed==2)*rmnorm(1, rep(-3.5,K), 0.5*exp(-0.1 * Distance))
phi[low] <- rmnorm(1, rep(0.5, K), 0.5*exp(-0.1 * Distance))[low]
phi[up] <- rmnorm(1, rep(-2.5, K), 0.5*exp(-0.1 * Distance))[up]
phi[1] <- 0
for(i in wind_seed.alt[-1]){
  mean.temp <- 0 + Sigma[i,-i]%*%solve(Sigma[-i,-i])%*%matrix(c(phi[-i]-c(mat)[-i]), ncol=1)
  sigma.temp <- Sigma[i,i] - Sigma[i,-i]%*%solve(Sigma[-i,-i])%*%t(t(Sigma[i,-i]))
  phi[i] <- rnorm(1, mean.temp, sqrt(sigma.temp))
}
phi <- phi - mean(phi)
    
theta <- rnorm(K, mean=rep(0,K), sd=0.1)
Y.true <- Y <- 2.5+1.5*x1 + 1.2*x2 - 1.25*x3 + 1.5*x4 -1.25*x5 + 0.25*x4*x5 + theta + phi

# randomly selected missing locations
missing <- c(1L, 2L, 3L, 5L, 6L, 7L, 8L, 9L, 11L, 12L, 13L, 14L, 15L, 16L, 
  17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L, 26L, 27L, 28L, 30L, 31L, 
  32L, 33L, 34L, 35L, 36L, 37L, 38L, 40L, 41L, 42L, 43L, 44L, 45L, 
  46L, 47L, 48L, 49L, 50L, 51L, 52L, 53L, 54L, 55L, 56L, 60L, 61L, 
  62L, 63L, 64L, 65L, 67L, 68L, 69L, 70L, 71L, 72L, 73L, 75L, 76L, 
  77L, 78L, 79L, 80L, 81L, 82L, 83L, 84L, 85L, 86L, 87L, 88L, 90L, 
  91L, 92L, 93L, 94L, 95L, 96L, 97L, 99L, 101L, 102L, 103L, 104L, 
  105L, 106L, 107L, 108L, 109L, 111L, 112L, 113L, 114L, 116L, 117L, 
  118L, 120L, 121L, 122L, 123L, 124L, 125L, 126L, 127L, 128L, 129L, 
  131L, 132L, 133L, 134L, 135L, 136L, 137L, 138L, 139L, 140L, 141L, 
  142L, 143L, 144L, 145L, 146L, 147L, 148L, 150L, 151L, 152L, 153L, 
  155L, 156L, 157L, 158L, 159L, 160L, 161L, 163L, 164L, 166L, 167L, 
  168L, 169L, 170L, 172L, 173L, 174L, 175L, 176L, 177L, 178L, 179L, 
  180L, 181L, 183L, 185L, 186L, 187L, 189L, 190L, 192L, 194L, 196L, 
  197L, 198L, 199L, 200L, 201L, 202L, 203L, 204L, 205L, 207L, 208L, 
  209L, 210L, 211L, 212L, 213L, 214L, 215L, 216L, 217L, 218L, 220L, 
  221L, 222L, 223L, 224L, 225L)
Y[missing] <- NA

P <- 50      # Num. of predictors
n <- 225

cov <- list()
for(i in 1:45){
  cov[[i]] <- rnorm(n,0,1) # draw x independently from a normal distribution
}
Xpred1 <- do.call(cbind, cov)

# Vector of X predictors
Xpred <- cbind(x1,x2,x3,x4,x5, Xpred1)


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

P <- dim(Xpred)[2]  # Num. of predictors
Xcut <- lapply(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t]))) # unique values of the predictors
shift <- mean(Y, na.rm=TRUE)
Y <- Y - shift      # Shifting the mean of Y

A.WEIGHT <- NULL    # uncertainty parameter 
n.full <- length(Y)  # Num. of the locations
n.complete <- length(which(!is.na(Y)))  # Num. of the locations with observation
m <- 100 # Num. of the trees in the BART model (default:100)
n.iter <- 17000 # Num. of MCMC iterations

# Create empty objects
Tree1 <- matrix(0,nrow=n.complete, ncol=m)  # tree for complete model
Effect <- matrix(nrow=n.full, ncol=n.iter)
TREE1 <- list()            # Predicted mean Values
TREE1[[1]] <- list(T=matrix(nrow=n.complete,ncol=m))
MU1 <- list()              # Mean parameters associated with terminal nodes
MU1[[1]] <- list()
I.nodes1 <- list()         # Place holder of internal nodes (predcitors) + terminal nodes (0's)
I.nodes1[[1]] <- list()
Rules1 <- list()           # Place holder of internal nodes (rules)
Rules1[[1]] <- list()
Prob <- list()
pr <- matrix(nrow=P, ncol=1)
Prob[[1]] <- pr 
Sigma2 <- NULL            # Variance parameter
ind <- matrix(nrow=1, ncol=P)

### New place holder for the next iteration
TREE1[[2]] <- TREE1[[1]]
TREE1[[2]] <- TREE1[[1]]
MU1[[2]] <- MU1[[1]]
I.nodes1[[2]] <- I.nodes1[[1]]
Rules1[[2]] <- Rules1[[1]]
Prob[[2]] <-list(pr=Prob[[1]])


# Set (hyper-)parameters and initial values
prior.tau2 <- c(1, 0.01)
tau2.posterior.shape <- prior.tau2[1] + 0.5*n.full
proposal.sd.rho <- 0.02
rho <- 0.5
spatial <- rnorm(length(Y), 0, 1)
tau2 <- 0.1
R <- Y  # Initial values of R


# Missing index & Obj for imputed outcome
mis.ind <- which(is.na(Y))
Y.da <- matrix(nrow=length(mis.ind), ncol=n.iter)
Y.da[,1] <- rnorm(length(mis.ind), mean(Y, na.rm=TRUE), 1)
Y[mis.ind] <- Y.da[,1]
count <- 0
sigma2 <- 1
Sigma2[1] <- sigma2
sigma_mu <- max((min(Y)/(-2*sqrt(m)))^2, (max(Y)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y

# Set 'nu' and 'q' parameters (Chipman et al. 2010)
nu <- 3                   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
lambda <- uniroot.all(f, c(0.1^5,10))

dir.alpha <- 1 # alpha parameter in the prior for selection probabilities
post.dir.alpha <- rep(dir.alpha, P)
alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
beta <- 2                 # default setting (alpha, beta) = (0.95, 2)
prop.prob <- rdirichlet(1, post.dir.alpha)
Prob[[2]]$pr[,1] <- prop.prob


# Prob. of selecting an alteration (GROW, PRUNE and CHANGE)
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE



### Initial Setup for the "m" trees (depth = 0 and 1)
for(t in 1:m){
  i.nodes1 <- list()
  for(i in 1:2){
    i.nodes1[[i]] <- rep(NA, 2^(i-1))
  }
  i.nodes1[[1]] <- 999
  I.nodes1[[1]][[t]] <- i.nodes1
  
  rules1 <- list()
  for(i in 1:2){
    rules1[[i]] <- rep(NA, 2^(i-1))
  }
  Rules1[[1]][[t]] <- rules1
}

# Weight matrix
a.weight <- 1 # initial value for a_w parameter in the weight
diag(distance) <- NA 
W2 <- 1/distance
diag(W2) <- 0
W2 <- ifelse(W2 == Inf, 0, W2)

W2.wind <- wind_mat*a.weight*W2+(1-wind_mat)*(1-a.weight)*W2
W2.wind <- W2.wind[-mis.ind, -mis.ind]
rownames(W2.wind) <- 1:(n.full-length(mis.ind))
colnames(W2.wind) <- 1:(n.full-length(mis.ind))
W2.post <- common.Wcheckformat(W2.wind)
W2.wind.full <- wind_mat*a.weight*W2+(1-wind_mat)*(1-a.weight)*W2
W2.post.full <- common.Wcheckformat(W2.wind.full)

det.Q <- 0
Wstar <- 0
Wstar.eigen <- 0
Wstar.val <- 0

#####################################################################
############# Run main MCMC #########################################
#####################################################################

for(j in 2:n.iter){
  count <- count + 1
  
  for(t in 1:m){
    
    R <- Y[-mis.ind] - rowSums(Tree1[,-t]) - spatial[-mis.ind]
    i.nodes <- I.nodes1[[1]][[t]]
    rules <- Rules1[[1]][[t]]
    
    ### Find the depth of the tree (0 or 1 or 2)
    tree.length <- length(which(!is.na(unlist(rules))))
    if(tree.length == 0){  # tree has no node yet
      grow.step <- GROW.first(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
      I.nodes1[[2]][[t]] <- grow.step$i.nodes
      Rules1[[2]][[t]] <- grow.step$rules
    } else {
      step <- sample(1:3, 1, prob=c(p.grow, p.prune, p.change))  # Pick a step
      if(step==1){  # GROW step
        grow.step <- GROW(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
        I.nodes1[[2]][[t]] <- grow.step$i.nodes
        Rules1[[2]][[t]] <- grow.step$rules
      }
      if(step==2){   # PRUNE step
        prune.step <- PRUNE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
        I.nodes1[[2]][[t]] <- prune.step$i.nodes
        Rules1[[2]][[t]] <- prune.step$rules
      }
      if(step==3){   # CHANGE step
        change.step <- CHANGE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
        I.nodes1[[2]][[t]] <- change.step$i.nodes
        Rules1[[2]][[t]] <- change.step$rules
      }
    }
    
    Mean <- Mean.Parameter(Sigma2[j-1], sigma_mu, i.nodes=I.nodes1[[2]][[t]], rules=Rules1[[2]][[t]],  R, ind=1)
    T <- Mean$T
    mu <- Mean$mu
    
    eval(parse(text=(paste("TREE1[[",2,"]]$T[,",t,"]=T",sep=""))))
    eval(parse(text=(paste("Tree1[,",t,"]=T",sep=""))))
    eval(parse(text=(paste("MU1[[",2,"]]$mu",t,"=mu",sep=""))))
  }
  
  # Sample variance parameter
  Sigma2[j] <- 1 / rgamma(1, 1+0.5*n.complete, scale=1/(0.01+0.5*sum((Y[-mis.ind]-rowSums(Tree1)-spatial[-mis.ind])^2)))

  #######################################
  # Start update the spatial random effect
  #######################################  
  offset <- (Y[-mis.ind]-rowSums(Tree1))
  spatial[-mis.ind] <- gaussiancarupdate(Wtriplet=W2.post$W.triplet, Wbegfin=W2.post$W.begfin, W2.post$W.triplet.sum, nsites=length(Y[-mis.ind]), phi=spatial[-mis.ind], tau2=tau2, rho=rho, nu2=Sigma2[j], offset=offset)

  temp2 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], rho)
  tau2.posterior.scale <- temp2 + prior.tau2[2] 
  tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
  
  Wstar <- diag(apply(W2.wind,1,sum)) - W2.wind
  Wstar.eigen <- eigen(Wstar)
  Wstar.val <- Wstar.eigen$values
  det.Q <- 0.5 * sum(log((rho * Wstar.val + (1-rho))))    
  
  # update rho parameter
  proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)  
  temp3 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], proposal.rho)
  tau3.posterior.scale <- temp3 + prior.tau2[2]
  tau3 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau3.posterior.scale))
  det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
  logprob.current <- det.Q - temp2 / tau2
  logprob.proposal <- det.Q.proposal - temp3 / tau3
  hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
  prob <- exp(logprob.proposal - logprob.current + hastings)
  #### Accept or reject the proposal
  if(prob > runif(1))
  {
    rho <- proposal.rho
    det.Q <- det.Q.proposal
  }              
  
  for(i in seq_along(mis.ind)){
    spatial[mis.ind[i]] <- rnorm(1, rho*sum(W2.wind.full[mis.ind[i], mis.ind[i]]*spatial[-mis.ind[i]])/(rho*W2.post.full$W.triplet.sum[mis.ind[i]]+1-rho),  sqrt(tau2 /(rho*W2.post.full$W.triplet.sum[mis.ind[i]]+1-rho) ) )
  }

  # update a_w parameter in the weights  
  proposal.a.weight <- rbeta(1, (a.weight-0.1^5)/(1-a.weight+0.1^5)*10, 10)
  
  W2.wind.proposal <- wind_mat*proposal.a.weight*W2+(1-wind_mat)*(1-proposal.a.weight)*W2
  W2.wind.proposal <- W2.wind.proposal[-mis.ind, -mis.ind]
  rownames(W2.wind.proposal) <- 1:(n.full-length(mis.ind))
  colnames(W2.wind.proposal) <- 1:(n.full-length(mis.ind))
  W2.post.proposal <- common.Wcheckformat(W2.wind.proposal)
  
  temp2 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], rho)
  temp3 <- quadform(as.matrix(W2.post.proposal$W.triplet),W2.post.proposal$W.triplet.sum , W2.post.proposal$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], rho)
  tau3.posterior.scale <- temp3 + prior.tau2[2]
  tau3 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau3.posterior.scale))
  
  Wstar.proposal <- diag(apply(W2.wind.proposal,1,sum)) - W2.wind.proposal
  Wstar.eigen.proposal <- eigen(Wstar.proposal)
  Wstar.val.proposal <- Wstar.eigen.proposal$values
  
  det.Q.proposal <- 0.5 * sum(log((rho * Wstar.val.proposal + (1-rho))))              
  logprob.current <- det.Q - temp2 / tau2
  logprob.proposal <- det.Q.proposal - temp3 / tau3
  hastings <- log(dbeta(proposal.a.weight, 2,1))-log(dbeta(a.weight, 2, 1))+log(dbeta(a.weight, (proposal.a.weight-0.1^5)/(1-proposal.a.weight+0.1^5)*10, 10)) - log(dbeta(x=proposal.a.weight, (a.weight-0.1^5)/(1-a.weight+0.1^5)*10, 10)) 
  prob <- exp(logprob.proposal - logprob.current + hastings)
  
  #### Accept or reject the proposal
  if(prob > runif(1))
  {
    a.weight <- proposal.a.weight
    det.Q <- det.Q.proposal
    W2.wind <- W2.wind.proposal
    W2.post <- W2.post.proposal
    Wstar <- Wstar.proposal
    Wstar.eigen <- Wstar.eigen.proposal
    Wstar.val <- Wstar.val.proposal
  }              
  A.WEIGHT[j] <- a.weight
  
  # Update alpha parameter in the Dirichlet prior for the selection probabilities
  add <- as.numeric(table(factor(unlist(I.nodes1[[2]]), levels = 1:P)))
  if(j<n.iter/10){
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
  TREE1[[1]] <- TREE1[[2]]
  MU1[[1]] <- MU1[[2]]
  I.nodes1[[1]] <- I.nodes1[[2]]
  Rules1[[1]] <- Rules1[[2]]
  Prob[[1]] <-Prob[[2]]
  
  # Predict the outcomes based on the current parameters
  Tree11 <- matrix(0,nrow=n.full, ncol=m)
  for(t in 1:m){
    i.nodes <- I.nodes1[[2]][[t]]
    rules <- Rules1[[2]][[t]]
    eval(parse(text=(paste("mu <- MU1[[",2,"]]$mu",t,sep=""))))
    T <- Mean.Parameter_pred(sigma2=Sigma2[j], sigma_mu=sigma_mu, i.nodes=i.nodes, rules=rules,  mu=mu, ind=3, xpred=Xpred)$T
    eval(parse(text=(paste("Tree11[,",t,"]=T",sep=""))))
  }
  
  # Data augmentation
  Y.da[,count] <- rnorm(length(mis.ind), c(rowSums(Tree11)+spatial)[mis.ind], sqrt(Sigma2[j]))
  Y[mis.ind] <- Y.da[,count]
  ind.temp<-ifelse(add > 0, 1, 0)
  ind <- rbind(ind, ind.temp)
}

# Post-processing; burn-in the first 7000 iterations and shift the estimates by 'shift'.
mat <- matrix(nrow=n.full, ncol=10000)
for(i in 1:10000){
  Y[mis.ind] <- Y.da[,i+6999]
  mat[,i] <- Y+shift
}
ind <- ind[7001:17000,]
mat <- rowMeans(mat)












