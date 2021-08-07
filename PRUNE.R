# Fun. of PRUNE alteration
PRUNE <-  function(sigma2, sigma_mu, i.nodes, rules,  R, prop.prob, ind=NULL){
    
    if(ind==1){
        xpred <- Xpred[-mis.ind,]
        xcut <- lapply(1:dim(xpred)[2], function(t) sort(unique(xpred[,t]))) # e.g. unique values of predictors
  }else{
    if(ind==0){
      xpred <- Xpred0; xcut <- Xcut0
    }else{
      xpred <- Xpred; xcut <- Xcut
    }
  }
  
  
    t.node <- c(0, 0)
    for(i in 2:length(i.nodes)){
        temp <- grep(999, i.nodes[[i]])
        if(length(temp) > 0){
            t.node <- rbind(t.node, cbind(i, temp))
        }
    }
    t.node <- t.node[-1,]
    
    # find nodes with two terminal child nodes (singly internal parent nodes)
    count <- 0
    count.ind <- NULL
    for(i in 2:dim(t.node)[1]){
        if(t.node[i,2] %% 2 == 0){if(t.node[i-1,2] == (t.node[i,2]-1) & t.node[i, 1] == t.node[(i-1),1]){count <- count + 1; count.ind[count] <- i}}
    }
  
    singly.inode <- count.ind[sample(1:count, 1, replace=FALSE, prob=rep(1/count, count))] # pick a singly internal parent node
    prop.i.nodes <- i.nodes # a new tree
    prop.i.nodes[[t.node[singly.inode,1]]][t.node[singly.inode,2]-1] <-NA
    prop.i.nodes[[t.node[singly.inode,1]]][t.node[singly.inode,2]] <-NA
    prop.i.nodes[[t.node[singly.inode,1]-1]][t.node[singly.inode,2]/2] <- 999
  
    prop.rules <- rules
    prop.rules[[t.node[singly.inode,1]-1]][t.node[singly.inode,2]/2] <- NA
  
  
    prop.t.node <- c(0, 0)
    for(i in 2:length(prop.i.nodes)){
        temp <- grep(999, prop.i.nodes[[i]])
        if(length(temp) > 0){
            prop.t.node <- rbind(prop.t.node, cbind(i, temp))
        }
    }
  
    prop.rule <- NULL
    if(length(prop.t.node)==2){
        prop.pred <- i.nodes[[1]]
        prop.rule <- rules[[1]]
        temp <- cbind(xpred, 1:dim(xpred)[1])
        temp.t <- subset(temp, temp[,prop.pred] < xcut[[prop.pred]][prop.rule])
        prop.prob_star <- prop.prob
        prop.prob <- prop.prob_star[prop.pred] #P(selecting the jth attribute to split on)
        temp.L <- temp.t[,dim(temp.t)[2]]
        temp.R <- setdiff(temp[,dim(temp)[2]], temp.L)
        unique.len <- length(unique(temp[,prop.pred]))
      }else{
          #        prop.t.node <- prop.t.node[-1,]
        prop.pred <- i.nodes[[t.node[singly.inode,1]-1]][t.node[singly.inode,2]/2] # predictor of the singly internal node
        prop.rule <- rules[[t.node[singly.inode,1]-1]][t.node[singly.inode,2]/2] # rule of the singly internal node

        lv <- NULL # hierarchical level of the selected node
        lv[1] <- t.node[singly.inode,1]-1
        h.lv <- NULL   # location in a given hierarchical level
        h.lv[1] <- t.node[singly.inode,2]/2
        ii <- 1
        while(lv[ii]!=1){
            ii <- ii+1
            if(h.lv[ii-1] %% 2 ==0){
                h.lv[ii] <- h.lv[ii-1]/2
                lv[ii] <- t.node[singly.inode,1]-1-ii+1
            }else{
                h.lv[ii] <- (h.lv[ii-1]+1)/2
                lv[ii] <- t.node[singly.inode,1]-1-ii+1
            }
        }
        h.lv <- rev(h.lv)
        lv <- rev(lv)
        
        temp <- cbind(xpred, 1:dim(xpred)[1])
        if(h.lv[2] %% 2 ==0 ){
            temp.t <- subset(temp, temp[,i.nodes[[1]]] >= xcut[[i.nodes[[1]]]][rules[[1]]])
            temp.R <- temp.t[,dim(temp.t)[2]]
            temp.L <- setdiff(1:dim(temp)[1],temp.R)
        } else {
            temp.t <- subset(temp, temp[,i.nodes[[1]]] < xcut[[i.nodes[[1]]]][rules[[1]]])
            temp.L <- temp.t[,dim(temp.t)[2]]
            temp.R <- setdiff(1:dim(temp)[1],temp.L)
        }
        
        if(length(h.lv) > 2){
            for(i in 3:length(h.lv)){
                if(h.lv[i] %% 2 ==0){
                    temp <- temp.t
                    temp.t <- subset(temp, temp[,i.nodes[[lv[i-1]]][h.lv[i-1]]] >= xcut[[i.nodes[[lv[i-1]]][h.lv[i-1]]]][rules[[lv[i-1]]][h.lv[i-1]]])
                    temp.R <- temp.t[,dim(temp.t)[2]]
                    temp.L <- setdiff(temp[,dim(temp)[2]], temp.R)
                } else {
                    temp <- temp.t
                    temp.t <- subset(temp, temp[,i.nodes[[lv[i-1]]][h.lv[i-1]]] < xcut[[i.nodes[[lv[i-1]]][h.lv[i-1]]]][rules[[lv[i-1]]][h.lv[i-1]]])
                    temp.L <- temp.t[,dim(temp.t)[2]]
                    temp.R <- setdiff(temp[,dim(temp)[2]], temp.L)
                }
            }
        }
        enough.unique <- which(sapply(1:P, function(x) length(unique(temp.t[,x]))) >= 2) # covariates having enough unique values (2 or more)
        prop.prob_star <- prop.prob  # P(selecting the jth attribute to split on)
        prop.prob <- prop.prob_star[prop.pred]/sum(prop.prob_star[enough.unique])
        unique.len <- length(unique(temp.t[,prop.pred]))

        temp <- temp.t
        temp.t <- subset(temp, temp[,prop.pred] < xcut[[prop.pred]][prop.rule])
        temp.L <- temp.t[,dim(temp.t)[2]]
        temp.R <- setdiff(temp[,dim(temp)[2]], temp.L)
    }
    
    # transition
    transition_forward <- p.prune * 1/count
    transition_back <- p.grow * 1 / (dim(t.node)[1] - 1) * (prop.prob) * (1 / unique.len)
  
    # Transition ratio
    TRANS <- log(p.grow) - (log(dim(t.node)[1] - 1)) + log(max(prop.prob, 0)) - log(unique.len) - log(p.prune) + log(count)
  
    # Likelihood ratio
    nlL <- length(temp.L)
    nlR <- length(temp.R)
  
    LH <- log(sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu))/sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu))) + (sigma_mu / (2*sigma2) * ( -(sum(R[temp.L]))^2/(sigma2+nlL*sigma_mu) -(sum(R[temp.R]))^2/(sigma2+nlR*sigma_mu) + (sum(R[union(temp.R, temp.L)]))^2/(sigma2+(nlR+nlL)*sigma_mu) ) )
  
  
    # Structure ratio
    d <- t.node[singly.inode,1]-1-1
    STR <- -log(alpha) - 2*log((1- alpha / (2 + d)^beta )) + log((1 + d)^beta - alpha) - log(max(prop.prob, 0)) + log(unique.len)
  
    r <- TRANS+LH+STR
    if(r > log(runif(1))){
        i.nodes <- prop.i.nodes;
        rules <- prop.rules;
    }
    
    t.node <- c(0, 0)
    for(i in 2:length(i.nodes)){
        temp <- grep(999, i.nodes[[i]])
        if(length(temp) > 0){
            t.node <- rbind(t.node, cbind(i, temp))
        }
    }
    
    if(i.nodes[[1]]==999){
        i.nodes[[3]] <- NULL
        rules[[3]] <- NULL
    }else{
        t.node <- t.node[-1,]  # Level and Location of each terminal node
        if(t.node[dim(t.node)[1],1]<(length(i.nodes)-1)){
        i.nodes[[length(i.nodes)]] <- NULL
        rules[[length(rules)]] <- NULL
        }
    }
    
  return(list(i.nodes=i.nodes, rules=rules))
  
}
