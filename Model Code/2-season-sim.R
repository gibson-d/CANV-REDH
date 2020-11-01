#################################################################################################
# Simulate a basic 2-season (breeding and winter banding regime), 2-age class recovery data and model
###################################################################################################

n.occasions <- 20                    # Number of occasions
marked.jb <- rpois(n.occasions,1000)  # Annual number of newly marked young  (breeding)
marked.jw <- rpois(n.occasions,1000)  # Annual number of newly marked young  (winter)
marked.ab <- rpois(n.occasions,100)   # Annual number of newly marked adults (breeding)
marked.aw <- rpois(n.occasions,100)   # Annual number of newly marked adults (winter)

mean.s1  <- 0.45
mean.s1w <- 0.75
sjuv   <- mean.s1 * mean.s1w           # Juvenile survival probability

mean.s2 <- 0.65
mean.s2w <- 0.85
sadu   <- mean.s2 * mean.s2w           # Adult survival probability


rjuv <- 0.5                          # Juvenile recovery probability (direct)
rad  <- 0.45                          # Adult recovery probability (direct)
rind <- 0.20                          # Indirect recovery

##########################################
# survival rates
##########################################
sj.b <- matrix(0,n.occasions,n.occasions)
diag(sj.b)      <- mean.s1
sj.b[upper.tri(sj.b, diag = FALSE)] <- mean.s2
sj.b1 <- sj.b
diag(sj.b1) <- 1

sj.w <- matrix(0,n.occasions,n.occasions)
diag(sj.w)      <- mean.s1w
sj.w[upper.tri(sj.w, diag = FALSE)] <-mean.s2w 

sa.b <- matrix(0,n.occasions,n.occasions)
diag(sa.b)      <- mean.s2 
sa.b[upper.tri(sa.b, diag = FALSE)] <-mean.s2 
sa.b1 <- sa.b
diag(sa.b1) <- 1
sa.w <- matrix(0,n.occasions,n.occasions)
diag(sa.w)      <- mean.s2w
sa.w[upper.tri(sa.w, diag = FALSE)] <- mean.s2w


##########################################
# reporting rates
##########################################
rj.b <- matrix(0,n.occasions,n.occasions)
diag(rj.b)      <- rjuv 
rj.b[upper.tri(rj.b, diag = FALSE)] <-rind

rj.w <- matrix(0,n.occasions,n.occasions)
diag(rj.w)      <- 0
rj.w[upper.tri(rj.w, diag = FALSE)] <-rind 
for(t in 1:(nrow(rj.w)-1)){
  for(j in t+1){
    rj.w[t,j] <- rad
  }
}

ra.b <- matrix(0,n.occasions,n.occasions)
diag(ra.b)      <- rad
ra.b[upper.tri(ra.b, diag = FALSE)] <-rind

ra.w <- matrix(0,n.occasions,n.occasions)
diag(ra.w)      <- 0 
ra.w[upper.tri(ra.w, diag = FALSE)] <-rind
for(t in 1:(nrow(ra.w)-1)){
  for(j in t+1){
    ra.w[t,j] <- rad
  }
}


# Define matrices with survival and recovery probabilities
SJ1  <- sj.b[rep(seq_len(nrow(sj.b)),marked.jb),]
SJ2  <- sj.w[rep(seq_len(nrow(sj.w)),marked.jb),]

SJW1 <- sj.b1[rep(seq_len(nrow(sj.b1)),marked.jw),]
SJW2 <- sj.w[rep(seq_len(nrow(sj.w)),marked.jw),]

SA1  <- sa.b[rep(seq_len(nrow(sa.b)),marked.ab),]
SA2  <- sa.w[rep(seq_len(nrow(sa.w)),marked.ab),]

SAW1 <- sa.b1[rep(seq_len(nrow(sa.b1)),marked.aw),]
SAW2 <- sa.w[rep(seq_len(nrow(sa.w)),marked.aw),]

RJ  <- rj.b[rep(seq_len(nrow(rj.b)),marked.jb),]
RJW <- rj.w[rep(seq_len(nrow(rj.w)),marked.jw),]
RA  <- ra.b[rep(seq_len(nrow(ra.b)),marked.ab),]
RAW <- ra.w[rep(seq_len(nrow(ra.w)),marked.aw),]


# Create Capture Histories
simul.mrb <- function(SS, SW, R, marked){
  n.occasions <- dim(SS)[2]
  MR <- matrix(NA, ncol = n.occasions+1, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:(n.occasions), marked) 
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    MR[i, mark.occ[i]] <- 1    # Write an 1 at the release occasion
    
    for (t in (mark.occ[i]):n.occasions){
      # Bernoulli trial: has individual survived occasion? 
      sur1 <- rbinom(1, 1, SS[i,t])
      sur2 <- rbinom(1, 1, SW[i,t])
      if (sur1== 1 & sur2 == 1) next    # If still alive, move to next occasion 
      if (sur1== 1 & sur2 == 0) {
        MR[i,t+1] <- 0
        break
      }
      # Bernoulli trial: has dead individual been recovered? 
      rp <- rbinom(1, 1, R[i,t] )
      if (rp==0 ){
        MR[i,t+1] <- 0
        break
      } 
      if (rp==1){
        MR[i,t+1] <- 1
        break
      } 
    } #t
  } #i
  # Replace the NA in the file by 0
  MR[which(is.na(MR))] <- 0
  return(MR)
}

# Create Capture Histories
simul.mrw <- function(SS, SW, R, marked){
  n.occasions <- dim(SS)[2]
  MR <- matrix(NA, ncol = n.occasions+1, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:(n.occasions-1), marked) 
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    MR[i, mark.occ[i]] <- 1    # Write an 1 at the release occasion
    for (t in (mark.occ[i]):n.occasions){
      
      # Bernoulli trial: has individual survived occasion? 
      sur1 <- rbinom(1, 1, SS[i,t])
      sur2 <- rbinom(1, 1, SW[i,t])
      if (sur1== 1 & sur2 == 1) next    # If still alive, move to next occasion 
      if (sur1== 1 & sur2 == 0) {
        MR[i,t+1] <- 0
        break
      }
      # Bernoulli trial: has dead individual been recovered? 
      rp <- rbinom(1, 1, R[i,t] )
      if (rp==0 ){
        MR[i,t+1] <- 0
        break
      } 
      if (rp==1){
        MR[i,t+1] <- 1
        break
      } 
    } #t
  } #i
  # Replace the NA in the file by 0
  MR[which(is.na(MR))] <- 0
  return(MR)
}
# Execute simulation function
MRjb <- simul.mrb(SJ1,SJ2, RJ, marked.jb)
MRab <- simul.mrb(SA1,SA2, RA, marked.ab)

MRjw <- simul.mrw(SJW1, SJW2, RJW, marked.jw[-n.occasions])
MRaw <- simul.mrw(SAW1, SAW2, RAW, marked.aw[-n.occasions])


marray.dead <- function(MR){
  nind <- dim(MR)[1]
  n.occasions <- dim(MR)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Create vector with occasion of marking 
  get.first <- function(x) min(which(x!=0))
  f <- apply(MR, 1, get.first)
  # Calculate the number of released individuals at each time period
  first <- as.numeric(table(f))
  for (t in 1:n.occasions){
    m.array[t,1] <- first[t]
  }
  # Fill m-array with recovered individuals
  rec.ind <- which(apply(MR, 1, sum)==2)
  rec <- numeric()
  for (i in 1:length(rec.ind)){
    d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
    rec[i] <- d + f[rec.ind[i]]
    m.array[f[rec.ind[i]],rec[i]] <- m.array[f[rec.ind[i]],rec[i]] + 1
  }
  # Calculate the number of individuals that are never recovered
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}
# Summarize data in m-arrays

marr.jb <- marray.dead( MRjb) 
marr.jw <- marray.dead( MRjw) 
marr.ab <- marray.dead( MRab) 
marr.aw <- marray.dead( MRaw) 

marr.main <- abind(marr.jb,marr.ab, along = 3)
marr.wint <- abind(marr.jw,marr.aw, along = 3)

marr.wint[n.occasions,n.occasions+1,] <- 0

rel.main <- cbind(rowSums(marr.jb),rowSums(marr.ab))
rel.wint <- cbind(rowSums(marr.jw),rowSums(marr.aw))


rel.wint[n.occasions,]<- 0

# Specify model in BUGS language

sink("2-season-2-age-rec.jags")
cat("
model {

# Model Priors
  mean.s1b ~ dbeta(1,1) # Feb-Aug (Juv)
  mean.s1w ~ dbeta(1,1) # Sept-Jan (Juv)
  mean.s2b ~ dbeta(1,1) # Feb-Aug (Adu)
  mean.s2w ~ dbeta(1,1) # Sept-Jan (Adu)
    mean.r ~ dbeta(1,1) # indirect recovery
  mean.rda ~ dbeta(1,1) # direct adult
  mean.rdj ~ dbeta(1,1) # direct juvenile
  
    for (t in 1:n.occasions){
    S1B[t] <- mean.s1b
    S1W[t] <- mean.s1w
    S2B[t] <- mean.s2b
    S2W[t] <- mean.s2w
      r[t] <- mean.r
   rd[1,t] <- mean.rdj
   rd[2,t] <- mean.rda
    }
    
# Model Likelihood
     for(k in 1:n.age){
         for (t in 1:n.occasions){
          marr[t,1:(n.occasions+1),k] ~ dmulti(    pr[k,t,1:(n.occasions+1)], rel[t,k])                              # Breeding Releases
     }
     for (t in 1:(n.occasions-1)){
      marr.win[t,1:(n.occasions+1),k] ~ dmulti(pr.win[k,t,1:(n.occasions+1)], rel.win[t,k])                          # Winter Releases
    }
  }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
# Define the cell probabilities of the winter banding m-array
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SY Birds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Main diagonal
      for (t in 1:(n.occasions-1)){
        pr.win[1,t,t] <- 0
        # Directly above diagonal
        pr.win[1,t,t+1] <- S1B[t] * (1-S2W[t+1]) * rd[2,t+1]  
        # Above main diagonal
        for (i in (t+2):n.occasions){
          pr.win[1,t,i] <-  S1B[t] * prod(S2W[(t+1):(i-1)]) *  prod(S2B[(t+1):(i-1)]) * (1-S2W[i]) * r[i] 
        } #i
        # Below main diagonal
        for (i in 1:(t-1)){
          pr.win[1,t,i] <- 0
        } #i
      } #t
      # Last column: probability of non-recovery
      for (t in 1:(n.occasions-1)){
        pr.win[1,t,n.occasions+1] <- 1-sum(pr.win[1,t,1:n.occasions])
      } #t
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ASY Birds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      for (t in 1:(n.occasions-1)){
        pr.win[2,t,t] <- 0
        # Directly above diagonal
        pr.win[2,t,t+1] <- S2B[t] * (1-S2W[t+1]) * rd[2,t+1]  
        # Above main diagonal
        for (i in (t+2):n.occasions){
          pr.win[2,t,i] <-   prod(S2B[t:(i-1)]) *  prod(S2W[(t+1):(i-1)]) * (1-S2W[i]) * r[i] 
        } #i
        # Below main diagonal
        for (i in 1:(t-1)){
          pr.win[2,t,i] <- 0
        } #i
      } #t
      # Last column: probability of non-recovery
      for (t in 1:(n.occasions-1)){
        pr.win[2,t,n.occasions+1] <- 1-sum(pr.win[2,t,1:n.occasions])
      } #t
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  # Define the cell probabilities of the breeding banding m-array
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # HY Birds
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Main diagonal
          #One above main diagonal
      for (t in 1:(n.occasions-1)){
        for (i in t+1){
          pr[1,t,i] <- S1W[t] * S1B[t] * (1-S2W[i])* r[i]
        } #i
      }#t
      # Main diaganol
      for (t in 1:n.occasions){
        pr[1,t,t] <- (1- S1W[t])* rd[1,t] 
         for (i in (t+2):n.occasions){
          pr[1,t,i] <- S1W[t] * S1B[t]  * prod(S2W[(t+1):(i-1)]) * prod(S2B[(t+1):(i-1)])  * (1-S2W[i])* r[i]
        } #i
        # Below main diagonal
        for (i in 1:(t-1)){
          pr[1,t,i] <- 0
        } #i
      } #t
      # Last column: probability of non-recovery
      for (t in 1:n.occasions){
        pr[1,t,n.occasions+1] <- 1-sum(pr[1,t,1:n.occasions])
      } #t
      
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # AHY Birds
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   
      # Main diagonal
      for (t in 1:n.occasions){
        pr[2,t,t] <- (1-S2W[t])* rd[2,t]
        # Above main diagonal
        for (i in (t+1):n.occasions){
          pr[2,t,i] <-  prod(S2W[t:(i-1)]) * prod(S2B[t:(i-1)]) * (1-S2W[i]) * r[i]
        } #i
        # Below main diagonal
        for (i in 1:(t-1)){
          pr[2,t,i] <- 0
        } #i
      } #t
      # Last column: probability of non-recovery
      for (t in 1:n.occasions){
        pr[2,t,n.occasions+1] <- 1-sum(pr[2,t,1:n.occasions])
      } #t
    
  


}
",fill = TRUE)
sink()


# Bundle data
jags.data <- list(marr = marr.main, marr.win = marr.wint, n.age = 2,
                  rel = rel.main, rel.win = rel.wint,
                  n.occasions = dim(marr.main)[2]-1)

# Initial values
inits <- function(){list(  mean.s1b= rbeta(1,1,1),
                           mean.s1w= rbeta(1,1,1),
                           mean.s2b= rbeta(1,1,1),
                           mean.s2w= rbeta(1,1,1),
                           mean.r  = rbeta(1,1,1),
                           mean.rda = rbeta(1,1,1),
                           mean.rdj = rbeta(1,1,1))}

# Parameters monitored
parameters <- c("mean.s1b", "mean.s1w",'mean.s2b','mean.s2w','mean.r','mean.rda', 'mean.rdj')

# MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3

# Call JAGS from R (BRT <1 min)
mr <- jags(jags.data, inits, parameters, "2-season-2-age-rec.jags",
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = 5000,
           n.cores = 3, parallel = TRUE)




refer <- c(mean.s1w,
           mean.s1,
           mean.s2w,
           mean.s2,
           rind, 
           rad,
           rjuv,NA)

library(MCMCvis)
MCMCtrace(mr, gvals = refer)


