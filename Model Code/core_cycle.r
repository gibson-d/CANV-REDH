
era <- as.numeric(as.factor(ordering.b$Era))
period <- c(rep(1,43),rep(2,14), rep(3,12))

st <- c(1,44,58)
start <- st[era]

rel.time1 <- as.numeric(as.factor(ordering.b$Month))
sex1      <- as.numeric(as.factor(ordering.b$sex))
species1  <- as.numeric(as.factor(ordering.b$species))
age1      <- as.numeric(as.factor(ordering.b$age))

rel.time2 <- as.numeric(as.factor(ordering.w$Month))
sex2      <- as.numeric(as.factor(ordering.w$sex))
species2  <- as.numeric(as.factor(ordering.w$species))
age2      <- as.numeric(as.factor(ordering.w$age))


order <- matrix(c(1,1,1,1,2,1,1,1,1,2,1,1,2,2,1,1,1,1,2,1,2,1,2,1,1,2,2,1,2,2,2,1,1,1,1,2,2,1,1,2,1,2,1,2,2,2,1,2,1,1,2,2,2,1,2,2,1,2,2,2,2,2,2,2), ncol = 4, byrow = TRUE)


marr.main[5,70,17] <- 0
 marr.main[5,5,17] <- 0


marr.off[34,34,16] <- 0
    rel.off[34,16] <- sum(marr.off[34,,16])


###############################################
model_1 <- nimbleCode( {
# Model Priors
  for(k in 1:n.species){
    for(l in 1:n.sex){
      for(j in 1:n.era){
            mean.r[j,k,l] ~ dbeta(1,1) # indirect recovery
              mu.r[j,k,l] <- logit(mean.r[j,k,l])
      }
          
               local[k,l] ~ dbeta(1,1) # Monthly pre-harvest season survival (July-August)
      
            mean.s1w[k,l] ~ dbeta(1,1) # Monthly HY  overwinter survival (August - Dec)
            mean.s2w[k,l] ~ dbeta(1,1) # Monthly AHY overwinter survival (August - Dec)
            
            mean.s1b[k,l] ~ dbeta(1,1) # Monthly SY  oversummer survival (Jan-July)
            mean.s2b[k,l] ~ dbeta(1,1) # Monthly ASY oversummer survival (Jan-July)
            
              mu.s1w[k,l] <- logit(mean.s1w[k,l])
              mu.s2w[k,l] <- logit(mean.s2w[k,l])
              
              mu.s1b[k,l] <- logit(mean.s1b[k,l])
              mu.s2b[k,l] <- logit(mean.s2b[k,l])
              
                HY[k,l]   <- mean.s1b[k,l]^6 * mean.s1w[k,l]^4 *  local[k,l]^2
               AHY[k,l]   <- mean.s2b[k,l]^6 * mean.s2w[k,l]^6
               
        for(i in 1:n.era){       
          for(j in 1:n.age){
           mean.rd[i,j,k,l] ~ dbeta(1,1) # direct recoveries 
             mu.rd[i,j,k,l] <- logit(mean.rd[i,j,k,l])
       }
     }

    for (t in 1:n.occasions){
      
             logit(S1Ba[k,l,t]) <- mu.s1b[k,l] + eps[1,1,k,l,t]
             logit(S1Wa[k,l,t]) <- mu.s1w[k,l] + eps[2,1,k,l,t]
             
             logit(S2Ba[k,l,t]) <- mu.s2b[k,l] + eps[1,2,k,l,t]
             logit(S2Wa[k,l,t]) <- mu.s2w[k,l] + eps[2,2,k,l,t]
        
                S1W[1,1,k,l,t] <- pow(S1Wa[k,l,t],4) * pow(local[k,l],2)   # Locals released by end of the July
                S1W[2,1,k,l,t] <- pow(S1Wa[k,l,t],4) * local[k,l]          # Locals released by end of the August
             
                S1W[1,2,k,l,t] <- pow(S1Wa[k,l,t],5)   # HY released by end of the August
                S1W[2,2,k,l,t] <- pow(S1Wa[k,l,t],4)   # HY released by end of the September
              
                  S2W[1,k,l,t] <- pow(S2Wa[k,l,t],6)   # AHY released by end of the July
                  S2W[2,k,l,t] <- pow(S2Wa[k,l,t],5)   # AHY released by end of the August
              
                  S1B[1,k,l,t] <- pow(S1Ba[k,l,t],6)   # SY released by end of the Jan
                  S1B[2,k,l,t] <- pow(S1Ba[k,l,t],5)   # SY released by end of the Feb
                  S1B[3,k,l,t] <- pow(S1Ba[k,l,t],4)   # SY released by end of the Mar
            
                  S2B[1,k,l,t] <- pow(S2Ba[k,l,t],6)   # ASY released by end of the Jan
                  S2B[2,k,l,t] <- pow(S2Ba[k,l,t],5)   # ASY released by end of the Feb
                  S2B[3,k,l,t] <- pow(S2Ba[k,l,t],4)   # ASY released by end of the Mar
               
              HYt[k,l,t] <- pow(S1Ba[k,l,t],6) * pow(S1Wa[k,l,t],4) * pow(local[k,l],2)
             AHYt[k,l,t] <- pow(S2Ba[k,l,t],6) * pow(S2Wa[k,l,t],6)
          
       logit(rd[1,k,l,t]) <- mu.rd[period[t], 1,k,l] + eps.r[1,k,l,t]
       logit(rd[2,k,l,t]) <- mu.rd[period[t], 2,k,l] + eps.r[2,k,l,t]
       logit(  r[k,l,t])  <-  mu.r[period[t],k,l]    + eps.r[2,k,l,t]
      
    }
  }
}       
                # Covariance matrices for survival and harvest rates
                for(i in 1:n.species){
                   mu.harv[i] <- 0
                       sig[i] ~ T(dt(0, pow(2.5, -2), 1),0,)
                  sig.harv[i] ~ T(dt(0, pow(2.5, -2), 1),0,)
                     rho[i,i] <- 1
                rho.harv[i,i] <- 1
                
                # Priors for cyclic-trend
                  nu[i] <- 0    #Intercept
                   K[i] ~ dgamma(20,1)              #Periodicity (informed prior)
                  for(j in 1:2){
                    beta[j,i] ~ dnorm(0, 0.001)    # Combined arc-phase and amplitude
                  }
                }  
                for (i in 2){
                  for (j in 1:(i - 1)) {
                        rho[i, j] ~ dunif(-1, 1)
                        rho[j, i] <-rho[i, j]
                    rho.harv[i,j]~ dunif(-1, 1)
                    rho.harv[j,i] <- rho.harv[i, j]
                  }
                }
  
                for (i in 1:(n.species)){  
                  for (j in 1:(n.species)){
                    Sigma.harv[i, j] <-  rho.harv[i,j] *  sig.harv[i] *  sig.harv[j]
                         Sigma[i, j] <-       rho[i,j] *  sig[i] *  sig[j]
                  }
                }
  
    for (t in 1:n.occasions){
      for (i in 1:(n.species)){  
           # Cyclic Trend
            mu[i,t] <- nu[i] + beta[1,i] * cos((2 *  3.141593 * t)/K[i]) + beta[2,i] * sin((2 * 3.141593 * t)/K[i])
        }
      
      # Variation around trend
      EPS[t,1:(n.species)] ~ dmnorm(mu[1:(n.species),t],   cov = Sigma[1:(n.species),1:(n.species)])
      EPS.harv[t,1:(n.species)] ~ dmnorm(mu.harv[1:(n.species)],  cov = Sigma.harv[1:(n.species),1:(n.species)])
      
      # Variance inflation factor on three of the four 'groups'. Similiar to Peron and Koons (2012)
      # Not certain this is the best approach to deal with temporal variation between sexes/ages
      for(i in 1:2){
        for(j in 1:3){
            vif[i,j] ~ dnorm(1, 0.01)
          vif.h[i,j] ~ dnorm(1, 0.01)
        }
      }
      
     # Age - Species - Sex - Time
      eps.r[1,1,1,t] <-  EPS.harv[t,1] 
      eps.r[2,1,1,t] <-  EPS.harv[t,1] * vif.h[1,1] 
      eps.r[1,1,2,t] <-  EPS.harv[t,1] * vif.h[1,2] 
      eps.r[2,1,2,t] <-  EPS.harv[t,1] * vif.h[1,3] 
      
      eps.r[1,2,1,t] <-  EPS.harv[t,2] 
      eps.r[2,2,1,t] <-  EPS.harv[t,2] * vif.h[2,1]
      eps.r[1,2,2,t] <-  EPS.harv[t,2] * vif.h[2,2]
      eps.r[2,2,2,t] <-  EPS.harv[t,2] * vif.h[2,3]
      
      # Season - Age - Species - Sex - Time
      # HY - F - CANV
      eps[1,1,1,1,t] <- EPS[t,1]  
      eps[2,1,1,1,t] <- EPS[t,1]  
      # AHY - F - CANV
      eps[1,2,1,1,t] <- EPS[t,1] * vif[1,1] 
      eps[2,2,1,1,t] <- EPS[t,1] * vif[1,1]  
      # HY - M - CANV
      eps[1,1,1,2,t] <- EPS[t,1] * vif[1,2]   
      eps[2,1,1,2,t] <- EPS[t,1] * vif[1,2]  
      # AHY - M - CANV
      eps[1,2,1,2,t] <- EPS[t,1] * vif[1,3] 
      eps[2,2,1,2,t] <- EPS[t,1] * vif[1,3] 
      # HY - F - REDH
      eps[1,1,2,1,t] <- EPS[t,2] 
      eps[2,1,2,1,t] <- EPS[t,2] 
      # AHY - F - REDH
      eps[1,2,2,1,t] <- EPS[t,2] * vif[2,1] 
      eps[2,2,2,1,t] <- EPS[t,2] * vif[2,1] 
      # HY - M - REDH
      eps[1,1,2,2,t] <- EPS[t,2] * vif[2,2] 
      eps[2,1,2,2,t] <- EPS[t,2] * vif[2,2] 
      # AHY - M - RED 
      eps[1,2,2,2,t] <- EPS[t,2] * vif[2,3] 
      eps[2,2,2,2,t] <- EPS[t,2] * vif[2,3] 
    }
                
  # Model Likelihood
  for(k in 1:24){
    for (t in 1:n.occasions){
      marr[t, 1:(n.occasions+1),k]    ~ dmulti(pr[rel.time1[k],age1[k],species1[k],sex1[k],t,1:(n.occasions+1)], rel[t,k])   
    }
    for (t in 1:(n.occasions-1)){
      marr.win[t,1:(n.occasions+1),k] ~ dmulti(pr.win[rel.time2[k], age2[k],species2[k],sex2[k],t,1:(n.occasions+1)], rel.win[t,k])                        
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  # Define the cell probabilities of the winter banding m-array                    #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # SY Birds                                 #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Main diagonal
  for(k in 1:n.species){
    for(l in 1:n.sex){
      for(j in 1:n.rel2){
        for (t in 1:(n.occasions-1)){
          pr.win[j,1,k,l,t,t] <- 0
          # Directly above diagonal
          pr.win[j,1,k,l,t,t+1] <- S1B[j,k,l,t] * (1-S2W[1,k,l,t+1]) * rd[2,k,l,t+1]  
          # Above main diagonal
          for (i in (t+2):n.occasions){
            pr.win[j,1,k,l,t,i] <-  S1B[j,k,l,t] * prod(S2W[1,k,l,(t+1):(i-1)]) *  prod(S2B[1,k,l,(t+1):(i-1)]) * (1-S2W[1,k,l,i]) * r[k,l,i] 
          } #i
          # Below main diagonal
          for (i in 1:(t-1)){
            pr.win[j,1,k,l,t,i] <- 0
          } #i
        } #t
        # Last column: probability of non-recovery
        for (t in 1:(n.occasions-1)){
          pr.win[j,1,k,l,t,n.occasions+1] <- 1-sum(pr.win[j,1,k,l,t,1:n.occasions])
        } #t
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # ASY Birds
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        for (t in 1:(n.occasions-1)){
          pr.win[j,2,k,l,t,t] <- 0
          # Directly above diagonal
          pr.win[j,2,k,l,t,t+1] <- S2B[j,k,l,t] * (1-S2W[1,k,l,t+1]) * rd[2,k,l,t+1]  
          # Above main diagonal
          for (i in (t+2):n.occasions){
            pr.win[j,2,k,l,t,i] <- S2B[j,k,l,t] * prod(S2B[1,k,l,(t+1):(i-1)]) *  prod(S2W[1,k,l,(t+1):(i-1)]) * (1-S2W[1, k,l,i]) * r[k,l,i] 
          } #i
          # Below main diagonal
          for (i in 1:(t-1)){
            pr.win[j,2,k,l,t,i] <- 0
          } #i
        } #t
        # Last column: probability of non-recovery
        for (t in 1:(n.occasions-1)){
          pr.win[j,2,k,l,t,n.occasions+1] <- 1-sum(pr.win[j,2,k,l,t,1:n.occasions])
        } #t
      } # j
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
      # Define the cell probabilities of the breeding banding m-array
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      for(j in 1:n.rel1){  
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # LOCAL Birds
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #One above main diagonal
        for (t in 1:(n.occasions-1)){
          for (i in t+1){
            pr[j,1,k,l,t,i] <- S1W[j,1,k,l,t] * S1B[1,k,l,t] * (1-S2W[1,k,l,i])* r[k,l,i]
          } #i
        }#t
        # Main diagonal
        for (t in 1:n.occasions){
          pr[j,1,k,l,t,t] <- (1- S1W[j,1,k,l,t])* rd[1,k,l,t] 
          for (i in (t+2):n.occasions){
            pr[j,1,k,l,t,i] <- S1W[j,1,k,l,t] * S1B[1,k,l,t]  * prod(S2W[1,k,l,(t+1):(i-1)]) * prod(S2B[1,k,l,(t+1):(i-1)])  * (1-S2W[1,k,l,i])* r[k,l,i]
          } #i
          # Below main diagonal
          for (i in 1:(t-1)){
            pr[j,1,k,l,t,i] <- 0
          } #i
        } #t
        # Last column: probability of non-recovery
        for (t in 1:n.occasions){
          pr[j,1,k,l,t,n.occasions+1] <- 1-sum(pr[j,1,k,l,t,1:n.occasions])
        } #t
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # HY Birds
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #One above main diagonal
        for (t in 1:(n.occasions-1)){
          for (i in t+1){
            pr[j,2,k,l,t,i] <- S1W[j,2,k,l,t] * S1B[1,k,l,t] * (1-S2W[1,k,l,i])* r[k,l,i]
          } #i
        }#t
        # Main diagonal
        for (t in 1:n.occasions){
          pr[j,2,k,l,t,t] <- (1- S1W[j,2,k,l,t])* rd[1,k,l,t] 
          for (i in (t+2):n.occasions){
            pr[j,2,k,l,t,i] <- S1W[j,2,k,l,t] * S1B[1,k,l,t]  * prod(S2W[1,k,l,(t+1):(i-1)]) * prod(S2B[1,k,l,(t+1):(i-1)])  * (1-S2W[1,k,l,i])* r[k,l,i]
          } #i
          # Below main diagonal
          for (i in 1:(t-1)){
            pr[j,2,k,l,t,i] <- 0
          } #i
        } #t
        # Last column: probability of non-recovery
        for (t in 1:n.occasions){
          pr[j,2,k,l,t,n.occasions+1] <- 1-sum(pr[j,2,k,l,t,1:n.occasions])
        } #t
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # AHY Birds
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # Main diagonal
        for (t in 1:n.occasions){
          pr[j,3,k,l,t,t] <- (1-S2W[j,k,l,t])* rd[2,k,l,t]
          # Above main diagonal
          for (i in (t+1):n.occasions){
            pr[j,3,k,l,t,i] <- S2W[j,k,l,t] *  prod(S2W[1,k,l,(t+1):(i-1)]) * prod(S2B[1,k,l,t:(i-1)]) * (1-S2W[1,k,l,i]) * r[k,l,i]
          } #i
          # Below main diagonal
          for (i in 1:(t-1)){
            pr[j,3,k,l,t,i] <- 0
          } #i
        } #t
        # Last column: probability of non-recovery
        for (t in 1:n.occasions){
          pr[j,3,k,l,t,n.occasions+1] <- 1-sum(pr[j,3,k,l,t,1:n.occasions])
        } #t
        
      }
    }
  }
  
})

##############################################    
# Current glitches in the data formatting    
marr.main[5,70,17] <- 0
marr.main[5,5,17]  <- 0

marr.off[34,34,16] <- 0
rel.off[34,16]     <- sum(marr.off[34,,16])
###############################################



# Bundle data
dat <- list(marr = marr.main, marr.win = marr.off, 
             rel = rel.main, rel.win = rel.off)


constants <- list(n.age = 2, n.species = 2, n.sex = 2, n.era = 3, n.season = 2, 
                  n.rel1 = 2, n.rel2 = 3,
                  n.occasions = dim(marr.main)[2]-1,
                  age1 = age1, age2 = age2,
                  species1 = species1, species2 = species2,
                  sex1 = sex1,   sex2 = sex2,
                  rel.time1 = rel.time1, rel.time2 = rel.time2,
                  period = period
                 
                  )

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

SigA <-  Posdef(n=2, ev=1:2)
rho  <- cov2cor(SigA)

SigB <-  Posdef(n=2, ev=1:2)
rhoh  <- cov2cor(SigB)

     


# Initial values
initsFunction <- function()list(  mean.s1b =  matrix(rbeta(4,5,2),2,2),
                                  mean.s1w =  matrix(rbeta(4,5,2),2,2),
                                  mean.s2b =  matrix(rbeta(4,5,2),2,2),
                                  mean.s2w =  matrix(rbeta(4,5,2),2,2),
                                  mean.r   =  array(rbeta(12,2,5), dim = c(3,2,2)),
                                  mean.rd  =  array(rbeta(24,2,5), dim = c(3,2,2,2)),
                                  rho      = rho, 
                                  rho.harv = rhoh, 
                                  vif = matrix(rnorm(6,1,1), 2, 3),
                                  vif.h = matrix(rnorm(6,1,1), 2, 3),
                                  sig = runif(2,0,1),
                                  sig.harv = runif(2,0,1), 
                                  
                                  local    = matrix(rbeta(4,5,2),2,2),
                                  beta = matrix(rnorm(4,1),2,2),
                                  K = c(15,15)
                                 )
                     
# Parameters monitored
pars <- c("mean.s1b", "mean.s1w",'mean.s2b','mean.s2w','mean.r','mean.rd','HY','AHY')

inits <- initsFunction() 


library(parallel)
library(coda)
set.seed(2)
nc <- 4 # number of chains

clduck<-makeCluster(nc,timeout=5184000)

clusterExport(clduck, c("model_1", "inits", "dat", "constants", "pars"))

for (j in seq_along(clduck)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(clduck[j], "inits")
}

out <- clusterEvalQ(clduck, {
  library(nimble)
  library(coda)
  model <- nimbleModel( code = model_1, constants = constants,  dat =  dat, inits = inits)
        
  model$simulate(c('S1Ba', 'S1Wa', 'S2Ba', 'S2Wa','local', 'S1W', 'S2W', 'S1B', 'S2B', 'HYt',  'AHYt', 
                   'rd', 'r', 'EPS', 'EPS.harv','eps.r', 'eps'))

  model$simulate(c('pr','pr.win'))
  
  modelConf  <- configureMCMC(model,
                              useConjugacy = FALSE,
                              monitors2 = c('AHY','HY','S1B','S1W','S2B','S2W', 'r','rd','K','beta','mu','vif','vif.h',
                                                  'AHYt','HYt','rho', 'local','EPS', 'EPS.harv'))
  modelMCMC  <- buildMCMC( modelConf  )
  
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  out1 <- runMCMC(CmodelMCMC,nburnin = 25000, niter = 75000, thin = 2)
  
  return(as.mcmc(out1))
})


# Build the mcmc.list
samples <- list( chain1 =  out[[1]]$samples, 
                 chain2 =  out[[2]]$samples, 
                 chain3 =  out[[3]]$samples,
                 chain4 =  out[[4]]$samples)

# prior-based parameters
mcmcList <- as.mcmc.list(lapply(samples, mcmc))

samples2 <- list( chain1 =  out[[1]]$samples2, 
                  chain2 =  out[[2]]$samples2, 
                  chain3 =  out[[3]]$samples2,
                  chain4 =  out[[4]]$samples2)
# all monitored parameters
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))


results <- list(model_1,mcmcList,mcmcList2 )

save(results, file = 'corrs_time-vars_11_17.rdata')


ages <- unique(ordering.w$age)
specs <- unique(ordering.b$species)
sexes <- unique(ordering.b$sex)
years <- c(1951:2019)
dr <- MCMCsummary(mcmcList2, 'rd')

dr.grid <- expand.grid(age = ages, species = specs, sex = sexes, year = years )
dr.add <- cbind.data.frame(dr, dr.grid)


dodge <- position_dodge(.5)
dir.rs <- ggplot(dr.add, aes( x = year, y = mean)) +
          geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), position = dodge) +
          facet_grid(sex ~ species+age) +
          coord_cartesian(ylim = c(0,1)) +
          scale_color_manual(values = c('dodgerblue', 'goldenrod','firebrick'))+ 
          theme_bw() + theme(legend.position = 'top')
dir.rs 


ir <- MCMCsummary(mcmcList2, 'r')
ir.grid <- expand.grid(species = specs, sex = sexes, year = years )
ir.add <- cbind.data.frame(ir, ir.grid)


dodge <- position_dodge(.75)
indir.rs <- ggplot(ir.add, aes( x = year, y = mean)) +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), position = dodge) +
  facet_grid(sex ~ species) +
  scale_color_manual(values = c('dodgerblue', 'goldenrod','firebrick'))+ 
  labs(y = 'Reporting Rate (Overwinter)')+
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() + theme(legend.position = 'top')
indir.rs


AHYs <- MCMCsummary(mcmcList2, 'AHYt')
AHS.grid <- expand.grid(species = specs, sex = sexes, year = years )
AHS.add <- cbind.data.frame(AHYs , AHS.grid)

AHYm <- MCMCsummary(mcmcList2, 'AHY')

AHYm$species <- rep(specs, 2)
AHYm$sex <- rep(sexes, each = 2)

dodge <- position_dodge(.25)
AH.Surv <- ggplot(AHS.add, aes( x = year, y = mean, color = species)) +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), position = dodge) +
  facet_wrap(~ sex, ncol = 1) +
  scale_color_manual(values = c('dodgerblue', 'goldenrod'))+ 
  geom_hline(data = AHYm, aes(yintercept = mean, color = species )) +
  labs(y = 'Annual adult survival', x = '')+
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() + theme(legend.position = 'top')
AH.Surv

HYs <- MCMCsummary(mcmcList2, 'HYt')
HS.grid <- expand.grid(species = specs, sex = sexes, year = years )
HS.add <- cbind.data.frame(HYs , HS.grid)

HYm <- MCMCsummary(mcmcList2, 'HY')

HYm$species <- rep(specs, 2)
HYm$sex <- rep(sexes, each = 2)

dodge <- position_dodge(.25)
H.Surv <- ggplot(HS.add, aes( x = year, y = mean, color = sex)) +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), position = dodge) +
  facet_wrap(~ species, ncol = 2) +
  scale_color_manual(values = c('dodgerblue', 'goldenrod'))+ 
  geom_hline(data = HYm, aes(yintercept = mean, color = sex )) +
  labs(y = 'Annual juvenile survival', x = '')+
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() + theme(legend.position = 'top')
H.Surv

S1B <- MCMCsummary(mcmcList2, 'S1B')
S2B <- MCMCsummary(mcmcList2, 'S2B')
S1W <- MCMCsummary(mcmcList2, 'S1W')
S2W <- MCMCsummary(mcmcList2, 'S2W')


HS1.grid <- expand.grid(release = c('Jan','Feb','Mar'), species = specs, sex = sexes, year = years )
HS2.grid <- expand.grid(release = c('Early','Late'), species = specs, sex = sexes, year = years )

S1B.add <- cbind.data.frame(S1B , HS1.grid)
S2B.add <- cbind.data.frame(S2B , HS1.grid)

S1W.add <- cbind.data.frame(S1W , HS2.grid)
S2W.add <- cbind.data.frame(S2W , HS2.grid)

S1B.add$Season <- rep('Breeding')
S2B.add$Season <- rep('Breeding')
S1W.add$Season <- rep('Winter')
S2W.add$Season <- rep('Winter')

S1 <- rbind.data.frame(S1B.add, S1W.add)
S2 <- rbind.data.frame(S2B.add, S2W.add)

HS.Surv <- ggplot(S1, aes( x = year, y = mean, color = sex)) +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), position = dodge) +
  facet_grid(Season + release~ species ) +
  scale_color_manual(values = c('dodgerblue', 'goldenrod'))+ 
  labs(y = 'Seasonal juvenile survival', x = '')+
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() + theme(legend.position = 'top')
HS.Surv

AHS.Surv <- ggplot(S2, aes( x = year, y = mean, color = sex)) +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), position = dodge) +
  facet_grid(Season + release~ species) +
  scale_color_manual(values = c('dodgerblue', 'goldenrod'))+ 
  labs(y = 'Seasonal adult survival', x = '')+
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() + theme(legend.position = 'top')
AHS.Surv

library(ggplot2)


trend <- MCMCsummary(mcmcList2, 'mu')
trend.grid <- expand.grid(species = specs, year = years )
trend.add <- cbind.data.frame(trend , trend.grid)


mu.Surv <- ggplot(trend.add , aes( x = year, y = mean, color = species)) +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), position = dodge) +
  scale_color_manual(values = c('dodgerblue', 'goldenrod'))+ 
  labs(y = 'Species trend in survival', x = '')+
  theme_bw() + theme(legend.position = 'top')
mu.Surv



