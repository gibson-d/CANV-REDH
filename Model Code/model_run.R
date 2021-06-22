library(nimble)
# Specify model in NIMBLE language
model_1 <- nimbleCode( {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Dead-Recovery model
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Priors and constraints for variance-covariance parameters
  
  rho[1,1] <- 1                # W_C vs W_C
  rho[1,2] ~ dunif(-1,1)       # W_C vs W_R
  rho[1,3] ~ dunif(-1,1)       # W_C vs H_C
  rho[1,4] <- 0                # W_C vs H_R
  rho[1,5] <- 0                # W_C vs S_C
  rho[1,6] <- 0                # W_C vs S_R
  
  rho[2,1] <- rho[1,2]         # W_R vs W_C 
  rho[2,2] <- 1                # W_R vs W_R
  rho[2,3] <- 0                # W_R vs H_C
  rho[2,4] ~ dunif(-1,1)       # W_R vs H_R
  rho[2,5] <- 0                # W_R vs S_C
  rho[2,6] <- 0                # W_R vs S_R
  
  rho[3,1] <- rho[1,3]         # H_C vs W_C
  rho[3,2] <- rho[2,3]         # H_C vs W_R
  rho[3,3] <- 1                # H_C vs H_C
  rho[3,4]  ~ dunif(-1,0.95)      # H_C vs H_R
  rho[3,5]  <- 0               # H_C vs S_C
  rho[3,6]  <- 0               # H_C vs S_R
  
  rho[4,1] <- rho[1,4]
  rho[4,2] <- rho[2,4]
  rho[4,3] <- rho[3,4]
  rho[4,4] <- 1
  rho[4,5] <- 0
  rho[4,6] <- 0
  
  rho[5,1] <-  rho[1,5]        # S_C vs W_C
  rho[5,2] <-  rho[2,5]        # S_C vs W_R
  rho[5,3] <-  rho[3,5]        # S_C vs H_C
  rho[5,4] <-  rho[4,5]        # S_C vs H_R
  rho[5,5] <- 1                # S_C vs S_C
  rho[5,6]  ~ dunif(-1,1)      # S_C vs S_R
  
  rho[6,1] <- rho[1,6]
  rho[6,2] <- rho[2,6]
  rho[6,3] <- rho[3,6]
  rho[6,4] <- rho[4,6]
  rho[6,5] <- rho[5,6]
  rho[6,6] <- 1
  
  for(i in 1:6){
    mu.dem[i] <- 0
    sig.dem[i] ~ T(dt(0, pow(2.5, -2), 1),0,2.5)
  }
  for(l in 1:6){
    for(j in 1:6){
      Sigma[l,j]  <-  rho[l,j] *  sig.dem[l] *  sig.dem[j]
    }
  }
  
  for (t in 1:n.occasions){
    for(l in 1:2){
      theta.eps[1:2,l,t] ~ dmnorm(mu.sex[1:2,1], cov = Sigma.sex[1,1:2,1:2]) # annual variation in age ratio
    }
    omega.eps[1:2,t] ~ dmnorm(mu.sex[1:2,2], cov = Sigma.sex[2,1:2,1:2]) # annual variation in sex ratio
  }
  
  for(i in 1:2){ # type (age or sex ratio)
    for(j in 1:2){ # species
      mu.sex[i,j] <- 0
      sig.sex[i,j]  ~ dunif(0,5)
      for(k in 1:2){
        Sigma.sex[i,j,k] <-  rho.sex[j,k,i] *  sig.sex[j,i] *  sig.sex[k,i]
      }
    }
    
    rho.sex[1,1,i] <- 1
    rho.sex[1,2,i] ~ dunif(-1,1)
    rho.sex[2,1,i] <- rho.sex[1,2,i]
    rho.sex[2,2,i] <- 1
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Priors and constraints for hazard/survival models
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for(l in 1:n.sex){
    for(k in 1:n.species){
      for(a in 1:n.age){
        slop[a,k,l] ~ dbeta(10,5)  # Garbage parameter that represents the breeding season survival in 1961 (before mark-recapture data)
        # Harvest hazard prior
        mu.kappa[a,k,l] ~ dnorm(-2, sd = 1)
      }
      # Natural hazard (j = age-season combinations)
      for(j in 1:3){
        mu.sur[k,l,j]  ~ dnorm(-3, sd =  1)
      }
    }
  }
  
  # Regression Coefficients
  for(k in 1:n.species){
    beta.winter[k] ~  dlogis(0, 1)
    for(j in 1:4){ # n.covars
      for(a in 1:2){ #n.age
        beta.surv[k,a,j] ~  dlogis(0, 1) # Winter hazards (natural and harvest)
      }
    }
    for(a in 1:2){
      beta.indirect[k,a] ~  dlogis(0, 1) # Difference between direct and indirect recoveries (Species-Sex)
      beta.surv3[k,a] ~  dlogis(0, 1) # Summer hazards (Species - Covariate Type)
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Hazard rate and Survival models 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for(t in 1:n.occasions){
    hazards.eps[1:6,t] ~ dmnorm( mu.dem[1:6], cov = Sigma[1:6,1:6])
    
    for(l in 1:n.sex){
      for(k in 1:n.species){
        for(a in 1:n.age){
          haz.winter[a,k,l,t] <- exp(l.haz.winter[a,k,l,t] + hazards.eps[k,t])
          haz.kappa[a,k,l,t] <- exp(l.haz.kappa[a,k,l,t] + hazards.eps[k+2,t])
        }
        haz.summer[k,l,t] <- exp(l.haz.summer[k,l,t] + hazards.eps[k+4,t])
      }
    }
  }
  
  for(k in 1:n.species){
    for(l in 1:n.sex){ 
      for (t in 1:(n.occasions-1)){
        l.haz.summer[k,l,t] <- mu.sur[k,l,3] + beta.surv3[k,1] * N.std[t+1,k] +  beta.surv3[k,2] * pond.std[t+1] 
      }
      l.haz.summer[k,l,n.occasions] <- mu.sur[k,l,3]
    }
  }
  
  for (t in 1:n.occasions){
    for(k in 1:n.species){
      for(l in 1:n.sex){
        # Natural mortality hazards
        l.haz.winter[1,k,l,t] <- mu.sur[k,l,1] + beta.surv[k,1,1] * N.std[t,k] + beta.surv[k,1,2] * eps.sal[t,k] 
        l.haz.winter[2,k,l,t] <- mu.sur[k,l,2] + beta.surv[k,2,1] * N.std[t,k] + beta.surv[k,2,2] * eps.sal[t,k]  
        
        # Harvest mortality hazards
        l.haz.kappa[1,k,l,t] <- mu.kappa[1,k,l] + beta.surv[k,1,3] * N.std[t,k] + beta.surv[k,1,4] * pond.std[t] 
        l.haz.kappa[2,k,l,t] <- mu.kappa[2,k,l] + beta.surv[k,2,3] * N.std[t,k] + beta.surv[k,2,4] * pond.std[t]  
        
        haz.kappa[3,k,l,t] <- exp(l.haz.kappa[2,k,l,t] + beta.indirect[k,l] + hazards.eps[k+2,t])
        haz.kappa[4,k,l,t] <- exp(l.haz.kappa[2,k,l,t] + beta.indirect[k,l] + beta.winter[k] + hazards.eps[k+2,t])
        for(a in 1:n.age){
          # Seasonal Risks
          winter_risk[a,k,l,t] <-  haz.kappa[a,k,l,t] + haz.winter[a,k,l,t]
          summer_risk[a,k,l,t] <-  haz.summer[k,l,t]
          annual_risk[a,k,l,t] <-  haz.kappa[a,k,l,t] + haz.winter[a,k,l,t] + haz.summer[k,l,t]
          
          # Seasonal Survival using log-log link function
          #std.winter[a,k,l,t] <-  (s.winter[a,k,l,t] - mean.winter[a,k,l])/ sd.winter[a,k,l]
          
          s.winter[a,k,l,t]  <- exp(-winter_risk[a,k,l,t])
          s.summer[a,k,l,t]  <- exp(-summer_risk[a,k,l,t])
          s.annual[a,k,l,t]  <- exp(-annual_risk[a,k,l,t])
          # Derive the kill rate as the proportion of winter mortality associated with harvest   
          kappa[a,k,l,t] <- ( 1 - s.winter[a,k,l,t]) * (haz.kappa[a,k,l,t]/ winter_risk[a,k,l,t])      # Direct Harvest rate
          natural.winter[a,k,l,t] <- ( 1 - s.winter[a,k,l,t]) * (haz.winter[a,k,l,t]/ winter_risk[a,k,l,t])     # Direct Harvest rate
          relative.summer[a,k,l,t] <-  s.winter[a,k,l,t] * (1 - s.summer[a,k,l,t])
          annual.kappa[a,k,l,t] <- ( 1 - s.annual[a,k,l,t]) * (haz.kappa[a,k,l,t]/ annual_risk[a,k,l,t]) 
          f[a,k,l,t] <- kappa[a,k,l,t] * r[t] * (1 - cr)
        }
        kappa[3,k,l,t] <- ( 1 - s.winter[2,k,l,t]) * (haz.kappa[3,k,l,t]/ winter_risk[2,k,l,t])        # Indirect Harvest rate
        f[3,k,l,t] <- kappa[3,k,l,t] * r[t] * (1 - cr)
        
        kappa[4,k,l,t] <- ( 1 - s.winter[2,k,l,t]) * (haz.kappa[4,k,l,t]/ winter_risk[2,k,l,t])        # Indirect Harvest rate
        f[4,k,l,t] <- kappa[4,k,l,t] * r[t] * (1 - cr)
      }
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Arnold et al. reporting rate prior
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for (t in 1:n.occasions){
    r[t] ~ dbeta(alpha[t], beta[t])
    alpha[t] <- pow(SD[t],-2) * RR[t]
    beta[t] <- pow(SD[t],-2) * (1 - RR[t])
  }
  cr    ~ dbeta(25,100)    # informative prior on crippling rate
  # Age-specific annual survival
  for(k in 1:n.species){
    for(l in 1:n.sex){
      
      for(t in 1:n.occasions){
        SJ[k,l,t] <- s.winter[1,k,l,t] * s.summer[1,k,l,t]
        SA[k,l,t] <- s.winter[2,k,l,t] * s.summer[2,k,l,t]
      }
    }
  }
  
  # Define the multinomial likelihood
  for (t in 1:n.occasions){
    # Band-recovery data: 1955 (breeding season) - 2020 (breeding season)
    for(z in 1:8){
      marr.fall[t,1:(n.occasions+1),z] ~ dmulti(pr.f[age[z],spec[z],sex[z],t,1:(n.occasions+1)], rel.fall[t,z])
    }
    for(z in 1:4){
      marr.wint[t,1:(n.occasions+1),z] ~ dmulti(pr.w[spec1[z],sex1[z],t,1:(n.occasions+1)], rel.wint[t,z]) 
    }
  }
  
  # Define the cell probabilities of the m-array
  # Main diagonal
  for(l in 1:n.species){
    for(k in 1:n.sex){
      for (t in 1:n.occasions){
        pr.f[1,l,k,t,t] <-   f[1,l,k,t] # reported in fall after banding
        pr.f[2,l,k,t,t] <-   f[2,l,k,t] # reported in fall after banding
        # the winter-releases are separated by month of release (i.e., jan, feb, mar)
        pr.w[l,k,t,t] <- 0 
        # pr.w[2,l,k,t,t] <- 0 
      }
      
      for (t in 1:(n.occasions-1)){
        # Above main diagonal
        for (j in (t+1)){
          pr.f[1,l,k,t,j] <- s.winter[1,l,k,t] * s.summer[1,l,k,t] * f[3,l,k,j]  # survived the winter and summer after banding, reported following fall
          pr.f[2,l,k,t,j] <- s.winter[2,l,k,t] * s.summer[2,l,k,t] * f[3,l,k,j]  # survived the winter and summer after banding, reported following fall
          pr.w[l,k,t,j] <- s.summer[1,l,k,t]  * f[4,l,k,j]                      # survived Jan-July after Jan banding, reported in fall (HY)
          #pr.w[2,l,k,t,j] <- s.summer[2,l,k,t]  * f[3,l,k,j]                     # survived Jan-July after Jan banding, reported in fall (AHY)
        }
      }
      for (t in 1:(n.occasions-cut)){
        for (j in (t+2):(t+cut)){
          pr.f[1,l,k,t,j] <- s.winter[1,l,k,t] * s.summer[1,l,k,t] * prod(s.winter[2,l,k,(t+1):(j-1)]) * prod(s.summer[2,l,k,(t+1):(j-1)]) * f[3,l,k,j]
          pr.f[2,l,k,t,j] <- prod(s.winter[2,l,k,t:(j-1)]) * prod(s.summer[2,l,k,t:(j-1)]) *  f[3,l,k,j]
          
          pr.w[l,k,t,j] <- prod(s.winter[2,l,k,(t+1):(j-1)]) * s.summer[1,l,k,t]  * prod(s.summer[2,l,k,(t+1):(j-1)]) * f[4,l,k,j]
          #  pr.w[2,l,k,t,j] <- prod(s.winter[2,l,k,(t+1):(j-1)])                      * prod(s.summer[2,l,k,t:(j-1)])     * f[3,l,k,j]
          
        }
        for (j in (t+cut+1):n.occasions){
          pr.f[1,l,k,t,j] <- 0
          pr.f[2,l,k,t,j] <- 0
          
          pr.w[l,k,t,j] <- 0
          #  pr.w[2,l,k,t,j] <- 0
        }
      }
      for (t in (n.occasions-(cut-1)):n.occasions){
        for (j in (t+2):n.occasions){
          pr.f[1,l,k,t,j] <- s.winter[1,l,k,t] * s.summer[1,l,k,t] * prod(s.winter[2,l,k,(t+1):(j-1)]) * prod(s.summer[2,l,k,(t+1):(j-1)]) * f[3,l,k,j]
          pr.f[2,l,k,t,j] <- prod(s.winter[2,l,k,t:(j-1)]) * prod(s.summer[2,l,k,t:(j-1)]) *  f[3,l,k,j]
          
          pr.w[l,k,t,j] <- prod(s.winter[2,l,k,(t+1):(j-1)]) * s.summer[1,l,k,t]  * prod(s.summer[2,l,k,(t+1):(j-1)]) * f[4,l,k,j]
          #pr.w[2,l,k,t,j] <- prod(s.winter[2,l,k,(t+1):(j-1)])                      * prod(s.summer[2,l,k,t:(j-1)])     * f[3,l,k,j]
        }
      }
      
      # Below main diagonal
      for (t in 1:n.occasions){
        for (j in 1:(t-1)){
          pr.f[1,l,k,t,j] <- 0
          pr.f[2,l,k,t,j] <- 0
          pr.w[l,k,t,j] <- 0
          #  pr.w[2,l,k,t,j] <- 0
        }
      }
      
      # Last column: probability of non-recovery
      for (t in 1:n.occasions){
        pr.f[1,l,k,t,n.occasions+1] <- 1-sum(pr.f[1,l,k,t,1:n.occasions])
        pr.f[2,l,k,t,n.occasions+1] <- 1-sum(pr.f[2,l,k,t,1:n.occasions])
        
        pr.w[l,k,t,n.occasions+1] <- 1-sum(pr.w[l,k,t,1:n.occasions])
        # pr.w[2,l,k,t,n.occasions+1] <- 1-sum(pr.w[2,l,k,t,1:n.occasions])
        
      } #t
    } #k
  } #l
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Duckling Production model
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  for(k in 1:n.species){
    for(j in 1:3){
      beta.np[k,j] ~ dnorm(0, 0.01)
    }
    for(l in 1:n.sex){
      nu.theta[k,l] ~ dlogis(0,1)
    }
    nu.omega[k] ~dlogis(0,1)
    sig.beta[k] ~ T(dt(0, pow(2.5, -2), 1),0,)
    for(i in 1:K){
      beta.spline[k,i] ~ dnorm(0, sd = sig.beta[k])
    }
    
    for(t in 1:n.occasions){
      fem.suscept[t,k]       <-  f[2,k,1,t]/ f[2,k,2,t]                 # relative risk of female to male recovery
      for(l in 1:n.sex){
        juv.suscept[t,l,k]   <-  f[1,k,l,t]/ f[2,k,l,t]                 # relative risk of juvenile to adult recovery (sex specific)
      }
    }
    for(l in 1:n.sex){
      l.theta[k,l,1]  <- nu.theta[k,l] +  beta.np[k,1] * N.std[1,k] + beta.np[k,2] * pond.std[1] + inprod(beta.spline[k,1:K], BM[1,1:K,k])
      logit(theta[k,l,1]) <- l.theta[k,l,1] +  theta.eps[k,l,1]
      
      for(t in 2:n.occasions){
        l.theta[k,l,t]  <- nu.theta[k,l] +  beta.np[k,1] * N.std[t,k] + beta.np[k,2] * pond.std[t] + beta.np[k,3] *  eps.sal[t-1,k] + inprod(beta.spline[k,1:K], BM[t,1:K,k])
        logit(theta[k,l,t]) <- l.theta[k,l,t] +  theta.eps[k,l,t]
      }
    }
    for(t in 1:n.occasions){
      for(l in 1:n.sex){
        logit(tall[k,l,t]) <- nu.theta[k,l]  + inprod(beta.spline[k,1:K], BM[t,1:K,k])
      }
      # Linear sex-ratio model
      l.omega[k,t] <- nu.omega[k]     # proportion of adults female (biased)
      logit(omega[k,t]) <- l.omega[k,t] +  omega.eps[k,t]
      for(l in 1:n.sex){
        # Harvest age-ratio
        age.ratio.prime[k,l,t]  <- theta[k,l,t]/((juv.suscept[t,l,k] * (1-theta[k,l,t])+ theta[k,l,t])) # proportion of ducks juveniles (unbiased)
        age.ratio[k,l,t]  <- age.ratio.prime[k,l,t] / (1 -   age.ratio.prime[k,l,t])              # unbiased age ratio
      }
      # Adult Sex Ratio
      sex.ratio.prime[k,t] <-  omega[k,t]/ (fem.suscept[t,k] * ( 1 - omega[k,t]) + omega[k,t])          # proportion of adults female (unbiased)
      sex.ratio[k,t] <-  sex.ratio.prime[k,t] / (1 -  sex.ratio.prime[k,t] )                      # unbiased sex ratio
    }
    
    for(t in 1:n.occasions){
      for(j in 1:flyway){
        # USA wing-bee data: 1961-2019
        IF[t,j,k] ~ dbin( theta[k,1,t ], F[t,j,k])             # HY Females ~ binom(p = age.rate, N = Females)
        IM[t,j,k] ~ dbin( theta[k,2,t ], M[t,j,k])             # HY Males   ~ binom(p = age.rate, N = Males)
        FA[t,j,k] ~ dbin( omega[k,t ],   T[t,j,k])             #    Females ~ binom(p = sex.rate, N = Adults)
      }
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Pond Abundance Model
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Pond time-series: 1955-2015
  # Strata-specific time-trend
  mu.trend ~ dnorm(0, 0.1)
  sig.trend ~ T(dt(0, pow(2.5, -2), 1),0,)
  
  # Global parameters for random-intercepts (strata)
  b0 ~ dnorm(0, 0.1)
  sig.b0 ~ T(dt(0, pow(2.5, -2), 1),0,)
  tau.b0 <- pow(sig.b0, -2)
  
  mu.sst ~ dnorm(0, 0.001)  
  sig.sst ~ T(dt(0, pow(2.5, -2), 1),0,)
  
  # Priors for covariates
  beta.area ~ dnorm(0, 0.01)
  beta.trend.area ~ dnorm(0, 0.01)
  for (i in 1:n.stratum) {
    nu.pond[i] ~ dnorm(b0,tau.b0)
    beta.sst[i] ~ dnorm(mu.sst,sd = sig.sst)
    phi[i] ~ dunif(-1,1)                      # residual autocorrelation (AR1)
    beta.trend[i] ~ dnorm(mu.trend, sd = sig.trend)  # temporal trend
  }
  sig ~ T(dt(0, pow(2.5, -2), 1),0,)     # spatial variance (partial sill)
  tau ~ T(dt(0, pow(2.5, -2), 1),0,)     # non-spatial variance (nugget)
  decay.pond ~ T(dt(0, pow(2.5, -2), 1),0,)     # spatial decay term
  p.sig ~ T(dt(0, pow(2.5, -2), 1),0,) 
  for (i in 1:n.stratum) {
    sigma[i,i] <- sig + tau
    for (j in (i+1):n.stratum) {
      # model spatial covariance via an exponential semivariogram
      sigma[i,j] <-  sig  * exp(-decay.pond * dist[i,j])
      sigma[j,i] <-  sigma[i,j]
    }
    
    # AR(1) model
    trend[1,i] <- nu.pond[i] +  beta.trend[i] * (1) + beta.sst[i] * pna[1] + beta.area * area[i] + beta.trend.area * area[i] * 1 
    mu.pond[i,1] <- trend[1,i]
    for (t in 2:n.occasions) {
      trend[t,i] <- nu.pond[i] + beta.trend[i] * (t) + beta.sst[i] * pna[t]  + beta.area * area[i]+ beta.trend.area * area[i] * t 
      mu.pond[i,t] <- trend[t,i] + phi[i] * ( pond[t-1,i] - trend[t-1,i] )
    }
  }
  
  for (t in 1:n.occasions) {
    pond[t,1:n.stratum] ~ dmnorm(mu.pond[1:n.stratum,t], cov = sigma[1:n.stratum,1:n.stratum])
    mu.all[t] <- sum(exp(mu.pond[1:n.stratum,t]))
    all.ponds[t] ~ dnorm(mu.all[t],sd = p.sig)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Winter habitat model
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  car_lgm[1:N.lgm] ~ dcar_normal(adj.lgm[1:L.lgm], weights.lgm[1:L.lgm], num.lgm[1:N.lgm], tau.car[1])
  car_chb[1:N.chb] ~ dcar_normal(adj.chb[1:L.chb], weights.chb[1:L.chb], num.chb[1:N.chb], tau.car[2])
  
  beta.flow ~ dnorm(0, 0.01)
  
  beta.depth[1]<- 0
  beta.depth[2]~ dnorm(0, 0.01)
  beta.depth[3]~ dnorm(0, 0.01)
  
  for(i in 1:2){
    MU.SAL[i] <- 0
    tau.car[i] ~ T(dt(0, pow(2.5, -2), 1),0,)
    sig.sal[i] ~ T(dt(0, pow(2.5, -2), 1),0,)
    sig.month[i] ~ T(dt(0, pow(2.5, -2), 1),0,) 
    sig.year[i] ~ T(dt(0, pow(2.5, -2), 1),0,)
    beta.inter[i] ~ dnorm(0, 0.01)
    for(j in 1:4){
      beta.sal[i,j] ~ dnorm(0, 0.01)
    }
    
    for(j in 1:12){
      eps.month[i,j] ~ dnorm(0, sd = sig.month[i])
    }
    rho.sal[i,i] <- 1
  }
  rho.sal[1,2] ~ dunif(-1,1)
  rho.sal[2,1] <- rho.sal[1,2] 
  
  for(i in 1:2){
    for(j in 1:2){
      sigma.sal[i,j] <- rho.sal[i,j] * sig.year[i] * sig.year[j]
    }
  }
  
  for (t in 1:n.occasions) {
    eps.year[t,1:2] ~ dmnorm(MU.SAL[1:2], cov = sigma.sal[1:2,1:2])
    mu.pred[1,t] <-  inprod(beta.sal[1,1:4],pred.covs[t,1:4,1]) + beta.flow * pred.covs[t,5,1] + eps.year[t,1] + beta.inter[1] * pred.covs[t,1,1] * pred.covs[t,4,1]
    mu.pred[2,t] <-  inprod(beta.sal[2,1:4],pred.covs[t,1:4,2]) + eps.year[t,2] + beta.inter[2] * pred.covs[t,1,2] * pred.covs[t,4,2]
    
    eps.sal[t,2] <- mu.pred[1,t] 
    eps.sal[t,1] <- mu.pred[2,t] 
  }
  
  for(i in 1:nobs.lgm){
    mu.lgm[i]    <-  car_lgm[lgm_breaks[i]] + inprod(beta.sal[1,1:4],covs.lgm[i,1:4]) +beta.flow * covs.lgm[i,5] + eps.month[1,month.lgm[i]] + eps.year[year.lgm[i],1] + beta.inter[1] * covs.lgm[i,1] * covs.lgm[i,4] 
    salinity.lgm[i] ~ dnorm(mu.lgm[i], sd = sig.sal[1])
  }
  
  for(i in 1:nobs.chb){
    mu.chb[i] <-   car_chb[chb_breaks[i]] + inprod(beta.sal[2,1:4],covs.chb[i,1:4])+ beta.depth[depths[i]] + eps.month[2,month.chb[i]] + eps.year[year.chb[i],2] + beta.inter[2] * covs.chb[i,1] * covs.chb[i,4] 
    salinity.chb[i] ~ dnorm(mu.chb[i], sd = sig.sal[2])
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Abundance model
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for(t in 1:n.occasions){
    for(k in 1:n.species){
      SF[k,t] <-  sex.ratio.prime[k,t] *  SA[k,1,t]                            # Survivng Females        
      SM[k,t] <-  (1 - sex.ratio.prime[k,t]) *  SA[k,2,t]                      # Surviving Males
      PF[k,t] <-  sex.ratio.prime[k,t] * age.ratio[k,1,t] *  SJ[k,1,t]         # Recruiting Females
      PM[k,t] <-  (1 - sex.ratio.prime[k,t]) * age.ratio[k,2,t] *  SJ[k,2,t]   # Recruiting Males
      lambda[k,t] <-  SF[k,t] +  SM[k,t] + PF[k,t] + PM[k,t]                   # Population growth
      rate[k,t] <- log(lambda[k,t])                                            # log(population growth)
    }
  }
  for(k in 1:n.species){
    N.start[k] ~ dunif(12,15)
    sigma.obs[k] ~ dunif(0,2)
    lN[k,1] <- N.start[k]
    for(t in 2:n.occasions){
      lN[k,t] <- lN[k,t-1] + rate[k,t-1]
    }
    for(t in 1:n.occasions){
      y[k,t] ~  dnorm(lN[k,t], sd = sigma.obs[k])
      N[k,t] <- exp(lN[k,t])
    }
  }
  
})


start_time <- Sys.time()
library(parallel)
library(coda)

nc <- 3 # number of chains

cl_ipm2<-makeCluster(nc,timeout=5184000)

clusterExport(cl_ipm2, c("model_1", "inits", "dat", "constants", "pars"))

for (j in seq_along(cl_ipm2)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(cl_ipm2[j], "inits")
}

out <- clusterEvalQ(cl_ipm2, {
  library(nimble)
  library(coda)
  model <- nimbleModel( code = model_1, constants = constants,  dat =  dat, inits = inits)
  
  model$simulate(c('slop', 'haz.summer', 'haz.winter', 'haz.kappa', 'winter_risk', 'summer_risk', 'annual_risk','s.annual',
                   'l.haz.summer', 'l.haz.winter','l.haz.kappa', 'theta.eps', 'omega.eps', 'hazards.eps',
                   'relative.summer','natural.winter','annual.kappa',
                   's.winter', 's.summer', 'kappa', 'f', 'r', 'SJ', 'SA', 'pr.f', 'pr.w', 'fem.suscept', 'juv.suscept', 'theta', 'omega', 'age.ratio.prime',
                   'age.ratio', 'sex.ratio.prime', 'sex.ratio', 'nu.pond', 'trend', 'mu.pond', 'all.ponds', 
                   'mu.pred','mu.all',
                   'SF', 'SM', 'PF', 'PM', 'lambda', 'rate','lN', 'N'))
  
  model$initializeInfo()
  model$calculate()
  
  modelConf  <- configureMCMC(model, useConjugacy = FALSE, 
                              thin2 = 5, 
                              monitors2 = c('SJ','SA',
                                            'r','cr', 
                                            's.winter','kappa', 
                                            'haz.winter','haz.summer','haz.kappa', 'mu.sur','mu.kappa',
                                            's.summer', 'natural.winter','annual.kappa','relative.summer',
                                            'age.ratio.prime', 'age.ratio', 'sex.ratio.prime', 'sex.ratio', 
                                            'theta', 'omega',  'nu.theta','nu.omega',
                                            'tall','beta.spline',
                                            'beta.surv','beta.surv3',
                                            'lambda','N','SF', 'SM', 'PF', 'PM',
                                            'beta.sst', 'nu.pond', 'phi', 
                                            'beta.trend', 'trend', 'mu.pond', 'all.ponds',
                                            'eps.month','eps.year','beta','car_lgm','car_chb','rho.sal','mu.pred','eps.sal'
                                            
                              ))
  
  modelConf$removeSamplers(c('mu.sur'))
  
  modelConf$addSampler(target = c('mu.sur[1,1,1]','mu.sur[1,1,2]', 'mu.sur[1,1,3]','mu.kappa[1,1,1]','mu.kappa[2,1,1]'), type = 'AF_slice')        
  modelConf$addSampler(target = c('mu.sur[1,2,1]','mu.sur[1,2,2]', 'mu.sur[1,2,3]','mu.kappa[1,1,2]','mu.kappa[2,1,2]'), type = 'AF_slice')       
  modelConf$addSampler(target = c('mu.sur[2,1,1]','mu.sur[2,1,2]', 'mu.sur[2,1,3]','mu.kappa[1,2,1]','mu.kappa[2,2,1]'), type = 'AF_slice')       
  modelConf$addSampler(target = c('mu.sur[2,2,1]','mu.sur[2,2,2]', 'mu.sur[2,2,3]','mu.kappa[1,2,2]','mu.kappa[2,2,2]'), type = 'AF_slice') 
  
  modelMCMC  <- buildMCMC( modelConf)
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  CmodelMCMC$run(500000, thin = 2, thin2 = 2, nburnin = 250000)
  
  return(list( as.mcmc(as.matrix(CmodelMCMC$mvSamples)),
               as.mcmc(as.matrix(CmodelMCMC$mvSamples2))))
  # return(as.mcmc(out1))
  
})
end_time <- Sys.time()
end_time - start_time

samples2 <- list( chain1 =   out[[1]][[1]], chain2 =   out[[2]][[1]], chain3 =   out[[3]][[1]])
samples1 <- list( chain1 =   out[[1]][[2]], chain2 =   out[[2]][[2]], chain3 =   out[[3]][[2]])
