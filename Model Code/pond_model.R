#------------------------------------------------------------------------------
# Data manipulation of banding data 
#------------------------------------------------------------------------------
rm(list=ls())

library(ggplot2)
library(dplyr)
library(reshape)
library(Imap)
library(sf)
library(spdep)

setwd("C:\\Users\\gibsond\\Desktop\\Ducks\\Data\\bpop")
pond_data <- read.csv('Ponds-stratum_bpop.csv')


pond_data <- subset(pond_data, STRATUM > 18)
pond_data <- subset(pond_data, !grepl('50', STRATUM))

# Read in shape file
setwd("C:\\Users\\gibsond\\Desktop\\Ducks\\Data\\bpop\\stratum")

bpop <- st_read('WBPHS_stratum_boundaries.shp')
code <- c(11,11,12,10,11,11,7,6,4,5,13,14,15,9,2,3,2,16,8,17,18,19,78,22,24,68,21,77,70,67,25,71,60,76,23,59,27,72,31,51,26,71,32,28,69,58,73,37,33,29,65,38,35,34,64,30,66,52,41,57,36,63,42,39,40,1,46,43,48,44,1,1,47,53,54,56,1,45,55,49,50)
code <- code - 1
bpop$code <- code

keep <- unique(pond_data$STRATUM)
bpop_sub <- bpop[bpop$code %in% keep,]
plot(bpop_sub )

W     = st_touches(bpop_sub$geometry , sparse=FALSE)
listW = mat2listw(W)

bpop_coords = bpop_sub %>% st_centroid() %>% st_coordinates()
plot(st_geometry(bpop_sub))
plot(listW, bpop_coords, add=TRUE, col='blue', pch=16)

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

bpop_coords.1 <- cbind.data.frame(strata = 1:2, lat = bpop_coords[,2], lon = bpop_coords[,1])

# Create Distance Matrix
dist <- GeoDistanceInMetresMatrix(bpop_coords.1)/100000
dist.prime <- dist[-c(1:2),-c(1:2)] 


strata <- cbind.data.frame(STRATUM = bpop_sub$code, STRATA =  1:length(bpop_sub$code))
pond_data <- left_join(pond_data, strata)

ponds <- cast(pond_data,  YEAR ~ STRATA, value = 'POP')
vcf   <- cast(pond_data,  YEAR ~ STRATA, value = 'VCF')
lat <- scale(bpop_coords.1$lat)
ponds $YEAR <- NULL
ponds  <- (data.matrix(ponds/1000))

vcf$YEAR <- NULL
vcf <- data.matrix(vcf)
vcf <- (vcf - mean(vcf, na.rm = TRUE))/ sd(vcf, na.rm = TRUE)
vcf[is.na(vcf)] <- 0

# Average (April-July) Pacific SST (ENO) 1955 - 2019
sst <- c(-0.7475,-0.535,1,0.72,0.07,0.0525,0.1825,-0.195,0.49,-0.5225,0.6725,0.3725,-0.16,0.095,0.5425,-0.175,
         -0.78,0.7825,-0.655,-0.81,-0.9175,-0.155,0.2825,-0.285,0.15,0.3925,-0.305,0.66,0.845,-0.4225,-0.67,-0.035,
         1.1625,-0.9475,-0.53,0.305,0.52,0.8625,0.565,0.395,0.0425,-0.3175,0.9625,0.1325,-1.035,-0.6775,-0.1975,0.5175,
         -0.095,0.2725,0.1925,-0.0425,-0.3675,-0.64,0.1675,-0.3175,-0.4875,-0.06,-0.3,0.1825,1.15,0.2925,0.305,-0.0925,0.55)


library(nimble)
code <- nimbleCode({

     # priors
     mu.trend  ~ dnorm(0, 0.0000001)
     sig.trend ~ T(dt(0, pow(2.5, -2), 1),0,)
     tau.trend <- pow(sig.trend, -2)

     b0  ~ dnorm(0, 0.0000001)
     sig.b0 ~ T(dt(0, pow(2.5, -2), 1),0,)
     tau.b0 <- pow(sig.b0, -2)
      
     beta.sst ~ dnorm(0, 0.000001)
     beta.lat ~ dnorm(0, 0.000001)
     beta.lat.sst ~ dnorm(0, 0.000001)
      for (i in 1:n.stratum) {
   
       nu.pond[i] ~ dnorm(b0,tau.b0)
           phi[i] ~ dunif(-1.1,1.1)             # residual autocorrelation
    beta.trend[i] ~ dnorm(mu.trend, tau.trend)  # temporal trend
      }
      
         sig ~ T(dt(0, pow(2.5, -2), 1),0,) # partial sill
         tau ~ T(dt(0, pow(2.5, -2), 1),0,) # nugget effect
       decay ~ T(dt(0, pow(2.5, -2), 1),0,)
       
          for (i in 1:n.stratum) {
              sigma[i,i] <- sig + tau
            for (j in (i+1):n.stratum) {
            # model spatial covariance via an exponential semivariogram 
              sigma[i,j] <-  sig  * exp(-decay * dist[i,j])
              sigma[j,i] <-  sigma[i,j]
            }
          
         trend[1,i] <- nu.pond[i] + beta.trend[i] * log(1) + 
                                    #beta.vcf[i] * vcf[1,i] + 
                                    beta.sst * sst[1] +
                                    beta.lat * lat[i] + 
                                    beta.lat.sst * sst[1] * lat[i] 
         
         mu.pond[i,1] <- trend[1,i]
               for (t in 2:n.year) {
                   trend[t,i] <- nu.pond[i] + beta.trend[i] * log(t) +
                                              #beta.vcf[i] * vcf[t,i] + 
                                              beta.sst * sst[t] + 
                                              beta.lat * lat[i] + 
                                              beta.lat.sst * sst[t] * lat[i] 
                 mu.pond[i,t] <- trend[t,i] + phi[i] * ( pond[t-1,i] - trend[t-1,i] )
      }
    }
    for (t in 1:n.year) {
          pond[t,1:n.stratum] ~ dmnorm(mu.pond[1:n.stratum,t], cov = sigma[1:n.stratum,1:n.stratum])
                 all.ponds[t] <- sum(mu.pond[1:n.stratum,t]) 
    }

})


initsFunction <- function()list( 
  mu.trend = 0, sig.trend = 10,  
  mu.sst = 0, sig.sst = 10,
  b0 = 200, sig.b0 = 100,
  mu.vcf = 0, sig.vcf = 10,
  decay = .5, sig = 10, tau = 10
  )

dat <- list(pond = ponds[,-c(1:2)], vcf = vcf[,-c(1:2)],dist = dist.prime, sst = sst[1:61],lat = lat[-c(1:2)])
constants <- list( n.year = 61, n.stratum = 24)

pars <- c('phi','decay','sig','tau')

inits <- initsFunction() 

library(parallel)
library(coda)
set.seed(2)
nc <- 4 # number of chains

clduck<-makeCluster(nc,timeout=5184000)

clusterExport(clduck, c("code", "inits", "dat", "constants", "pars"))

for (j in seq_along(clduck)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(clduck[j], "inits")
}

out <- clusterEvalQ(clduck, {
  library(nimble)
  library(coda)
  model <- nimbleModel( code = code, constants = constants,  dat =  dat, inits = inits)
  model$simulate()
  modelConf  <- configureMCMC(model,
                              useConjugacy = FALSE,
                              monitors2 = c('phi','nu.pond','mu.pond','sigma','all.ponds',
                                            'beta.trend', 'beta.sst','beta.lat','beta.lat.sst',
                                            'trend','mu.trend','sig.trend',
                                             'sig.b0',
                                            'decay','sig','tau',
                                            'b0'))
 
  
  modelMCMC  <- buildMCMC( modelConf  )
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  out1 <- runMCMC(CmodelMCMC,nburnin = 5000, niter = 10000, thin = 2)
  
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

  pond.est <- MCMCsummary(mcmcList2,'mu.pond')

  pond_data_sub <- subset(pond_data, STRATUM < 74)
  pond_data_sub$STRATA <- pond_data_sub$STRATA -2
  pond.est$YEAR    <- rep(seq(from = 1955, to = 2015), each = 24) 
  pond.est$STRATA <- rep(c(1:24), times = 61)
  
  pond.est_me <- left_join(pond.est, pond_data_sub)

  ggplot(data = pond.est_me, aes(x = YEAR, y =  mean)) +
  geom_pointrange(aes(x = YEAR, y = mean, ymin = `2.5%`, ymax = `97.5%`), color = 'dodgerblue') +
  facet_wrap(~ STRATUM, ncol = 6, scale = 'free_y') +
  geom_point(aes(x = YEAR, y = POP/1000), color = 'goldenrod') +
  theme_bw()

  all.pond.est <- MCMCsummary(mcmcList2,'all.ponds')
  all.pond.est $YEAR <- c(1955:2015)
  all.pond.est $Count <- rowSums(ponds[,-(1:2)])
  
  
    ggplot(data =  all.pond.est, aes(x = YEAR, y =  mean)) +
    geom_pointrange(aes(x = YEAR, y = mean, ymin = `2.5%`, ymax = `97.5%`), color = 'dodgerblue') +
 
    geom_point(aes(x = YEAR, y = Count ), color = 'goldenrod') +
    theme_bw()
  
  
  
