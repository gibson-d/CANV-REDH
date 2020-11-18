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
duck_data <- read.csv('CANV_REDH-stratum_bpop.csv')



pond_data <- subset(pond_data, STRATUM > 18)
pond_data <- subset(pond_data, !grepl('50', STRATUM))

duck_data <- subset(duck_data, STRATUM > 18)
duck_data <- subset(duck_data, !grepl('50', STRATUM))


# Read in dataset
setwd("C:\\Users\\gibsond\\Desktop\\Ducks\\Data\\bpop\\stratum")

bpop <- st_read('WBPHS_stratum_boundaries.shp')

code <- c(11,11,12,10,11,11,7,6,4,5,13,14,15,9,2,3,2,16,8,17,18,19,78,22,24,68,21,77,70,67,25,71,60,76,23,59,27,72,31,51,26,71,32,28,69,58,73,37,33,29,65,38,35,34,64,30,66,52,41,57,36,63,42,39,40,1,46,43,48,44,1,1,47,53,54,56,1,45,55,49,50)
code <- code - 1
bpop$code <- code

keep <- unique(pond_data$STRATUM)

bpop_sub <- bpop[bpop$code %in% keep,]
plot(bpop_sub )

keep2 <- unique(duck_data$STRATUM)
bpop_sub2 <- bpop[bpop$code %in% keep2,]
plot(bpop_sub2 )

W     = st_touches(bpop_sub2$geometry , sparse=FALSE)
listW = mat2listw(W)

adj  = W * 1L
nadj = (rowSums(W))


bpop_coords = bpop_sub2 %>% st_centroid() %>% st_coordinates()
plot(st_geometry(bpop_sub2))
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


bpop_coords.1 <- cbind.data.frame(strata = 1:33, lat = bpop_coords[,2], lon = bpop_coords[,1])

dist <- GeoDistanceInMetresMatrix(bpop_coords.1)/100000

dist.prime <- dist[-c(1:2),-c(1:2)] 

strata <- cbind.data.frame(STRATUM = bpop_sub2$code, STRATA =  1:length(bpop_sub2$code))

duck_data <- left_join(duck_data, strata)


ducks <- cast(duck_data, SPECIES ~ YEAR ~ STRATA, value = 'POP')
vcf   <- cast(duck_data, SPECIES ~ YEAR ~ STRATA, value = 'VCF')
lat <- scale(bpop_coords.1$lat)

duck.p <- cast(duck_data,  YEAR+SPECIES ~ STRATA, value = 'POP')

vcf.means <- apply(vcf, 1, mean, na.rm = TRUE)
vcf.sd    <- apply(vcf, 1, sd, na.rm = TRUE)
vcf_std <-  vcf

for(i in 1:dim(vcf)[2]){
  for(j in 1:dim(vcf)[3]){
    for(k in 1:dim(vcf)[1]){
      vcf_std[k,i,j] <- (vcf[k,i,j] - vcf.means[k]) / vcf.sd[k]
    }
  }
}
vcf_std[is.na(vcf_std)] <- 0


# Average (April-July) Pacific SST (ENO) 1955 - 2019
sst <- c(-0.7475,-0.535,1,0.72,0.07,0.0525,0.1825,-0.195,0.49,-0.5225,0.6725,0.3725,-0.16,0.095,0.5425,-0.175,
         -0.78,0.7825,-0.655,-0.81,-0.9175,-0.155,0.2825,-0.285,0.15,0.3925,-0.305,0.66,0.845,-0.4225,-0.67,-0.035,
         1.1625,-0.9475,-0.53,0.305,0.52,0.8625,0.565,0.395,0.0425,-0.3175,0.9625,0.1325,-1.035,-0.6775,-0.1975,0.5175,
         -0.095,0.2725,0.1925,-0.0425,-0.3675,-0.64,0.1675,-0.3175,-0.4875,-0.06,-0.3,0.1825,1.15,0.2925,0.305,-0.0925,0.55)



library(nimble)
code <- nimbleCode({
  
  # priors
  for(k in 1:n.species){
        mu.fec[k] ~ dnorm(0, 0.01)
      beta.sst[k] ~ dnorm(0, 0.01)
      beta.lat[k] ~ dnorm(0, 0.01)
      beta.lat.sst[k] ~ dnorm(0, 0.01)
      omega[k] ~ dbeta(1,1)
      kappa[k] ~ dbeta(1,1)
     mean.p[k] ~ dbeta(1,1)
       mu.p[k] <- logit(mean.p[k])
   beta.vis[k] ~ dnorm(0, 0.01)
  
  for (i in 1:n.stratum) {
         mu.lambda0[k,i] ~ dunif(0,20)
            lambda0[k,i] <- exp(mu.lambda0[k,i])
                N[k,1,i] ~ dpois(lambda0[k,i])
    for (t in 2:n.year) {
            Ilambda[k,t-1,i] <- sum(E[k,t-1,1:n.stratum]/nadj[1:n.stratum] * adj[i,1:n.stratum])
          fecundity[k,t-1,i] <-exp( mu.fec[k] + beta.sst[k] * sst[t-1] + beta.lat[k] * lat[i] + beta.lat.sst[k] * lat[i] * sst[t-1]) 
                  S[k,t-1,i] ~ dbin(omega[k], N[k,t-1,i])
                  E[k,t-1,i] ~ dbin(kappa[k], S[k,t-1,i])
                  I[k,t-1,i] ~ dpois(Ilambda[k,t-1,i])
                  R[k,t-1,i] ~ dpois(0.5 * fecundity[k,t-1,i] * N[k,t-1,i])
            
                  N[k,t,i] <-  S[k,t-1,i] - E[k,t-1,i] + R[k,t-1,i] + I[k,t-1,i]
    }
    
    for (t in 1:n.year) {
        logit(pobs[k,t,i]) <- mu.p[k] + beta.vis[k] * vcf[k,t,i]
                 y[k,t,i] ~ dbin(pobs[k,t,i], N[k,t,i] )
      }
   }
    
    for (t in 1:n.year) {
      tN[k,t] <- sum(N[k,t,1:n.stratum])
    }
  }
  
})


N.init <- round(ducks + 10000)

S.init <- N.init[,-61,] - 2000
R.init <- pmax(N.init[,-1,] - S.init,0)
I.init <- S.init - 7000
E.init <- I.init



initsFunction <- function()list( 
                      mu.lambda0 = log(ducks[,1,] + 1),
                      beta.sst = c(0,0), beta.lat = c(0,0), beta.lat.sst = c(0,0), beta.vis = c(0,0),
                      mu.fec = c(0,0), omega = c(.7,.7), kappa = c(.2,.2), mean.p = c(.9,.9),
                      N = N.init, S = S.init, I = I.init, E = E.init, R = R.init
                    )

dat <- list(y =  round(ducks,0), 
            dist = dist, 
            sst = sst[1:61],
            vcf = vcf_std,
            lat = lat[,1],
            adj = adj, nadj = nadj
            )
constants <- list( n.year = 61, n.stratum = 33, n.species = 2)

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
  
  
  modelConf  <- configureMCMC(model,
                              useConjugacy = FALSE,
                              monitors2 = c( 'N','omega','kappa','tN', 'E','I','fecundity',
                                             'beta.sst','beta.lat','beta.lat.sst', 'beta.vis'
                                             
                                            ))
  
  
  modelMCMC  <- buildMCMC( modelConf  )
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  out1 <- runMCMC(CmodelMCMC,nburnin = 2500, niter = 5000, thin = 2)
  
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

N.est <- MCMCsummary(mcmcList2,'N')

N.est$YEAR    <- rep(seq(from = 1955, to = 2015), each = 2, times = 33) 
N.est$STRATA  <- rep(c(1:33), each = 2 * 61)
N.est$SPECIES <- rep(c('Canvasback','Redhead'), times = 33 * 61)

N.est_me <- left_join(N.est, duck_data)

ggplot(data = N.est_me, aes(x = YEAR, y =  mean)) +
  geom_pointrange(aes(x = YEAR, y = mean, ymin = `2.5%`, ymax = `97.5%`, color = SPECIES)) +
  facet_wrap(~ STRATUM, ncol = 6, scale = 'free_y') +
  scale_color_manual(values = c('dodgerblue','goldenrod')) +
  geom_point(aes(x = YEAR, y = POP, fill = SPECIES), shape = 21) +
  scale_fill_manual(values = c('black','white')) +
  theme_bw()

tN.est        <- MCMCsummary(mcmcList2,'tN')
tN.est$Year   <- rep(c(1955:2015), each = 2)
tN.est$Count  <- rowSums(duck.p[,-c(1:2)])
tN.est$Species <- duck.p$SPECIES

inits <-  t(apply(N.init, c(1,2), sum))
inits <- c(inits[,1], inits[,2])
inits.df<- cbind.data.frame(Year = rep(c(1955:2015), times = 2), Species = rep(unique(tN.est$Species ), each= 61), Count = inits)


ggplot(data =  tN.est, aes(x = Year, y =  mean)) +
  geom_pointrange(aes(x = Year, y = mean, ymin = `2.5%`, ymax = `97.5%`), color = 'dodgerblue') +
  geom_point(aes(x = Year, y = Count ), color = 'goldenrod') + facet_wrap(~Species) +
  geom_point(data = inits.df, aes(X = Year, y = Count), color = 'firebrick') +
  theme_bw()


