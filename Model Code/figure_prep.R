#load('model_6_18_2021_not_fixed.rdata')

mcmcList1 <- file_list[[1]]
mcmcList2 <- file_list[[2]]
###########################################################

library(MCMCvis)

theta_inter <- MCMCsummary(file_list[[2]],'nu.theta')
p_beta      <- MCMCsummary(file_list[[2]],'beta.np', hpd_prob = .95, HPD = TRUE)
s_beta      <- MCMCsummary(file_list[[1]],'beta.spline', hpd_prob = .95, HPD = TRUE)

haz_inter <- MCMCsummary(file_list[[2]],'mu.sur')
kap_inter <- MCMCsummary(file_list[[2]],'mu.kappa')
l_beta    <- MCMCsummary(file_list[[2]],'beta.surv', hpd_prob = .95, HPD = TRUE)
q_beta    <- MCMCsummary(file_list[[2]],'beta.surv3', hpd_prob = .95, HPD = TRUE)

sal       <- MCMCsummary(file_list[[1]],'mu.pred')
sal2       <- MCMCsummary(file_list[[1]],'eps.sal')
all.ponds <- MCMCsummary(file_list[[1]],'all.ponds')


mu.theta <- matrix(theta_inter$mean, 2,2)
sd.theta <- matrix(theta_inter$sd, 2,2)

beta.theta<- matrix(p_beta$mean, 2,3)
sd.b.theta<- matrix(p_beta$sd, 2,3)

mu.sur <- array(haz_inter$mean, dim = c(2,2,3))
sd.sur <- array(haz_inter$sd, dim = c(2,2,3))

mu.kappa <- array(kap_inter$mean, dim = c(2,2,2))
sd.kappa <- array(kap_inter$sd, dim = c(2,2,2))

beta.surv <- array(l_beta$mean, dim = c(2,2,4))
beta.sd   <- array(l_beta$sd, dim = c(2,2,4))
beta.surv2 <- matrix(q_beta$mean, 2,2) 
beta.sd2   <- matrix(q_beta$sd, 2,2)
sal_spl <- matrix(sal2$mean, ncol = 2, byrow = TRUE)

red_n  <- seq(from = min(scale(REDS)), to = max(scale(REDS)), length.out = 100)
can_n  <- seq(from = min(scale(CANS)), to = max(scale(CANS)), length.out = 100)
pond_n <- seq(from = min(pond.std), to = max(pond.std), length.out = 100)
sal_chb <- seq(from = min(sal_spl[,1]), to = max(sal_spl[,1]), length.out = 100)
sal_lgm <- seq(from = min(sal_spl[,2]), to = max(sal_spl[,2]), length.out = 100)

n <- 50000
plot.sal <- plot.nw  <- plot.nk <- plot.pondk <- array(NA, dim = c(2,2,2,100,n))
plot.ns <- plot.pond <- array(NA, dim = c(2,2,100,n))
plot.pn <- plot.pp <- plot.pds  <- array(NA, dim = c(2,2,100,n))

sim.sal  <- rbind(sal_chb,sal_lgm)
sim.n    <- rbind(can_n,red_n)
sim.pond <- pond_n 

mu_sa<-mu_nw<-mu_nk<-mu_pk <- array(NA, dim = c(2,2,2,100,3))
mu_ns<-mu_ps <- array(NA, dim = c(2,2,100,3))

mu_pp <- mu_pn <-mu_pds <- array(NA, dim = c(2,2,100,3))

for(k in 1:2){
  for(l in 1:2){
    for(i in 1:100){
      plot.sal[k,l,1,i,]  <-  exp(rnorm(n, mu.sur[k,l,1],sd.sur[k,l,1]) + rnorm(n, beta.surv[k,1,2], beta.sd[k,1,2]) * sim.sal[k,i] )
      plot.sal[k,l,2,i,]  <-  exp(rnorm(n, mu.sur[k,l,2],sd.sur[k,l,2]) + rnorm(n, beta.surv[k,2,2], beta.sd[k,2,2]) * sim.sal[k,i] )
      
      plot.nw[k,l,1,i,]   <-  exp(rnorm(n, mu.sur[k,l,1],sd.sur[k,l,1]) + rnorm(n, beta.surv[k,1,1], beta.sd[k,1,1]) * sim.n[k,i])
      plot.nw[k,l,2,i,]   <-  exp(rnorm(n, mu.sur[k,l,2],sd.sur[k,l,2]) + rnorm(n, beta.surv[k,2,1], beta.sd[k,2,1]) * sim.n[k,i])
      
      plot.ns[k,l,i,]     <-  exp(rnorm(n, mu.sur[k,l,3],sd.sur[k,l,3]) + rnorm(n, beta.surv2[k,1], beta.sd2[k,1]) * sim.n[k,i])
      plot.pond[k,l,i,]   <-  exp(rnorm(n, mu.sur[k,l,3],sd.sur[k,l,3]) + rnorm(n, beta.surv2[k,2], beta.sd2[k,2]) * sim.pond[i])
      
      plot.nk[k,l,1,i,]   <-  exp(rnorm(n, mu.kappa[1,k,l], sd.kappa[1,k,l]) + rnorm(n,beta.surv[k,1,3],beta.sd[k,1,3]) * sim.n[k,i])
      plot.nk[k,l,2,i,]   <-  exp(rnorm(n, mu.kappa[2,k,l], sd.kappa[2,k,l]) + rnorm(n,beta.surv[k,2,3],beta.sd[k,2,3]) * sim.n[k,i])
      
      plot.pondk[k,l,1,i,] <- exp( rnorm(n, mu.kappa[1,k,l], sd.kappa[1,k,l]) + rnorm(n, beta.surv[k,1,4],beta.sd[k,1,4]) * sim.pond[i])
      plot.pondk[k,l,2,i,] <- exp( rnorm(n, mu.kappa[2,k,l], sd.kappa[2,k,l]) + rnorm(n, beta.surv[k,2,4],beta.sd[k,2,4]) * sim.pond[i])
      
      plot.pn[k,l,i,]   <-  plogis(rnorm(n, mu.theta[k,l],sd.theta[k,l]) + rnorm(n, beta.theta[k,1], sd.b.theta[k,1]) * sim.n[k,i])
      plot.pp[k,l,i,]   <-  plogis(rnorm(n, mu.theta[k,l],sd.theta[k,l]) + rnorm(n, beta.theta[k,2], sd.b.theta[k,2]) * sim.pond[i])
      plot.pds[k,l,i,]   <-  plogis(rnorm(n, mu.theta[k,l],sd.theta[k,l]) + rnorm(n, beta.theta[k,3], sd.b.theta[k,3]) * sim.sal[k,i])
      
      for(a in 1:2){
        mu_sa[k,l,a,i,1:3] <- quantile(     plot.sal[k,l,a,i,], c(0.025,.5,.975)) 
        mu_nw[k,l,a,i,1:3] <- quantile(      plot.nw[k,l,a,i,], c(0.025,.5,.975)) 
        mu_nk[k,l,a,i,1:3] <- quantile(      plot.nk[k,l,a,i,], c(0.025,.5,.975)) 
        mu_pk[k,l,a,i,1:3] <- quantile(   plot.pondk[k,l,a,i,], c(0.025,.5,.975)) 
      }
      mu_ps[k,l,i,1:3] <- quantile(    plot.pond[k,l,i,], c(0.025,.5,.975)) 
      mu_ns[k,l,i,1:3] <- quantile(      plot.ns[k,l,i,], c(0.025,.5,.975)) 
      
      mu_pp[k,l,i,1:3] <- quantile(    plot.pp[k,l,i,], c(0.025,.5,.975)) 
      mu_pn[k,l,i,1:3] <- quantile(      plot.pn[k,l,i,], c(0.025,.5,.975)) 
      mu_pds[k,l,i,1:3] <- quantile(      plot.pds[k,l,i,], c(0.025,.5,.975)) 
    }
  }
}

mu_sa.melt <- cbind.data.frame(lower = melt(mu_sa)[1:800,6], mean = melt(mu_sa)[801:1600,6], upper = melt(mu_sa)[1601:2400,6])
mu_nw.melt <- cbind.data.frame(lower = melt(mu_nw)[1:800,6], mean = melt(mu_nw)[801:1600,6], upper = melt(mu_nw)[1601:2400,6])
mu_nk.melt <- cbind.data.frame(lower = melt(mu_nk)[1:800,6], mean = melt(mu_nk)[801:1600,6], upper = melt(mu_nk)[1601:2400,6])
mu_pk.melt <- cbind.data.frame(lower = melt(mu_pk)[1:800,6], mean = melt(mu_pk)[801:1600,6], upper = melt(mu_pk)[1601:2400,6])

mu_ns.melt <- cbind.data.frame(lower = melt(mu_ns)[1:400,5], mean = melt(mu_ns)[401:800,5], upper = melt(mu_ns)[801:1200,5])
mu_ps.melt <- cbind.data.frame(lower = melt(mu_ps)[1:400,5], mean = melt(mu_ps)[401:800,5], upper = melt(mu_ps)[801:1200,5])

mu_pn.melt <- cbind.data.frame(lower = melt(mu_pn)[1:400,5], mean = melt(mu_pn)[401:800,5], upper = melt(mu_pn)[801:1200,5])
mu_pp.melt <- cbind.data.frame(lower = melt(mu_pp)[1:400,5], mean = melt(mu_pp)[401:800,5], upper = melt(mu_pp)[801:1200,5])
mu_pds.melt <- cbind.data.frame(lower = melt(mu_pds)[1:400,5], mean = melt(mu_pds)[401:800,5], upper = melt(mu_pds)[801:1200,5])

species = c('CANVASBACK','REDHEAD')
sex     = c('F','M')
age     = c('HY','AHY')
age2     = c('HY:D','AHY:D','AHY:I')

sal_df   <- cbind.data.frame(expand.grid(species = species, sex = sex, age = age, range = 1:100),   mean = mu_sa.melt[,2], lower = mu_sa.melt[,1], upper = mu_sa.melt[,3])
nw_df    <- cbind.data.frame(expand.grid(species = species, sex = sex, age = age, range = 1:100),   mean = mu_nw.melt[,2], lower = mu_nw.melt[,1], upper = mu_nw.melt[,3])
nk_df    <- cbind.data.frame(expand.grid(species = species, sex = sex, age = age, range = 1:100),   mean = mu_nk.melt[,2], lower = mu_nk.melt[,1], upper = mu_nk.melt[,3])
pondk_df <- cbind.data.frame(expand.grid(species = species, sex = sex, age = age, Ponds = (sim.pond * 1585.194) + 5076.196),   mean = mu_pk.melt[,2], lower = mu_pk.melt[,1], upper = mu_pk.melt[,3])

ns_df    <- cbind.data.frame(expand.grid(species = species, sex = sex, range = 1:100),   mean = mu_ns.melt[,2], lower = mu_ns.melt[,1], upper = mu_ns.melt[,3])
pond_df  <- cbind.data.frame(expand.grid(species = species, sex = sex,Ponds = (sim.pond * 1585.194) + 5076.196),   mean = mu_ps.melt[,2], lower = mu_ps.melt[,1], upper = mu_ps.melt[,3])


prod.n_df  <- cbind.data.frame(expand.grid(species = species, sex = sex, range = 1:100),   mean = mu_pn.melt[,2], lower = mu_pn.melt[,1], upper = mu_pn.melt[,3])
prod.p_df  <- cbind.data.frame(expand.grid(species = species, sex = sex,Ponds = (sim.pond * 1585.194) + 5076.196),   mean = mu_pp.melt[,2], lower = mu_pp.melt[,1], upper = mu_pp.melt[,3])
prod.ds_df  <- cbind.data.frame(expand.grid(species = species, sex = sex, range = 1:100),   mean = mu_pds.melt[,2], lower = mu_pds.melt[,1], upper = mu_pds.melt[,3])

library(ggplot2)
library(see)

mean_sal <- c((mean(MCMCsummary(mcmcList1, 'car_chb')$mean)),(mean(MCMCsummary(mcmcList1, 'car_lgm')$mean))) #c(alpha$mean[2],alpha$mean[1]) #<- c(9.483431,39.22134)
std_sal  <- c(1,1) #<- c(5.051736,12.64485)

sal_normal <- rbind.data.frame(cbind.data.frame(species = 'CANVASBACK', range = 1:100, Salinity = exp(sal_chb * std_sal[1] + mean_sal[1])),
                               cbind.data.frame(species = 'REDHEAD', range = 1:100, Salinity = exp(sal_lgm * std_sal[2] + mean_sal[2])))

N_normal <- rbind.data.frame(cbind.data.frame(species = 'CANVASBACK', range = 1:100,  N = can_n *  117381.7 +  467738),
                             cbind.data.frame(species = 'REDHEAD', range = 1:100,N = red_n * 358077.4 + 790629.4))

nw_df     <- left_join(nw_df , N_normal) 
ns_df     <- left_join(ns_df , N_normal) 
nk_df     <- left_join(nk_df , N_normal) 
prod.n_df <- left_join(prod.n_df, N_normal)

prod.ds_df <- left_join(prod.ds_df, sal_normal)
sal_df <- left_join(sal_df, sal_normal)

len <- 225000
median_hdi = function(x, ...) {
  HPDinterval(mcmc(x), prob = .9, ...) %>% 
    data.frame() %>% 
    select(ymin = lower, ymax = upper) %>% 
    cbind(y = median(x, ...))
}

#simple posterior violin plot
ggposterior = function(.data, .aes) {
  ggplot(
    .data,
    .aes
  ) + 
    geom_violin(linetype=0, fill="skyblue") + 
    stat_summary(fun.data=median_hdi) +
    coord_flip() 
}

covs.mort    <- MCMCpstr(mcmcList2, params = c('beta.surv'), type = 'chains')
HPD_1        <- MCMCsummary(mcmcList2, params = c('beta.surv'), hpd_prob = .9, HPD= TRUE)
covs.mort2   <- MCMCpstr(mcmcList2, params = c('beta.surv3'), type = 'chains')
HPD_2        <- MCMCsummary(mcmcList2, params = c('beta.surv3'), hpd_prob = .9, HPD= TRUE)
covs.prod    <- MCMCpstr(mcmcList2, params = c('beta.np'), type = 'chains')
HPD_3        <- MCMCsummary(mcmcList2, params = c('beta.np'), hpd_prob = .9, HPD= TRUE)
covs.mort.m  <- melt(covs.mort)
covs.mort.m2 <- melt(covs.mort2)
covs.prod.p  <- melt(covs.prod)

covm1           <- cbind.data.frame(covs.mort.m,  expand.grid(species = c('CANVASBACK','REDHEAD'), age =c('HY','AHY'), Covar = c('W:Abundance','W:Salinity','H:Abundance','H:Ponds'), iter = 1:len))
covm1$Vital     <- rep(c('Natural Winter','Harvest'), each = 8)
covm1$Covariate <- rep(c('Abundance','Salinity','Abundance','Ponds'), each = 4)

covm2 <- cbind.data.frame(covs.mort.m2, expand.grid(species = c('CANVASBACK','REDHEAD'), age =c('Both'), Covar = c('Abundance','Ponds'), iter = 1:len))
covm2$Vital <- 'Natural Summer'
covm3 <- cbind.data.frame(covs.prod.p, expand.grid(species = c('CANVASBACK','REDHEAD'), Covar = c('Abundance','Ponds','Salinity'), iter = 1:len))
covm3$Vital <- 'Production'

labs_m1 <- covm1 %>%
  group_by(species, age,Vital,Covariate) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs_m1$perc <- ifelse(labs_m1$perc < .5, 1 - labs_m1$perc,labs_m1$perc)

covm1_plot <- ggplot() + 
  geom_violin(data = subset(covm1, Vital =='Harvest'),aes(x = Covariate, y = value, fill = age),linetype=0) +  
  stat_summary(data = subset(covm1, Vital =='Harvest'),aes(x = Covariate, y = value, group = age),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  #geom_text(data = subset(labs_m1, Vital =='Harvest'), aes(x = Covariate,  y = max, label=round(perc,2),group = age), position = position_dodge2(1)) +
  #geom_pointrange2(aes(x = age, y = mean, ymin = `90%_HPDL`, ymax = `90%_HPDU`, fill = age), shape = 21, size = 0.75) + 
  facet_wrap( Vital~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
covm1_plot

covm1b_plot <- ggplot() + 
  geom_violin(data = subset(covm1, Vital !='Harvest'),aes(x = Covariate, y = value, fill = age),linetype=0) +
  stat_summary(data = subset(covm1, Vital !='Harvest'),aes(x = Covariate, y = value, group = age), color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  #geom_text(data = subset(labs_m1, Vital !='Harvest'), aes(x = Covariate,  y = max, label=round(perc,2),group = age), position = position_dodge2(1)) +
  #geom_pointrange2(aes(x = age, y = mean, ymin = `90%_HPDL`, ymax = `90%_HPDU`, fill = age), shape = 21, size = 0.75) + 
  facet_wrap( Vital~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  coord_cartesian(ylim = c(-2.5, 3))+
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
covm1b_plot

labs_m2 <- covm2 %>%
  group_by(species, age,Covar) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs_m2$perc <- ifelse(labs_m2$perc < .5, 1 - labs_m2$perc,labs_m2$perc)

covm2_plot <- ggplot(data = covm2) + 
  geom_violin(aes(x = Covar, y = value), fill= '#8e44ad', linetype = 0) +
  stat_summary(aes(x = Covar, y = value), color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  #geom_text(data = labs_m2, aes(x = Covar,  y = max, label=round(perc,2)), position = position_dodge2(1)) +
  #geom_pointrange2(aes(x = age, y = mean, ymin = `90%_HPDL`, ymax = `90%_HPDU`, fill = age), shape = 21, size = 0.75) + 
  facet_wrap(Vital ~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
covm2_plot


labs_m3 <- covm3 %>%
  group_by(species,Covar) %>%
  summarize(total = n(),
            max = max(value) + .05,
            above = sum(value > 0),
            perc = above/total)

labs_m3$perc <- ifelse(labs_m3$perc < .5, 1 - labs_m3$perc,labs_m3$perc)

covm3_plot <- ggplot(data = covm3) + 
  geom_violin(aes(x = Covar, y = value), fill= '#8e44ad',linetype = 0) +
  stat_summary(aes(x = Covar, y = value), color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  #geom_text(data = labs_m3, aes(x = Covar,  y = max, label=round(perc,2)), position = position_dodge2(1)) +
  #geom_pointrange2(aes(x = age, y = mean, ymin = `90%_HPDL`, ymax = `90%_HPDU`, fill = age), shape = 21, size = 0.75) + 
  facet_wrap(Vital ~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
covm3_plot

library(ggpubr)
plots <- ggarrange(covm1_plot,covm1b_plot, covm2_plot,covm3_plot, ncol = 2, nrow = 2, labels = c('A','B','C','D'), common.legend = TRUE)
#ggsave(plots, file = 'beta_coefs.tiff', width = 14, height = 8, dpi = 320)
###########################################################

winter_nat  <- MCMCsummary(file_list[[1]], 'haz.winter')
summer_nat  <- MCMCsummary(file_list[[1]], 'haz.summer')
winter_hunt <- MCMCsummary(file_list[[1]], 'haz.kappa')
theta       <- MCMCsummary(file_list[[1]], 'theta')

winter_nat1 <- cbind.data.frame(winter_nat , expand.grid(age = c('HY','AHY'),species = c('CANVASBACK','REDHEAD'), sex = c('F','M'), year = 1961:2019, type = 'Winter:Natural'))
summer_nat1 <- cbind.data.frame(summer_nat , expand.grid(species = c('CANVASBACK','REDHEAD'), sex = c('F','M'), year = 1961:2019, type = 'Summer:Natural'))
winter_hunt1<- cbind.data.frame(winter_hunt, expand.grid(age = c('HY','AHY','AHY:I','AHY:W'),species = c('CANVASBACK','REDHEAD'), sex = c('F','M'), year = 1961:2019, type = 'Winter:Harvest'))
theta1       <- cbind.data.frame(theta , expand.grid(species = c('CANVASBACK','REDHEAD'), sex = c('F','M'), year = 1961:2019, type = 'Production'))

N_Year <- rbind.data.frame(cbind.data.frame(year = 1961:2019, species = 'CANVASBACK', N = CANS),cbind.data.frame(year = 1961:2019, species = 'REDHEAD', N = REDS))

salin_d <- rbind(cbind.data.frame(year = 1961:2019, species = 'CANVASBACK', salinity = c(exp(sal_spl[,1] * std_sal[1] + mean_sal[1]))),
                 cbind.data.frame(year = 1961:2019, species = 'REDHEAD', salinity = c(exp(sal_spl[,2] * std_sal[2] + mean_sal[2]))))

winter_nat1  <- left_join(winter_nat1 , salin_d)

summer_nat1  <- left_join(summer_nat1 , N_Year)
winter_nat1  <- left_join(winter_nat1 , N_Year)

winter_hunt2 <- subset(winter_hunt1, !grepl(":", age) )

all.ponds$year <- 1961:2019

ponds_l <- cbind.data.frame(year = all.ponds$year, Ponds = all.ponds$mean)
winter_hunt1  <- left_join(winter_hunt2  , N_Year)
winter_hunt1  <- left_join(winter_hunt1  , ponds_l )

summer_nat1  <- left_join(summer_nat1 , ponds_l)

salin_d_off <- subset(salin_d, year != 1961)
salin_d_off$year <- salin_d_off$year-1

theta1  <- left_join(theta1  , N_Year)
theta1  <- left_join(theta1  , ponds_l )
theta1  <- left_join(theta1  , salin_d_off)

wint2A <- ggplot(data = subset(winter_nat1 ,sex == 'F'), aes(x = salinity, y = mean, fill = age)) +
  geom_pointrange(aes(x = salinity, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21, linetype = 'dotted',size = 0.25,position = position_dodge(.5))+
  # geom_smooth(aes(x = salinity, y = mean), method = 'gam') +
  geom_ribbon(data = subset(sal_df, sex == 'F'), aes(x = Salinity, ymax = upper, ymin  = lower), color = NA, alpha =.25) +
  geom_smooth(data = subset(sal_df, sex == 'F'), aes(x = Salinity, y  = upper, color = age),linetype = 'solid') +
  geom_smooth(data = subset(sal_df, sex == 'F'), aes(x = Salinity, y  = lower, color = age), linetype = 'solid') +
  geom_line(data = subset(sal_df, sex == 'F'), aes(x = Salinity, y  = mean,  color = age),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free_x') + coord_cartesian(ylim = c(0,.8)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10) + labs(y = 'Natural winter mortality rate', x = 'Total Dissolved Solid (ppt)', color = 'Age', fill = 'Age')

covs2w <- MCMCpstr(mcmcList2, params = c('beta.surv'), type = 'chains')$beta.surv[,,2,]
covs2w <- melt(covs2w)
covs2w <- cbind.data.frame(covs2w, expand.grid(species = c('CANVASBACK','REDHEAD'), age =c('HY','AHY'), iter = 1:len))

labs1 <- covs2w %>%
  group_by(species, age) %>%
  summarize(total = n(),
            max = max(value) + .15 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs1$perc <- ifelse(labs1$perc < .5, 1 - labs1$perc,labs1$perc)


wint2B <-  ggplot() + 
  geom_violin(data = covs2w, aes(x = age, y = value, fill = age), linetype=0) +  
  stat_summary(data = covs2w, aes(x = age, y = value, group = age),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
wint2B 

wint2 <- ggarrange(wint2B,wint2A, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
wint2

#ggsave(wint2, file = 'Figures\\wint2.tiff',height = 8, width = 8, dpi = 320)

wint1A <-ggplot(data = subset(winter_nat1 ,sex == 'F'), aes(x = N/1000, y = mean, fill = age)) +
  geom_pointrange(aes(x = N/1000, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25,position = position_dodge(.5))+
  geom_ribbon(data = subset(nw_df, sex == 'F'), aes(x = N/1000, ymax = upper, ymin  = lower),color =  NA, alpha =.25) +
  geom_smooth(data = subset(nw_df, sex == 'F'),   aes(x = N/1000, y  = upper, color = age),linetype = 'solid') +
  geom_smooth(data = subset(nw_df, sex == 'F'),   aes(x = N/1000, y  = lower, color = age), linetype = 'solid') +
  geom_smooth(data = subset(nw_df, sex == 'F'),   aes(x = N/1000, y  = mean,  color = age),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free_x') + coord_cartesian(ylim = c(0,1.25)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10) + labs(y = 'Natural winter mortality rate', x = 'Conspecific Abundance (x1000)', color = 'Age', fill = 'Age')

covs1w <- MCMCpstr(mcmcList2, params = c('beta.surv'), type = 'chains')$beta.surv[,,1,]
covs1w <- melt(covs1w)
covs1w <- cbind.data.frame(covs1w, expand.grid(species = c('CANVASBACK','REDHEAD'), age =c('HY','AHY'), iter = 1:len))

labs2 <- covs1w %>%
  group_by(species, age) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs2$perc <- ifelse(labs2$perc < .5, 1 - labs2$perc,labs2$perc)

wint1B <-  ggplot() + 
  geom_violin(data = covs1w, aes(x = age, y = value, fill = age), linetype=0) +  
  stat_summary(data = covs1w, aes(x = age, y = value, group = age),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
wint1B 
wint1 <- ggarrange(wint1B,wint1A, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
wint1

#ggsave(wint1, file = 'Figures\\wint1.tiff',height = 8, width = 8, dpi = 320)

#####################################################################################
har1a <- ggplot(data = subset(winter_hunt1 ,sex == 'F'), aes(x = N/1000, y = mean, fill = age)) +
  geom_pointrange(aes(x = N/1000, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25,position = position_dodge(.5))+
  geom_ribbon(data = subset(nk_df, sex == 'F'), aes(x = N/1000, ymax = upper, ymin  = lower), color = NA, alpha =.25) +
  geom_smooth(data = subset(nk_df, sex == 'F'),   aes(x = N/1000, y  = upper, color = age),linetype = 'solid') +
  geom_smooth(data = subset(nk_df, sex == 'F'),   aes(x = N/1000, y  = lower, color = age), linetype = 'solid') +
  geom_smooth(data = subset(nk_df, sex == 'F'),   aes(x = N/1000, y  = mean,  color = age),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free') + 
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10) + labs(y = 'Harvest mortality hazard rate', x = 'Conspecific Abundance (x1000)', color = 'Age', fill = 'Age')

covs1h <- MCMCpstr(mcmcList2, params = c('beta.surv'), type = 'chains')$beta.surv[,,3,]
covs1h <- melt(covs1h)
covs1h <- cbind.data.frame(covs1h, expand.grid(species = c('CANVASBACK','REDHEAD'), age =c('HY','AHY'), iter = 1:len))

labs3 <- covs1h %>%
  group_by(species, age) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs3$perc <- ifelse(labs3$perc < .5, 1 - labs3$perc,labs3$perc)

har1b <-  ggplot() + 
  geom_violin(data = covs1h, aes(x = age, y = value, fill = age), linetype=0) +  
  stat_summary(data = covs1h, aes(x = age, y = value, group = age),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
har1b 

harv1 <- ggarrange(har1b,har1a, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
harv1

#ggsave(harv1, file = 'Figures\\harv1.tiff',height = 8, width = 8, dpi = 320)


har2a <-ggplot(data =  subset(winter_hunt1 ,sex == 'F'), aes(x = Ponds, y = mean, fill = age)) +
  geom_pointrange(aes(x = Ponds, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25,position = position_dodge(.5))+
  geom_ribbon(data = subset(pondk_df, sex == 'F'), aes(x = Ponds, ymax = upper, ymin  = lower), color = NA, alpha =.25) +
  geom_smooth(data = subset(pondk_df, sex == 'F'),   aes(x = Ponds, y  = upper, color = age),linetype = 'solid') +
  geom_smooth(data = subset(pondk_df, sex == 'F'),   aes(x = Ponds, y  = lower, color = age), linetype = 'solid') +
  geom_smooth(data = subset(pondk_df, sex == 'F'),   aes(x = Ponds, y  = mean,  color = age),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free') +  coord_cartesian(xlim = c(3000,6500)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10) + labs(y = 'Harvest mortality hazard rate', x = 'Number of ponds in PPR (x1000)', color = 'Age', fill = 'Age')

covs2h <- MCMCpstr(mcmcList2, params = c('beta.surv'), type = 'chains')$beta.surv[,,4,]
covs2h <- melt(covs2h)
covs2h <- cbind.data.frame(covs2h, expand.grid(species = c('CANVASBACK','REDHEAD'), age =c('HY','AHY'), iter = 1:len))

labs4 <- covs2h %>%
  group_by(species, age) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs4$perc <- ifelse(labs4$perc < .5, 1 - labs4$perc,labs4$perc)

har2b <-  ggplot() + 
  geom_violin(data = covs2h, aes(x = age, y = value, fill = age), linetype=0) +  
  stat_summary(data = covs2h, aes(x = age, y = value, group = age),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
har2b 

harv2 <- ggarrange(har2b,har2a, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
harv2 

#ggsave(harv2, file = 'Figures\\harv2.tiff',height = 8, width = 8, dpi = 320)
################################################################################

sum1a <- ggplot(data = subset(summer_nat1, sex == 'F'), aes(x = Ponds, y = mean, fill = species)) +
  geom_pointrange(aes(x = Ponds, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25,)+
  geom_ribbon(data = subset(pond_df, sex == 'F'),aes(x = Ponds, ymax = upper, ymin  = lower), color = NA, alpha =.25) +
  geom_smooth(data = subset(pond_df, sex == 'F'),   aes(x = Ponds, y  = upper, color = species),linetype = 'solid') +
  geom_smooth(data = subset(pond_df, sex == 'F'),   aes(x = Ponds, y  = lower, color = species), linetype = 'solid') +
  geom_smooth(data = subset(pond_df, sex == 'F'),   aes(x = Ponds, y  = mean,  color = species),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free') +  coord_cartesian(xlim = c(3000,6500)) +
  scale_fill_see_d()+ scale_color_see_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10) + labs(y = 'Summer mortality hazard rate', x = 'Number of ponds in PPR (x1000)', color = 'Species', fill = 'Species')

covs1s <- MCMCpstr(mcmcList2, params = c('beta.surv3'), type = 'chains')$beta.surv3[,2,]
covs1s <- melt(covs1s)
covs1s <- cbind.data.frame(covs1s, expand.grid(species = c('CANVASBACK','REDHEAD'), iter = 1:len))

labs5 <- covs1s %>%
  group_by(species) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs5$perc <- ifelse(labs5$perc < .5, 1 - labs5$perc,labs5$perc)

sum1b <-  ggplot() + 
  geom_violin(data = covs1s, aes(x = species, y = value, fill = species), linetype=0) +  
  stat_summary(data = covs1s, aes(x = species, y = value, group = species),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
sum1b 



summer1 <- ggarrange(sum1b,sum1a, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
summer1

#ggsave(summer1, file = 'Figures\\summer1.tiff',height = 8, width = 8, dpi = 320)

sum2a <-ggplot(data = subset(summer_nat1,sex == 'F'), aes(x = N/1000, y = mean, fill = species)) +
  geom_pointrange(aes(x = N/1000, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = species), shape = 21,linetype = 'dotted',size = 0.25,)+
  geom_ribbon(data = subset(ns_df, sex == 'F'), aes(x = N/1000, ymax = upper, ymin  = lower, fill = species), color = NA, alpha =.25) +
  geom_smooth(data = subset(ns_df, sex == 'F'),   aes(x = N/1000, y  = upper, color = species),linetype = 'solid') +
  geom_smooth(data = subset(ns_df, sex == 'F'),   aes(x = N/1000, y  = lower, color = species), linetype = 'solid') +
  geom_smooth(data = subset(ns_df, sex == 'F'),   aes(x = N/1000, y  = mean,  color = species),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free') + 
  scale_fill_see_d()+ scale_color_see_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10) + labs(y = 'Summer mortality hazard rate', x = 'Conspecific Abundance (x1000)', color = 'Species', fill = 'Species')

covs2s <- MCMCpstr(mcmcList2, params = c('beta.surv3'), type = 'chains')$beta.surv3[,1,]
covs2s <- melt(covs2s)
covs2s <- cbind.data.frame(covs2s, expand.grid(species = c('CANVASBACK','REDHEAD'), iter = 1:len))

labs6 <- covs2s %>%
  group_by(species) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs6$perc <- ifelse(labs6$perc < .5, 1 - labs6$perc,labs6$perc)

sum2b <-  ggplot() + 
  geom_violin(data = covs2s, aes(x = species, y = value, fill = species), linetype=0) +  
  stat_summary(data = covs2s, aes(x = species, y = value, group = species),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
sum2b 

summer2 <- ggarrange(sum2b,sum2a, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
summer2

#ggsave(summer2, file = 'Figures\\summer2.tiff',height = 8, width = 8, dpi = 320)

####################################################################################################

prod2a <- ggplot(data = subset(theta1, sex == 'F'), aes(x = Ponds, y = mean, fill = species)) +
  geom_pointrange(aes(x = Ponds, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted', size = .25)+
  geom_ribbon(data = subset(prod.p_df, sex == 'F'),aes(x = Ponds, ymax = upper, ymin  = lower), color = NA, alpha =.25) +
  geom_smooth(data =   subset(prod.p_df, sex == 'F'),   aes(x = Ponds, y  = upper, color = species),linetype = 'solid') +
  geom_smooth(data =   subset(prod.p_df, sex == 'F'),   aes(x = Ponds, y  = lower, color = species), linetype = 'solid') +
  geom_smooth(data =   subset(prod.p_df, sex == 'F'),   aes(x = Ponds, y  = mean,  color = species), linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free') + coord_cartesian(ylim = c(0,1), xlim = c(3500,7000)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend = 'top', legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) + labs(y = 'Duckling Production', x = 'Number of ponds in PPR (x1000)', color = 'Species', fill = 'Species')

covs2p <- MCMCpstr(mcmcList2, params = c('beta.np'), type = 'chains')$beta.np[,2,]
covs2p <- melt(covs2p)
covs2p <- cbind.data.frame(covs2p, expand.grid(species = c('CANVASBACK','REDHEAD'), iter = 1:len))

labs7 <- covs2p %>%
  group_by(species) %>%
  summarize(total = n(),
            max = max(value) + .1 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs7$perc <- ifelse(labs7$perc < .5, 1 - labs7$perc,labs7$perc)

prod2b <-  ggplot() + 
  geom_violin(data = covs2p, aes(x = species, y = value, fill = species), linetype=0) +  
  stat_summary(data = covs2p, aes(x = species, y = value, group = species),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
prod2b 


product2 <- ggarrange(prod2b,prod2a, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
product2

#ggsave(product2, file = 'Figures\\product2.tiff',height = 8, width = 8, dpi = 320)

prod1a <-ggplot(data = subset(theta1, sex == 'F'), aes(x = N/1000, y = mean, fill = species)) +
  geom_pointrange(aes(x = N/1000, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25)+
  geom_ribbon(data= subset(prod.n_df, sex == 'F'),aes(x = N/1000, ymax = upper, ymin  = lower, color = species), alpha =.25) +
  geom_smooth(data = subset(prod.n_df, sex == 'F'),   aes(x = N/1000, y  = upper, color = species),linetype = 'solid') +
  geom_smooth(data = subset(prod.n_df, sex == 'F'),   aes(x = N/1000, y  = lower, color = species), linetype = 'solid') +
  geom_smooth(data = subset(prod.n_df, sex == 'F'),   aes(x = N/1000, y  = mean,  color = species),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free') + coord_cartesian(ylim = c(0,1)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend = 'top', legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) + labs(y = 'Duckling Production', x = 'Conspecific Abundance (x1000)',  color = 'Species', fill = 'Species')

covs1p <- MCMCpstr(mcmcList2, params = c('beta.np'), type = 'chains')$beta.np[,1,]
covs1p <- melt(covs1p)
covs1p <- cbind.data.frame(covs1p, expand.grid(species = c('CANVASBACK','REDHEAD'), iter = 1:len))

labs8 <- covs1p %>%
  group_by(species) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs8$perc <- ifelse(labs8$perc < .5, 1 - labs8$perc,labs8$perc)

labs8$max <- c(.275,.05)

prod1b <-  ggplot() + 
  geom_violin(data = covs1p, aes(x = species, y = value, fill = species), linetype=0) +  
  stat_summary(data = covs1p, aes(x = species, y = value, group = species),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
prod1b 


product1 <- ggarrange(prod1b,prod1a, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
product1

#ggsave(product1, file = 'Figures\\product1.tiff',height = 8, width = 8, dpi = 320)


prod0a <-ggplot(data = subset(theta1, sex == 'F'), aes(x = salinity, y = mean, fill = species)) +
  geom_pointrange(aes(x = salinity, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25)+
  geom_ribbon(data= subset(prod.ds_df, sex == 'F'),aes(x = Salinity, ymax = upper, ymin  = lower, color = species), alpha =.25) +
  geom_smooth(data = subset(prod.ds_df, sex == 'F'),   aes(x = Salinity, y  = upper, color = species),linetype = 'solid') +
  geom_smooth(data = subset(prod.ds_df, sex == 'F'),   aes(x = Salinity, y  = lower, color = species), linetype = 'solid') +
  geom_smooth(data = subset(prod.ds_df, sex == 'F'),   aes(x = Salinity, y  = mean,  color = species),linetype = 'dashed', size = .5) +
  facet_wrap(~species, scales= 'free') + coord_cartesian(ylim = c(0,1)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend = 'top', legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) + labs(y = 'Duckling Production', x = 'Estimated Total Dissolved Solids\non Wintering Grounds (mg/L)',  color = 'Species', fill = 'Species')

covs0p <- MCMCpstr(mcmcList2, params = c('beta.np'), type = 'chains')$beta.np[,3,]
covs0p <- melt(covs0p)
covs0p <- cbind.data.frame(covs0p, expand.grid(species = c('CANVASBACK','REDHEAD'), iter = 1:len))

labs9 <- covs0p %>%
  group_by(species) %>%
  summarize(total = n(),
            max = max(value) + .25 * max(value),
            above = sum(value > 0),
            perc = above/total)

labs9$perc <- ifelse(labs9$perc < .5, 1 - labs9$perc,labs9$perc)


prod0b <-  ggplot() + 
  geom_violin(data = covs0p, aes(x = species, y = value, fill = species), linetype=0) +  
  stat_summary(data = covs0p, aes(x = species, y = value, group = species),color = 'white', fun.data=median_hdi, position = position_dodge(.9))  +
  facet_wrap(~species, scales = 'free') + geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_flat_d() +  theme_modern(legend.position = 'top', axis.title.size = 11,legend.text.size = 10, legend.title.size = 10,) + 
  labs(y = 'Posterior distribution\nfor slope coefficients', x = '', color = 'Age', fill = 'Age')
prod0b 


product0 <- ggarrange(prod0b,prod0a, common.legend = TRUE, labels = c('A','B'), nrow = 2, heights = c(1,2))
product0

#ggsave(product0, file = 'Figures\\product0.tiff',height = 8, width = 8, dpi = 320)

congen <- MCMCsummary(mcmcList1, params = 'tall')
congen1 <- cbind.data.frame(congen, expand.grid(species = c('CANVASBACK','REDHEAD'), sex = c('F','M'), year = 1961:2019, type = 'Production'))


RN_Year <- rbind.data.frame(cbind.data.frame(year = 1961:2019, species = 'CANVASBACK', RN = REDS),cbind.data.frame(year = 1961:2019, species = 'REDHEAD', RN = CANS))
congen1  <- left_join(congen1 , RN_Year)
theta2 <- left_join(theta1, RN_Year)

plot_2 <- ggplot(data = subset(theta2, sex == 'F'), aes(x = year, y = mean,fill = species)) +
  geom_pointrange(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25)+
  geom_ribbon(data =subset(congen1,sex == 'F'),aes(x = year, ymax = `97.5%`, ymin  = `2.5%`), alpha =.25) +
  geom_line(data =subset(congen1,sex == 'F'),aes(x = year, y = `50%`, color = species), linetype = 'solid', size =.5) +
  facet_wrap(~species, scales= 'free') +  coord_cartesian(ylim = c(0,1)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend = 'top', legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) + labs(y = 'Duckling Production', x = 'Year',fill = 'Species',color = 'Species')

plot_1 <- ggplot(data = subset(theta2, sex == 'F'), aes(x = RN/1000, y = mean, fill = species)) +
  geom_pointrange(aes(x = RN/1000, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted',size = 0.25)+
  geom_ribbon(data =subset(congen1,sex == 'F'),aes(x = RN/1000, ymax = `97.5%`, ymin  = `2.5%`), alpha =.25) +
  geom_line(data =subset(congen1,sex == 'F'),aes(x = RN/1000, y = `50%`, color = species), linetype = 'solid', size =.5) +
  facet_wrap(~species, scales= 'free') + coord_cartesian(ylim = c(0,1)) +
  scale_fill_flat_d()+ scale_color_flat_d()+ theme_modern(legend = 'top', legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) + labs(y = 'Duckling Production', x = 'Interspecific Abundance (x1000)', fill = 'Species',color = 'Species')


product3 <- ggarrange(plot_1,plot_2, common.legend = TRUE, labels = c('A','B'), nrow = 2)
product3
#ggsave(product3, file = 'Figures\\product3.tiff',height = 8, width = 8, dpi = 320)
library(ggpubr)



SJ <- MCMCsummary(mcmcList1, 'SJ')
SA <- MCMCsummary(mcmcList1, 'SA')


RS <- MCMCsummary(mcmcList1, 'relative.summer')
NW <- MCMCsummary(mcmcList1, 'natural.winter')
AK <- MCMCsummary(mcmcList1, 'annual.kappa')
WK <- MCMCsummary(mcmcList1, 'kappa')

RS1 <- cbind.data.frame(RS, expand.grid(age = c('HY','AHY'), species = c('Canvasback','Redhead'), sex = c('Female','Male'), year = 1961:2019))
NW1 <- cbind.data.frame(NW, expand.grid(age = c('HY','AHY'), species = c('Canvasback','Redhead'), sex = c('Female','Male'), year = 1961:2019), Type = 'Natural')
AK1 <- cbind.data.frame(AK, expand.grid(age = c('HY','AHY'), species = c('Canvasback','Redhead'), sex = c('Female','Male'), year = 1961:2019), Type = 'Harvest')
WK1 <- cbind.data.frame(WK, expand.grid(age = c('HY','AHY','AHY:I','AHY:W'), species = c('Canvasback','Redhead'), sex = c('Female','Male'), year = 1961:2019), Type = 'Harvest')
WK1 <- subset(WK1, !grepl(":",age))
winter <- rbind.data.frame(NW1,WK1)


winter.mort1 <- ggplot(data =winter,aes(x = year, y = mean, fill = interaction(Type,species) )) +
  geom_pointrange2(aes(x = year, y = mean, ymin =`2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = interaction(Type,species)), method = 'gam', se = FALSE)+
  facet_wrap(age~ sex, ncol = 2 ) + coord_cartesian(ylim = c(0,1)) +
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top') + labs(y = 'Probability of cause-specific winter mortality', fill= 'Cause/Species', color = 'Cause/Species')
winter.mort1

#ggsave(winter.mort1, file = 'Figures\\winter.mort1.tiff', dpi = 320, height = 8, width = 12)

annuals <- rbind(cbind.data.frame(SJ , expand.grid(species =c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019,age = 'HY')),
                 cbind.data.frame(SA , expand.grid(species =c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019,age = 'AHY')))                                     


annuals <- left_join(annuals, salin_d)


A_C <- subset(annuals, year < 2019 & species == 'CANVASBACK')
A_R <- subset(annuals, year < 2019 & species != 'CANVASBACK')

colnames(A_C) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c', 'species.c', 'sex', 'year', 'age','salinity.c')
colnames(A_R) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r', 'species.r', 'sex', 'year', 'age','salinity.r')

A_ <- left_join(A_C,A_R)

splot <- ggplot(data = A_, aes(x = mean.c, y = mean.r)) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), alpha = .5, linetype = 'dotted')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), alpha = .5, linetype = 'dotted')+  
  geom_point(aes(x = mean.c, y = mean.r, color = sex), size = 2, alpha = .75) +
  geom_smooth(aes(group = sex), color = 'black', se = FALSE, method = 'lm', formula= y ~x) +
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1,linetype = 'dashed') + facet_wrap(~age) +
  labs(y = 'Annual Survival\n Redheads',x = 'Annual Survival\n Canvasvacks', color= 'Sex') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)
splot

#ggsave(splot, file = 'surv_plot.tiff', width =12, height = 6, dpi = 320)

annual_surv <-  ggplot(data = annuals , aes(x = year, y = mean, fill = species)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21, size = .5, linetype = 'dotted', position = position_dodge2(.5))+
  #geom_smooth(aes(color = age), alpha = .25, se = FALSE)+
  facet_wrap( age ~ sex, scales = 'free') + coord_cartesian(ylim = c(0,1))+
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top') + labs(y = 'Survival probabilities', fill = 'Age', color = 'Age', x = '')


#ggsave(annual_surv, file = 'Figures\\annual_surv.tiff',height = 8, width = 8, dpi = 320)

HK <- MCMCsummary(mcmcList1, 'haz.kappa')
HW <- MCMCsummary(mcmcList1, 'haz.winter')
HS <- MCMCsummary(mcmcList1, 'haz.summer')

annual.kappa <- MCMCsummary(mcmcList1, 'annual.kappa')
K <- MCMCsummary(mcmcList1, 'kappa')
s.winter <- MCMCsummary(mcmcList1, 's.winter')
s.summer <- MCMCsummary(mcmcList1, 's.summer')
natural.winter <- MCMCsummary(mcmcList1, 'natural.winter')
relative.summer <- MCMCsummary(mcmcList1, 'relative.summer')



HK_risk <- cbind.data.frame(HK, expand.grid(age = c('HY:Direct','AHY:Direct','AHY:Indirect','AHY:Winter'), species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
KR      <- cbind.data.frame(K, expand.grid(age = c('HY:Direct','AHY:Direct','AHY:Indirect','AHY:Winter'), species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
AKR     <- cbind.data.frame(annual.kappa, expand.grid(age = c('HY:Direct','AHY:Direct'), species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
HW_risk <- cbind.data.frame(HW, expand.grid(age = c('HY','AHY'),species =c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
WS      <- cbind.data.frame(s.winter, expand.grid(age = c('HY','AHY'),species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
SS      <- cbind.data.frame(s.summer, expand.grid(age = c('HY','AHY'),species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
RS      <- cbind.data.frame(relative.summer, expand.grid(age = c('HY','AHY'),species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
HS_risk <- cbind.data.frame(HS, expand.grid(species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
natural.winter <- cbind.data.frame(natural.winter, expand.grid(age = c('HY','AHY'),species =c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))

summer.mort1 <- ggplot(data =RS,aes(x = year, y = mean, fill = sex)) +
  geom_pointrange2(aes(x = year, y = mean, ymin =`2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = sex))+
  facet_wrap(age~ species ) + coord_cartesian(ylim = c(0,.6)) +
  scale_fill_flat_d()+scale_color_flat_d()+ 
  theme_modern(legend.position = 'top', legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) +
  labs(y = 'Proportion of annual mortality that occurred during summer', fill= 'Sex', color = 'Sex', x = '')
summer.mort1

#ggsave(summer.mort1 , file = 'Figures\\summer.mort .tiff',height = 8, width = 8, dpi = 320)

kr_sub <- subset(AKR, grepl(':Direct', age))

harvest.mort <- ggplot(data = kr_sub,aes(x = year, y = mean, fill = sex)) +
  geom_pointrange2(aes(x = year, y = mean, ymin =`2.5%`, ymax = `97.5%`, fill = sex), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = sex))+
  facet_wrap(age~ species ) + coord_cartesian(ylim = c(0,1)) +
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) + 
  labs(y = 'Probability of being harvested (kill rate)', fill= 'Sex', color = 'Sex', x = '')
harvest.mort


#ggsave(harvest.mort , file = 'Figures\\harvest.mort .tiff',height = 8, width = 8, dpi = 320)

winter.mort1 <- ggplot(data =natural.winter,aes(x = year, y = mean, fill = sex)) +
  geom_pointrange2(aes(x = year, y = mean, ymin =`2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = sex))+
  facet_wrap(age~ species ) + coord_cartesian(ylim = c(0,1)) +
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 10, axis.title.size = 10, axis.text.size = 10) +
  labs(y = 'Probability of natural winter mortality', fill= 'Sex', color = 'Sex', x= '')
winter.mort1

#ggsave(winter.mort1 , file = 'Figures\\winter.mort .tiff',height = 8, width = 8, dpi = 320)


HS_C <- subset(HS_risk, year < 2019 & species == 'CANVASBACK')
HS_R <- subset(HS_risk, year < 2019 & species != 'CANVASBACK')

colnames(HS_C) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c ','species.c','sex','year')
colnames(HS_R) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r ','species.r','sex','year')

HS_ <- left_join(HS_C,HS_R)

HK_risk <- subset(HK_risk, grepl('Direct', age))

HK_C <- subset(HK_risk, year < 2019 & species == 'CANVASBACK')
HK_R <- subset(HK_risk, year < 2019 & species != 'CANVASBACK')

colnames(HK_C) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c ','age',  'species.c','sex','year')
colnames(HK_R) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r ','age',  'species.r','sex','year')

HK_ <- left_join(HK_C,HK_R)


HW_C <- subset(HW_risk, year < 2019 & species == 'CANVASBACK')
HW_R <- subset(HW_risk, year < 2019 & species != 'CANVASBACK')

colnames(HW_C) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c ','age',  'species.c','sex','year')
colnames(HW_R) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r ','age',  'species.r','sex','year')

HW_ <- left_join(HW_C,HW_R)


M1F <- ggplot(data = subset(HW_, sex == 'Female') , aes(x = mean.c, y = mean.r)) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), alpha = .5, linetype = 'dotted')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), alpha = .5, linetype = 'dotted')+  
  geom_point(aes(x = mean.c, y = mean.r, color = age), size = 2, alpha = .75) +
  geom_smooth(aes(group = age), method = 'lm', formula= y ~x, size = 1, color = 'black', se = FALSE) + coord_cartesian(ylim =  c(0,.75), xlim = c(0,.75))+
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + #facet_wrap(~age, scales = 'free') +
  labs(y = 'Risk of Natural Winter Mortality\n Female Redheads',x = 'Risk of Natural Winter Mortality\n Female Canvasbacks', color = 'Age') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)

M1M <- ggplot(data = subset(HW_, sex == 'Male') , aes(x = mean.c, y = mean.r)) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), alpha = .5, linetype = 'dotted')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), alpha = .5, linetype = 'dotted')+  
  geom_point(aes(x = mean.c, y = mean.r, color = age), size = 2, alpha = .75) +
  geom_smooth(aes(group = age), method = 'lm', formula= y ~x, size = 1, color = 'black', se = FALSE) + coord_cartesian(ylim =  c(0,.75), xlim = c(0,.75))+
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + #facet_wrap(~age, scales = 'free') +
  labs(y = 'Risk of Natural Winter Mortality\n Male Redheads',x = 'Risk of Natural Winter Mortality\n Male Canvasbacks', color = 'Age') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)


HK_$age <- recode(HK_$age, 'HY:Direct' = 'HY', 'AHY:Direct' = 'AHY')

M2F <- ggplot(data = subset(HK_, sex == 'Female') , aes(x = mean.c, y = mean.r)) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), alpha = .5, linetype = 'dotted')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), alpha = .5, linetype = 'dotted')+  
  geom_point(aes(x = mean.c, y = mean.r, color = age), size = 2, alpha = .75) +
  geom_smooth(aes(group = age), method = 'lm', formula= y ~x, size = 1, color = 'black', se = FALSE) + coord_cartesian(ylim =  c(0,.75), xlim = c(0,.75))+
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + #facet_wrap(~age, scales = 'free') +
  labs(y = 'Risk of Harvest Mortality\n Female Redheads',x = 'Risk of Harvest Mortality\n Female Canvasbacks', color = 'Age') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)

M2M <- ggplot(data = subset(HK_, sex == 'Male') , aes(x = mean.c, y = mean.r)) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), alpha = .5, linetype = 'dotted')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), alpha = .5, linetype = 'dotted')+  
  geom_point(aes(x = mean.c, y = mean.r, color = age), size = 2, alpha = .75) +
  geom_smooth(aes(group = age), method = 'lm', formula= y ~x, size = 1, color = 'black', se = FALSE) + coord_cartesian(ylim =  c(0,.75), xlim = c(0,.75))+
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + #facet_wrap(~age, scales = 'free') +
  labs(y = 'Risk of Harvest Mortality\n Male Redheads',x = 'Risk of Harvest Mortality\n Male Canvasbacks', color = 'Age') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)

M3F <- ggplot(data = subset(HS_, sex == 'Female'), aes(x = mean.c, y = mean.r)) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), alpha = .5, linetype = 'dotted')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), alpha = .5, linetype = 'dotted')+  
  geom_point(aes(x = mean.c, y = mean.r),color = '#9b59b6', size = 2, alpha = .75) +
  geom_smooth(method = 'lm', formula= y ~x, size = 1,se = FALSE, color = 'black') +
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + coord_cartesian(ylim =  c(0,.75), xlim = c(0,.75))+
  labs(y = 'Risk of Natural Summer Mortality\n Female Redheads',x = 'Risk of Natural Summer Mortality\n Female Canvasbacks') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)

M3M <- ggplot(data = subset(HS_, sex == 'Male'), aes(x = mean.c, y = mean.r)) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), alpha = .5, linetype = 'dotted')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), alpha = .5, linetype = 'dotted')+  
  geom_point(aes(x = mean.c, y = mean.r),color = '#9b59b6', size = 2, alpha = .75) +
  geom_smooth(method = 'lm', formula= y ~x, size = 1,color = 'black',se = FALSE) +
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + coord_cartesian(ylim =  c(0,.75), xlim = c(0,.75))+
  labs(y = 'Risk of Natural Summer Mortality\n Male Redheads',x = 'Risk of Natural Summer Mortality\n Male Canvasbacks') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)


risks <- ggarrange(M2F,M1F, M3F,M2M,M1M, M3M, labels = c('A','B','C','D','E','F'),common.legend = TRUE, ncol = 3, nrow = 2)
#ggsave(risks, file = 'risk_plot.tiff', width =12, height = 9, dpi = 320)


S1 <- ggplot(data = subset(HS_risk, year < 2019 & sex == 'Female'),aes(x = year, y = mean)) +
  geom_smooth(aes(color = species), alpha = .25)+
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = species), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  theme_modern(legend.position = 'top') + 
  scale_fill_flat_d() +  scale_color_flat_d() +  labs(y = 'Natural summer mortality risk', fill= 'Species', color = 'Species', x ='')

HK_risk <- subset(HK_risk, grepl('Direct', age))

H1 <- ggplot(data = subset(HK_risk ,sex == 'Female'), aes(x = year, y = mean, fill = age)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = age))+
  facet_wrap(~ species , ncol = 1) +  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_flat_d(labels = c('HY','AHY'))+scale_color_flat_d(labels = c('HY','AHY'))+ theme_modern(legend.position = 'top') + labs(y = 'Harvest mortality risk', fill = 'Age', color = "Age", x ='')

W1 <- ggplot(data = subset(HW_risk ,sex == 'Female'), aes(x = year, y = mean, fill = age)) +
  geom_pointrange2(aes(x = year, y = mean, ymin = `2.5%`, ymax = `97.5%`), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = age))+
  facet_wrap(~ species , ncol = 1) +  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top') + labs(y = 'Natural winter mortality risk', fill = 'Age', color = "Age", x ='')


test_risk <- ggarrange(W1,H1,S1, labels = c('A','B','C'), ncol = 3, common.legend = TRUE)

#ggsave(test_risk , file = 'Figures\\risks .tiff',height = 8, width = 12, dpi = 320)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df.SF      <- MCMCsummary(mcmcList1, 'SF')
df.SM      <- MCMCsummary(mcmcList1, 'SM')
df.PF      <- MCMCsummary(mcmcList1, 'PF')
df.PM      <- MCMCsummary(mcmcList1, 'PM')
df.N       <- MCMCsummary(mcmcList1, 'N')
df.lambda  <- MCMCsummary(mcmcList1, 'lambda')

df1.SF     <- cbind.data.frame(df.SF     , expand.grid(species = c('CANVASBACK','REDHEAD'), year = 1961:2019), type = 'Surviving Female')
df1.SM     <- cbind.data.frame(df.SM     , expand.grid(species = c('CANVASBACK','REDHEAD'), year = 1961:2019), type = 'Surviving Male')
df1.PF     <- cbind.data.frame(df.PF     , expand.grid(species = c('CANVASBACK','REDHEAD'), year = 1961:2019), type = 'Recruiting Female')
df1.PM     <- cbind.data.frame(df.PM     , expand.grid(species = c('CANVASBACK','REDHEAD'), year = 1961:2019), type = 'Recruiting Male')
df1.N      <- cbind.data.frame(df.N      , expand.grid(species = c('CANVASBACK','REDHEAD'), year = 1961:2019))
df1.lambda <- cbind.data.frame(df.lambda , expand.grid(species = c('CANVASBACK','REDHEAD'), year = 1961:2019))

df.age_rat <- MCMCsummary(mcmcList1, 'age.ratio')

df1.age_rat <- cbind.data.frame(df.age_rat, expand.grid(species = c('CANVASBACK','REDHEAD'), sex = c('Female','Male'), year = 1961:2019))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

df_types <- rbind.data.frame(df1.SF,df1.SM,df1.PF,df1.PM)

div <- cbind.data.frame(species = df1.lambda$species, year = df1.lambda$year, div = df1.lambda$mean)

df_types1 <- left_join(df_types, div)

dc<- ggplot(data = df_types1,aes(x = year, y = (mean/div))) +
  geom_pointrange2(aes(x = year, y = (mean/div), ymin = (`2.5%`/div), ymax = (`97.5%`/div), fill = type), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = type), se = FALSE, method = 'lm', formula = y~x)+ coord_cartesian(ylim = c(0,0.75))+
  facet_wrap( ~ species , scales = 'free') + 
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top') + labs(y = 'Proportional contribution\n to population growth',x = '', fill = 'Contributor', color  = 'Contributor')
dc


df1.type_r <- subset(df_types1, species == 'REDHEAD')
df1.type_c <- subset(df_types1, species != 'REDHEAD')

colnames(df1.type_c) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c','species.c','year','type','div.c')
colnames(df1.type_r) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r','species.r','year','type','div.r')

types_ <- left_join(df1.type_c,df1.type_r)

c.lam_com <- ggplot(data = types_, aes(x = mean.c/div.c, y = mean.r/div.r)) +
  geom_point(aes(x = mean.c/div.c, y = mean.r/div.r, color = type), size = 2, alpha = .75) +
  geom_smooth(aes(color = type), method = 'lm', formula= y ~x, size = 1.5) +
  geom_errorbar(aes(y = mean.c/div.c, ymin=low.r/div.r, ymax=upp.r/div.r, color= type), alpha = .5, linetype = 'dashed')+  
  geom_errorbarh(aes(y = mean.r/div.r,xmin=low.c/div.c, xmax=upp.c/div.c, color= type), alpha = .5, linetype = 'dashed')+  
  scale_color_flat_d(guide = guide_legend(nrow = 2)) + geom_abline(intercept = 0, slope = 1) +
  labs(x = 'Proportional Composition of\nPopulation Growth Rate\n Canvasbacks', y = 'Proportional Composition of\nPopulation Growth Rate\n Redheads', color = 'Contributor') +
  theme_modern(legend.position = 'top', legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)

c.lam_com 


dc1<- ggplot(data = df_types1,aes(x = year, y = (mean))) +
  geom_pointrange2(aes(x = year, y = (mean), ymin = (`2.5%`), ymax = (`97.5%`), fill = type), shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = type), se = FALSE, method = 'lm', formula = y~x)+ coord_cartesian(ylim = c(0,0.75))+
  facet_wrap( ~ species , scales = 'free') + 
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top') + labs(y = 'Relative numerical contribution\n to population growth',x = '', fill = 'Contributor', color  = 'Contributor')
dc1


df1.lambda_r <- subset(df1.lambda, species == 'REDHEAD')
df1.lambda_c <- subset(df1.lambda, species != 'REDHEAD')

colnames(df1.lambda_c) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c','species.c','year')
colnames(df1.lambda_r) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r','species.r','year')

LAM_ <- left_join(df1.lambda_c,df1.lambda_r)

lam_com <- ggplot(data = LAM_, aes(x = mean.c, y = mean.r)) +
  geom_point(aes(x = mean.c, y = mean.r), color = '#16a085', size = 2, alpha = .75) +
  geom_smooth(method = 'lm', formula= y ~x,color = '#16a085',size = 1.5) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r), color = '#16a085', alpha = .5, linetype = 'dashed')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c), color = '#16a085', alpha = .5, linetype = 'dashed')+  
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1) +
  labs(x = 'Population Growth Rate\n Canvasbacks', y = 'Population Growth Rate\n Redheads\n') +
  theme_modern(legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)
lam_com

lc <- ggplot(data = df1.lambda,aes(x = year, y = (mean))) +
  geom_pointrange2(aes(x = year, y = (mean), ymin = (`2.5%`), ymax = (`97.5%`)), fill = 'black', shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(color = 'red', se = FALSE, method = 'lm', formula = y~x)+ coord_cartesian(ylim = c(0.5,1.65))+
  geom_hline(yintercept = 1.0, color = 'black')+
  facet_wrap( ~ species,scales = 'free' ) + 
  theme_modern(legend.position = 'none') + labs(y = 'Rate of population growth\n', x = '')



growth <- ggarrange(dc,lc, common.legend = FALSE, labels = c('A','B'), nrow = 2)
growth

pc<- ggplot(data = df1.age_rat,aes(x = year, y = (mean))) +
  geom_pointrange2(aes(x = year, y = (mean), ymin = (`2.5%`), ymax = (`97.5%`), fill = sex), size = .65, shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  geom_smooth(aes(color = sex), se = FALSE, method = 'gam')+ 
  facet_wrap( ~ species , scales = 'free') +
  geom_hline(yintercept = 0, color = NA) +
  geom_hline(yintercept = 2.5, color = NA) +
  scale_fill_flat_d()+scale_color_flat_d()+ theme_modern(legend.position = 'top') + labs(y = 'Per-capita duck production',x = '', fill = 'Sex', color  = 'Sex')


df1.age_rat_r <- subset(df1.age_rat, species == 'REDHEAD')
df1.age_rat_c <- subset(df1.age_rat, species != 'REDHEAD')

colnames(df1.age_rat_c) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c','species.c','sex','year')
colnames(df1.age_rat_r) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r','species.r','sex','year')

SR_ <- left_join(df1.age_rat_c,df1.age_rat_r)


prod.comp <- ggplot(data = SR_, aes(x = mean.c, y = mean.r)) +
  geom_point(aes(x = mean.c, y = mean.r, color = sex), size = 2, alpha = .75) +
  geom_smooth(aes(color = sex), method = 'lm', formula= y ~x, size = 1.5) +
  geom_errorbar(aes(y = mean.c, ymin=low.r, ymax=upp.r, color= sex), alpha = .5, linetype = 'dashed')+  
  geom_errorbarh(aes(y = mean.r,xmin=low.c, xmax=upp.c, color= sex), alpha = .5, linetype = 'dashed')+  
  scale_color_flat_d(guide = guide_legend(nrow = 2)) + geom_abline(intercept = 0, slope = 1) +
  labs(x = 'Sex-Specific Per-Capita\nDuck Production\n Canvasbacks', y = 'Sex-Specific Per-Capita\nDuck Production\n Redheads', color = 'Sex') +
  theme_modern(legend.position = 'top',legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)
prod.comp

df1.N_r <- subset(df1.N, species == 'REDHEAD')
df1.N_c <- subset(df1.N, species != 'REDHEAD')

colnames(df1.N_c) <- c('mean.c', 'sd.c','low.c','med.c','upp.c','Rhat.c', 'n.eff.c','species.c','year')
colnames(df1.N_r) <- c('mean.r', 'sd.r','low.r','med.r','upp.r','Rhat.r', 'n.eff.r','species.r','year')

N_ <- left_join(df1.N_c,df1.N_r)

N.comp <-  ggplot(data = N_, aes(x = mean.c/1000, y = mean.r/1000)) +
  geom_point(aes(x = mean.c/1000, y = mean.r/1000), color = '#16a085', size = 2, alpha = .75) +
  geom_smooth(method = 'lm', formula= y~x,color = '#16a085', size = 1.5) +
  geom_errorbar(aes(y = mean.c/1000, ymin=low.r/1000, ymax=upp.r/1000), color = '#16a085', alpha = .5, linetype = 'dashed')+  
  geom_errorbarh(aes(y = mean.r/1000,xmin=low.c/1000, xmax=upp.c/1000), color = '#16a085', alpha = .5, linetype = 'dashed')+  
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1) +
  labs(x = 'Breeding Abundance (in thousands)\n Canvasbacks', y = 'Breeding Abundance (in thousands)\n Redheads\n') +
  theme_modern(legend.text.size = 10, legend.title.size = 11, axis.title.size = 11, axis.text.size = 10)
N.comp 


lN.comp <-  ggplot(data = N_, aes(x = log(mean.c), y = log(mean.r))) +
  geom_point(aes(x = log(mean.c), y = log(mean.r)), color = '#16a085', size = 2) +
  geom_smooth(method = 'lm', formula= y~x,color = '#16a085',) +
  geom_errorbar(aes(y = log(mean.c), ymin= log(low.r), ymax= log(upp.r)), color = '#16a085', alpha = .5, linetype = 'dashed')+  
  geom_errorbarh(aes(y = log(mean.r),xmin=log(low.c), xmax=log(upp.c)), color = '#16a085', alpha = .5, linetype = 'dashed')+  
  scale_color_flat_d() + geom_abline(intercept = 0, slope = 1) +
  labs(x = 'Breeding Abundance (log-scale)\n Canvasbacks', y = 'Breeding Abundance (log-scale)\n Redheads') +
  theme_modern()
lN.comp 




nc <- ggplot(data = df1.N,aes(x = year, y = (mean/1000))) +
  geom_pointrange2(aes(x = year, y = (mean/1000), ymin = (`2.5%`/1000), ymax = (`97.5%`/1000)), size = .65, fill = 'black', shape = 21,linetype = 'dotted', position = position_dodge2(.5))+
  facet_wrap( ~ species,scales = 'free' ) + 
  geom_hline(yintercept = 0, color = NA) +
  geom_hline(yintercept = 2250, color = NA) +
  theme_modern(legend.position = 'none') + labs(y = 'Breeding abundance in PPR (x1000)', x = '') + theme(axis.text.y = element_text(size = 11, angle = 90))


popdyn <- ggarrange(pc,nc, common.legend = FALSE, labels = c('A','B'), nrow = 2)
popdyn

#ggsave(popdyn, file = 'Figures\\popdyn.tiff',height = 10, width = 10, dpi = 320)





popdyn2  <- ggarrange(lam_com,prod.comp, common.legend = FALSE, labels = c('A','B'), nrow = 1)
growth2  <- ggarrange(prod.comp, c.lam_com, lam_com,N.comp, common.legend = FALSE, labels = c('A','B','C','D'), ncol=2, nrow = 2)

#ggsave(growth2, file = 'growth2.tiff', dpi = 320, height = 10, width = 10)


rho <- MCMCsummary(file_list[[2]],'rho')
rho_plot <- matrix(rho$mean, 4, 4)
row.names(rho_plot) <- colnames(rho_plot) <- c('CAN:W','RED:W','CAN:H','RED:H')
corrplot(rho_plot, type = 'upper',  addCoef.col = "black", tl.srt = 35)

########################################################################################################

setwd("C:\\Users\\gibsond\\Desktop\\Ducks\\Data\\bpop\\stratum")

load('..\\..\\..\\model_environment.rdata')
file_list <- list(mcmcList1, mcmcList2)

b.tr <- MCMCsummary(file_list[[1]], 'beta.trend')
b.tr$Lat <- (lat.p[3:26] * 2.985406) + 49.70753
nu.pond<- MCMCsummary(file_list[[1]], 'nu.pond')
mu.pond<- MCMCsummary(file_list[[1]], 'mu.all')
phi.tr <- MCMCsummary(file_list[[2]], 'phi')
b.sst <- MCMCsummary(file_list[[2]], 'beta.sst')


bpop_plot        <- bpop_sub[-c(1:2),]
bpop_plot$Growth <- round((exp(b.tr$`50%`) -1 ) * 100 ,2)
bpop_plot $AR    <- round(phi.tr$`50%`,2)
bpop_plot $inter <- round(nu.pond$`50%`,2)
bpop_coords.test <- bpop_coords.1[-c(1:2), ]
bpop_coords.test$Strata <- bpop_sub2$code[-c(1:2)]


trends <- cbind.data.frame(Strata = bpop_sub$code[-c(1:2)], Area = bpop_sub$area[-c(1:2)], Lat = b.tr$Lat, 
                           Intercept = nu.pond$mean, Inter.Low = nu.pond$`2.5%`,Inter.Upp = nu.pond$`97.5%`,
                           Trend = b.tr$mean, Trend.Low = b.tr$`2.5%`, Trend.Upp = b.tr$`97.5%`,
                           SST = b.sst$mean, SST.Low = b.sst$`2.5%`, SST.Upp = b.sst$`97.5%`,
                           Phi = phi.tr$mean, Phi.Low = phi.tr$`2.5%`, Phi.Upp = phi.tr$`97.5%` )


trends$Area <- as.numeric(trends$Area)
trends <- trends[order(trends$Lat),]

s1 <- ggplot(data = trends) +
  geom_pointrange2(aes(x = reorder(Strata, Lat), y = (exp(Trend) - 1) * 100, ymin = (exp(Trend.Low)-1)*100, ymax = (exp(Trend.Upp)-1)*100),  size = .65)+
  geom_hline(yintercept = 0) +
  labs(y = 'Annual Percent Increase in\n Pond Abundance', x = 'Stratum Number Arranged By Increasing Latitude') +
  theme_modern() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
s1

s2 <- ggplot(data = trends) +
  geom_pointrange2(aes(x = reorder(Strata, Lat),  y = Phi, ymin = Phi.Low, ymax = Phi.Upp), size = .65)+
  geom_hline(yintercept = 0) +
  labs(y = 'Residual Temporal \nAutocorrelation in Pond Density', x = 'Stratum Number Arranged By Increasing Latitude') +
  theme_modern() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
s2

s3 <- ggplot(data = trends) +
  geom_pointrange2(aes(x = reorder(Strata, Lat), y = exp(Intercept), ymin = exp(Inter.Low), ymax = exp(Inter.Upp)), size = .65)+
  labs(y = 'Spatial Variation in\n Mean Pond Density (x 1000)', x = 'Stratum Number (Arranged By Increasing Latitude)') +
  theme_modern() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
s3

s4 <- ggplot(data = trends) +
  geom_pointrange2(aes(x = reorder(Strata, Lat), y = (exp(SST) - 1) * 100, ymin = (exp(SST.Low) - 1) * 100, ymax = (exp(SST.Upp) - 1) * 100),  size = .65)+
  geom_hline(yintercept = 1) +
  labs(y = 'Percent decrease in pond abundance per standard deviation\n increase in Pacific-North American pattern', x = 'Stratum Number Arranged By Increasing Latitude') +
  theme_modern() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
s4

test3 <-  ggarrange( s1,s4,s2,s3, ncol =2, nrow = 2, labels = c('A','B','C','D'))


#ggsave(test3, file = 'pond_trends.tiff' ,height = 12, width = 12, dpi = 320)



#########################################################

bpop_plota        <- bpop_sub2

library("rnaturalearth")
library("rnaturalearthdata")
library("devtools") 
library("sf")
devtools::install_github("ropensci/rnaturalearthhires")
library("rnaturalearth")

world  <- ne_countries(scale = "medium", returnclass = "sf")
world  <- subset(world, admin == 'Mexico' | admin == 'Canada' | admin == 'United States of America')
states <- ne_states(returnclass = "sf")
states <- subset(states, admin == 'Canada'| admin == 'United States of America')

bpop_plot$C.growth <- ifelse(bpop_plot$Growth > 0, 'A','B')

setwd('C://Users//gibsond//Desktop//Ducks')
pp_shp <- st_read('gmannppr//gmannppr.shp')


bplot <- ggplot(data = bpop_plot) +
  geom_sf(data = states) +
  geom_sf(data = bpop_plota, fill = 'red', alpha = .15) +
  geom_sf(data = bpop_plot, aes(fill = exp(inter))) +
  geom_sf(data = pp_shp, fill = NA, color = 'red', size = .75, linetype = 'dashed') +
  # ggrepel::geom_label_repel(data = bpop_plot , aes(label = Growth, geometry = geometry), size = 5,
  #               stat = "sf_coordinates",
  #               min.segment.length = 0) +
  scale_fill_gradient2(midpoint = mean(exp(bpop_plot$inter)),low = 'dodgerblue', mid = 'white', high = 'goldenrod',  breaks=c(0, 100, 200, 300)) +
  coord_sf(crs = st_crs(2163), xlim = c(-1250000, 750000), ylim = c(-500000,  1750000))+
  labs(fill = 'Mean number of \nponds (x 1000)  ', y ='', x = '', color = 'Annual growth in \nregional pond numbers') +
  theme_modern(legend.position = 'top',axis.text.size = 10,legend.title.size = 10,legend.text.size = 10) +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

aplot <- ggplot(data = bpop_plot) +
  geom_sf(data = states) +
  geom_sf(data = bpop_plota, fill = 'red', alpha = .15) +
  geom_sf(data = bpop_plot, aes(fill = Growth)) +
  geom_sf(data = pp_shp, fill = NA, color = 'red', size = .75, linetype = 'dashed') +
  # ggrepel::geom_label_repel(data = bpop_plot , aes(label = Growth, geometry = geometry), size = 5,
  #               stat = "sf_coordinates",
  #               min.segment.length = 0) +
  scale_fill_gradient2(midpoint = mean(bpop_plot$Growth),low = 'dodgerblue', mid = 'white', high = 'goldenrod', breaks=c(0,1,2)) +
  coord_sf(crs = st_crs(2163), xlim = c(-1250000, 750000), ylim = c(-500000,  1750000))+
  labs(fill = 'Percent\nGrowth   ', y ='', x = '', color = 'Annual growth in \nregional pond numbers') +
  theme_modern(legend.position = 'top',axis.text.size = 10,legend.title.size = 10,legend.text.size = 10) +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
aplot


all.ponds <- MCMCsummary(file_list[[1]], 'mu.all')
all.ponds$Year <- 1961:2019
# Average (Nov-Apr)  Pacific North American Index 1961 - 2020
pna <- c(0.393333333,-0.605,-0.078333333,0.121666667,-0.778333333,-0.431666667,-1,-0.243333333,-0.348333333,0.658333333,-0.461666667,-1.028333333,0.191666667,
         -0.685,-0.051666667,-0.348333333,0.938333333,0.221666667,-0.541666667,0.538333333,0.743333333,-0.615,1.065,1.063333333,-0.341666667,0.318333333,0.881666667,
         0.988333333,-0.476666667,-0.35,-0.27,0.5,0.388333333,-0.07,0.171666667,0.025,-0.206666667,0.971666667,0.28,0.32,0.591666667,-0.295,0.87,0.236666667,0.481666667,
         0.186666667,0.421666667,-0.046666667,-0.241666667,0.99,-0.635,0.196666667,-0.438333333,-0.33,0.361666667,0.948333333,0.383333333,-0.6,0.081666667,-0.55)
n <- 25000
climate_ponds <- climate <- array(NA, dim = c(24,59,n))
mu.ponts <- matrix(NA, 24,n)
obs_ponds<- matrix(NA, nrow = 59, ncol = n)

beta.sst <- MCMCsummary(file_list[[1]], params = 'beta.sst')
beta.area <- MCMCsummary(file_list[[2]], params = 'beta.area')
beta.area.trend <- MCMCsummary(file_list[[2]], params = 'beta.trend.area')

for(i in 1:24){
  mu.ponts[i,] <- rnorm(n, nu.pond$mean[i],nu.pond$sd[i]) 
  for(j in 1:59){
    climate[i,j,] <- rnorm(n, beta.sst$mean[i],beta.sst$sd[i]) *  pna[j] #- exp(nu.pond$mean[i])
    for(k in 1:n){
      climate_ponds[i,j,k] <- exp(mu.ponts[i,k] + climate[i,j,k]) - exp(mu.ponts[i,k])
    }
  }
}

all.climate_ponds <- apply(climate_ponds,  3,colSums)

all.clim.ponds    <- cbind.data.frame(Year = 1961:2019, mean = apply(all.climate_ponds,1, mean), sd= apply(all.climate_ponds,1, sd))


beta_ponds <- MCMCsummary(file_list[[2]], params = c('mu.sst', 'beta.area','beta.trend.area'))

beta_ponds $Type <- c('PNA','\nArea\n','Area\nx\nTrend')

cn <-ggplot(data = beta_ponds, aes(x = (Type), y = mean)) +
  geom_pointrange2(aes(x = Type, y = mean,ymax = `97.5%`, ymin = `2.5%`), size = .75)+
  labs(y = '\nSlope Coefficients\n', x = '')+ geom_hline(yintercept = 0)+
  theme_modern() +theme(axis.text.y = element_text(size = 11, angle = 90))

pn <-  ggplot(data = all.ponds, aes(x = Year, y = mean)) +
  geom_pointrange2(aes(x = Year, y = mean,ymax = `97.5%`, ymin = `2.5%`), size = .75)+
  # geom_pointrange2(data = all_obs, aes(x = Year, y = Ponds, ymax = Ponds + SE*1.96, ymin = Ponds - SE*1.96), color = 'red', size = .75)+
  labs(y = 'Annual number of ponds within \nPrarie Potholes (x1000)', x = '')+
  theme_modern() + theme(axis.text.y = element_text(size = 11, angle = 90))

cp <-  ggplot(data = all.clim.ponds, aes(x = Year, y = mean)) +
  geom_pointrange2( aes(x = Year, y = mean, ymin = mean - sd*1.96, ymax = mean + sd*1.96), color = 'black', size = .75)+
  geom_smooth(aes(x = Year, y = mean), method = 'gam', se = FALSE) +
  geom_hline(yintercept = 0) +
  labs(y = 'Variation in pond numbers \nassociated with shifts in\nPacific-North American index (x1000)', x = '')+
  theme_modern() + theme(axis.text.y = element_text(size = 11, angle = 90))

sst.df <- cbind.data.frame(Year= 1961:2019, pna[-60])
st<- ggplot(data = sst.df, aes(x = Year, y = pna)) +
  geom_line( aes(x = Year, y = sst), color = 'black')+
  labs(y = 'Variation in de-trended \nsea surface temperatures\n', x = '')+
  theme_modern()



test1 <- ggarrange(bplot, aplot, labels = c('A','B'), nrow = 2, ncol = 1)
test2 <- ggarrange(pn, labels = 'C', ggarrange( cn,cp, ncol =2, nrow = 1, labels = c('D','E'), widths = c(1,1.25)), ncol = 1, nrow = 2)

test <- ggarrange(test1, test2, ncol = 2, widths = c(.75,1.25))

#ggsave(test, file = 'pond_numbers.tiff' ,height = 10, width = 12, dpi = 320)


setwd('C:\\Users\\gibsond\\Desktop\\Ducks')
df1 <- read.csv('CHB_SAL.csv')
df1$SampleDepth <- round_any( df1$SampleDepth, 0.5)
df1$lat_round <- round_any(df1$Latitude,.08333)
df1$long_round <- round_any(df1$Longitude,.08333)
df1$Observation <- paste(df1$lat_round,"-",df1$long_round,'-',df1$SampleDate,'-', df1$SampleDepth)
df1$Unique1 <-  paste(df1$lat_round,"-",df1$long_round)
df1 <- df1%>% group_by(Observation) %>% sample_n(size = 1)

df1 <- subset(df1, Year >= 1961 & Year < 2020 & Longitude < -76.15 &  Longitude > -76.5)
df1 <- df1[complete.cases(df1$ReportedValue),]
df1 <- subset(df1, SampleDepth <= 1 )
df1 <- subset(df1,ReportedValue > .1 & ReportedValue<30 )

df_unique1 <- df1[!duplicated(df1$Unique1),]
study_points1 <- cbind.data.frame(long = df_unique1$long_round, lat = df_unique1$lat_round)
blockcoord1 <- cbind(study_points1$long,study_points1$lat)
neigh1      <- dnearneigh(blockcoord1, d1 = 0, d2 = 20, longlat = TRUE)
winnb1 <- nb2WB(neigh1)          

distance <- as.matrix(dist(study_points1))

covs1 <- df1[,c('ACE','AMO','NAO','PREC','SampleDepth')]
covs1$SampleDepth <- round_any(covs1$SampleDepth,0.5)
covs_scale1 <- scale(covs1)
covs_scale1[is.na(covs_scale1)] <- 0

covs_scale1[,5] <- as.factor(covs_scale1[,5])

year1       <- df1$Year - 1960
month1      <- df1$Month
system.ind  <- ifelse(year1 < 24,2,3)

df <- read.csv('texas_salinity.csv')
df <- subset(df, Depth == 0.3)
df <- subset(df, !grepl('RIVER', SiteName))
df <- subset(df, !grepl('CREEK', SiteName))

df$lat_round <- round_any(df$Latitude,.08333)
df$long_round <- round_any(df$Longitude,.08333)
df$Unique1 <-  paste(df$lat_round,"-",df$long_round)
df$Observation <- paste(df$lat_round,"-",df$long_round,'-',df$SiteName,'-',df$Year,'-',df$Month)
df <- df%>% group_by(Observation) %>% sample_n(size = 1)

df_unique <- df[!duplicated(df$Unique1),]
study_points <- cbind.data.frame(long = df_unique$long_round, lat = df_unique$lat_round)
blockcoord <- cbind(study_points$long,study_points$lat)
neigh      <- dnearneigh(blockcoord, d1 = 0, d2 = 20, longlat = TRUE)
winnb <- nb2WB(neigh)          

covs <- df[,c('ACE','AMO6','NAO6','PREC6','FLOW6')]
covs_scale <- scale(covs)
covs_scale[is.na(covs_scale )] <- 0

year  <- df$Year - 1960
month <- df$Month

preds <- read.csv("fall_predictions.csv")
preds <- subset(preds, ï..Year > 1960)
preds_lgm <- cbind.data.frame(Ace = preds$Ace, AMO = preds$AMO, NAO = preds$NAO, PREC = preds$TEX.PREC, FLOW = preds$TEX.FLOW)
preds_chb <- cbind.data.frame(Ace = preds$Ace, AMO = preds$AMO, NAO = preds$NAO, PREC = preds$CHB.PREC, FLOW = rep(NA))

preds_lgm_sc <-  preds_chb_sc <- matrix(NA, 59,5)

for(i in 1:dim(preds_lgm)[1]){
  for(j in 1:dim(preds_lgm)[2]){
    preds_lgm_sc[i,j] <- (preds_lgm[i,j] - attr(covs_scale,"scaled:center")[j])/ attr(covs_scale,"scaled:scale")[j]
  }
  for(j in 1:(dim(preds_chb)[2] - 1)){
    preds_chb_sc[i,j] <- (preds_chb[i,j] - attr(covs_scale1,"scaled:center")[j])/ attr(covs_scale1,"scaled:scale")[j]
  }
}

preds_chb_sc[is.na(preds_chb_sc)] <- 0
preds_lgm_sc[is.na(preds_lgm_sc)] <- 0

library(abind)
pred_covs <- abind(preds_lgm_sc,preds_chb_sc, along = 3)

sal.std <- MCMCsummary(file_list[[1]], 'eps.sal')

salinity <- cbind.data.frame(expand.grid(Year = 1961:2019,Region = c('Chesapeake Bay','Laguna Madre')), sal.std)
salinity$mu <- rep(c(mean(car_chb$mean),mean(car_lgm$mean)), each = 59)

salinity$Line <- ifelse(salinity$Region == 'Chesapeake Bay' & salinity$Year == '1965' | salinity$Year == '1983', 1, 
                        ifelse(salinity$Region != 'Chesapeake Bay' & salinity$Year > 2015,1,2))

sp <-  ggplot(data =salinity, aes(x = Year, y = exp(mu + mean) )) +
  geom_pointrange2(data = subset(salinity, Region == 'Chesapeake Bay'| Year > 1968 & mean > -5),
                   aes(x = Year, y =exp(mu + mean), ymin = exp(mu + `2.5%`), ymax = exp(mu + `97.5%`), linetype = as.factor(Line),  fill = Region), shape = 21, size = .75)+
  geom_ribbon(data = subset(salinity, Region != 'Chesapeake Bay'& Year < 1969 ), aes(x = Year, y =exp(mu + mean),  ymin =exp(mu + `2.5%`),  ymax = exp(mu + `97.5%`),  fill = Region), alpha = .25)+
  geom_line(data = subset(salinity, Region != 'Chesapeake Bay'& Year < 1969), aes(x = Year, y =exp(mu + mean)), color = 'black')+
  scale_linetype_manual(values = c('dashed','solid')) +
  labs(y = 'Annual variation in total dissolved solids (mg/L) during the fall', x = '')+
  scale_fill_flat_d() + facet_wrap(~Region, scales = 'free_y', ncol = 1) + #coord_cartesian(ylim = c(0, 80))+
  theme_modern(legend.position = 'none')
sp 


car_lgm <- MCMCsummary(file_list[[1]], 'car_lgm')
car_lgm <- cbind.data.frame(car_lgm, study_points)
car_chb <- MCMCsummary(file_list[[1]], 'car_chb')
car_chb <- cbind.data.frame(car_chb, study_points1)


data(coastlineWorldFine, package = "ocedata")
coast <- as.data.frame(coastlineWorldFine@data)

CHB_fig <-
  ggplot() + 
  geom_polygon(data = coast, aes(x = longitude, y = latitude), alpha = .5) +  coord_sf(xlim = c(-78, -75), ylim = c(37, 40), expand = FALSE) +
  geom_tile(data = car_chb , aes(x = long, y = lat, fill = exp(mean))) + scale_fill_metro_c() +
  labs(y ='',x='',fill = 'Average Salinity (mg/L)\n in Chesapeake Bay  ') +
  theme_modern(legend.position = 'top',axis.text.size = 10,legend.title.size = 10,legend.text.size = 10) +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
CHB_fig 

LGM_fig <- ggplot() +
  geom_polygon(data = coast, aes(x = longitude, y = latitude), alpha = .5) +
  geom_tile(data = car_lgm , aes(x = long, y = lat, fill = exp(mean))) + scale_fill_metro_c() +
  coord_sf(xlim = c(-98, -96.5), ylim = c(26.5, 28), expand = FALSE)+
  labs(y ='',x='',fill = 'Average Salinity (mg/L)\nin Laguna Madre  ') +
  theme_modern(legend.position = 'top',axis.text.size = 10,legend.title.size = 10,legend.text.size = 10) +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
LGM_fig


covars <- rbind.data.frame(MCMCsummary(file_list[[2]], 'beta.sal'),
                           MCMCsummary(file_list[[2]], 'beta.flow'),
                           MCMCsummary(file_list[[2]], 'beta.inter')
                           
)

covars$Region <- c(rep(c('Laguna Madre','Chesapeake Bay'), times = 4),'Laguna Madre','Laguna Madre','Chesapeake Bay')

covars$Parameter <- c(rep(c('ACE','AMO','NAO','PREC'), each = 2), 'FLOW','ACExPREC','ACExPREC')

sp1 <-  ggplot(data = covars, aes(x = Parameter, y = exp(mean))) +
  geom_pointrange2( aes(x = Parameter, y =exp(mean), ymin =exp(`2.5%`), ymax = exp(`97.5%`), fill = Region), shape = 21, size = .75, position = position_dodge2(.5))+
  labs(y = 'Odds Ratios: Climate/weather effects on \nregional salinity', x = '')+
  geom_hline(yintercept = 1) +
  scale_fill_flat_d()+
  theme_modern(legend.position = 'top')
sp1


MCMCsummary(file_list[[2]],'rho.sal')

test1 <-  ggarrange(sp, labels = 'A', ggarrange(sp1, labels = 'B', ggarrange(CHB_fig, LGM_fig, labels = c('C','D'), nrow = 1, ncol = 2), ncol = 1, nrow = 2), nrow = 1, ncol = 2, widths = c(.75,1))

#ggsave(test1, file = 'salinity.tiff' ,height = 10, width = 14, dpi = 320)




sig <- MCMCsummary(file_list[[2]],'sig')
tau <- MCMCsummary(file_list[[2]],'tau')
decay <- MCMCsummary(file_list[[2]],'decay.pond')

dists <- seq(from = 0, to = 15, length.out = 100)
n <- 10000
sigma <- matrix(NA,100,n)
sigma[1,] <- rnorm(n, sig$mean, sig$sd) + rnorm(n, tau$mean, tau$sd)
for(i in 2:100){
  sigma[i,] <- rnorm(n, sig$mean, sig$sd) * exp(-rnorm(n, decay$mean, decay$sd) * dists[i])
}

sigma_mean <- apply( sigma,1,mean)
sigma_sd   <- apply( sigma,1,sd)


spatial_error <- cbind.data.frame(Distance = dists * 1000, mu = sigma_mean, sd = sigma_sd)

sp3 <-  ggplot(data =spatial_error , aes(x = Distance, y = sqrt(mu))) +
  geom_line( aes(x = Distance, y = mu), size = .75)+
  geom_segment(x = 0, xend =  log(.5)/-decay$mean * 1000, y =  0.1042339, yend = 0.1042339, linetype = 'dashed') +
  geom_segment(x = log(.5)/-decay$mean * 1000, xend =  log(.5)/-decay$mean * 1000, y = 0, yend = 0.1042339, linetype = 'dashed') +
  geom_ribbon( aes(x = Distance, y = mu, ymin = mu - sd*1.96, ymax = mu + sd*1.96), size = .75, alpha = .25)+
  labs(y = 'Residual similarity in pond numbers (log-scale) among strata', x = 'Distance between strata centroids (km)')+
  theme_modern(axis.text.size = 10,legend.title.size = 10,legend.text.size = 10)
sp3

#ggsave(sp3, file = 'decay_curv.tiff', height = 6.5, width = 6.5, dpi = 320 )
############################################################

mean_sur   <- MCMCsummary(mcmcList2, 'mu.sur')
mean_kappa <- MCMCsummary(mcmcList2, 'mu.kappa')

rho_sur    <- MCMCpstr(mcmcList2, params = 'rho', type = 'chains')

rhos <- cbind.data.frame('CAN-RED\nWinter' = rho_sur$rho[1, 2,], "CAN\nWinter-Harvest" = rho_sur$rho[1, 3,],
                         'RED\nWinter-Harvest' = rho_sur$rho[2, 4,], "CAN-RED\nHarvest" =rho_sur$rho[3, 4,],
                         "CAN-RED\nSummer" = rho_sur$rho[5, 6,])
rhos_m <- melt(rhos)

mean_sur   <- cbind.data.frame(mean_sur, expand.grid(Species = c('Canvasback','Redhead'), Sex = c('Female','Male'), Type = c('HY\nWinter', 'AHY\nWinter','HY & AHY\nSummer')))
mean_kappa <- cbind.data.frame(mean_kappa, expand.grid(Type = c('HY\nHarvest', 'AHY\nHarvest'), Species = c('Canvasback','Redhead'), Sex = c('Female','Male')))

means <- rbind.data.frame(mean_sur, mean_kappa)

means$Risk     <- exp(means$mean)
means$RiskLow  <- exp(means$`2.5%`)
means$RiskHigh <- exp(means$`97.5%`)

colors <- flat_colors(4,12)
m1 <- ggplot() + 
  geom_pointrange2(data = means, aes(x = Type, y = Risk, ymin = RiskLow, ymax = RiskHigh, fill = Sex), size = .75, shape = 21, position = position_dodge2(.5)) +
  facet_wrap(~Species) +
  theme_modern(legend.position = 'top', axis.title.size = 11,  axis.text.size = 10,legend.title.size = 10,legend.text.size = 10) + scale_fill_manual(values =c("#8e44ad","#f39c12")) + labs(y = 'Average risk of mortality', x = 'Cause')

rhos_m$Type <- rep(c('Between','Within','Within','Between','Between'), each = 562500 /5)

library(Ipaper)
m2 <- ggplot() +  geom_boxplot2(data = rhos_m, aes(x = variable, y = value, fill = Type), width.errorbar = .25, width = .5) +
  theme_modern(axis.title.size = 11, axis.text.size = 10,legend.title.size = 10,legend.text.size = 10,legend.position = 'top') + geom_hline(yintercept = 0) + 
  scale_fill_metro_d() + coord_cartesian(ylim = c(-1,1))+
  labs(y = 'Posterior Distribution\nResidual Annual Correlation', x = 'Correlation Type')



fig_means <- ggarrange(m1, m2, labels = c('A','B'), ncol = 1)
#ggsave(fig_means, file = 'fig_means.tiff', dpi = 320, height = 8, width = 8)




SR <- MCMCsummary(mcmcList1, 'sex.ratio')

SR1 <- cbind.data.frame(SR, expand.grid(Species = c('Canvasback','Redhead'),Year = 1961:2019))


SR_1 <- ggplot(data = SR1, aes(x = Year, y = mean)) + 
  geom_pointrange2(data = SR1, aes(x = Year, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = Species), size = .5, shape = 21, position = position_dodge2(1)) +
  geom_smooth(aes(color = Species)) +
  theme_modern(axis.title.size = 11, axis.text.size = 10,legend.title.size = 10,legend.text.size = 10,legend.position = 'top') + scale_color_flat_d()+scale_fill_flat_d()+ geom_hline(yintercept = 1) +
  labs(y = 'Harvest-Corrected Fall Sex Ratio (Females/Male)', x = '')
SR_1
#ggsave(SR_1, file = 'SR_1.tiff', dpi = 320, height = 8, width = 8)

