# We load the R packages that are needed. You have to install them first if you do not have them.
library(nlme)
library(rstan)
library(JMbayes)
library(MASS)
library(nlme)
library(rstan)
library(JMbayes)
library(splines)

# we load the R script extractFrames_new. This R script creates each time 
# the matrices required for the fitted model

source("U://R codes modified//code for the new proposed model//extractFrames_new.R")

# we define the parameters required 
# N: the total number of subjects 
# K: the total number of measurements that each subject has 
# t.max: the maximum time 
# beta1,beta2,beta3,beta4: real coefficients of the variables 
# sigma2.v: the 2 options of the variance of the error terms 
# var.v: the 2 options of variance of the random effects 
# corr.v: the 2 options of the correlation of the random effects 
# optionHC.v: the different parametrizations of the hierarchical models
# scaling.v: the different variable transformations 

N = 200; K = 6; t.max = 15;
beta1 = 15; beta2 = 1.8; beta3 = 2.1; beta4 = 0.5; 
sigma2.v = c(0.5, 5); var.v = c(0.5, 5); corr.v = c(0.1, 0.5)
optionHC.v = c("HC", "HNC")
scaling.v = c("standardized", "centred", "Non")

# we create a data frame with 4800 rows and 47 columns 
summary.table = data.frame(matrix(vector(), 4800, 47))
colnames(summary.table) = c("simulation","sigma2", "var", "corr", "optionHC", "scaling","time", 
                            "mean_beta1","se_mean_beta1","sd_beta1", "n_eff_beta1", "Rhat_beta1",
                            "mean_beta2","se_mean_beta2","sd_beta2", "n_eff_beta2", "Rhat_beta2",
                            "mean_beta3","se_mean_beta3","sd_beta3", "n_eff_beta3", "Rhat_beta3",
                            "mean_beta4","se_mean_beta4","sd_beta4", "n_eff_beta4", "Rhat_beta4",
                            "mean_sigma2","se_mean_sigma2","sd_sigma2", "n_eff_sigma2", "Rhat_sigma2",
                            "mean_D11","se_mean_D11","sd_D11", "n_eff_D11", "Rhat_D11",
                            "mean_D12","se_mean_D12","sd_D12", "n_eff_D12", "Rhat_D12",
                            "mean_D22","se_mean_D22","sd_D22", "n_eff_D22", "Rhat_D22")


# we repeat each of the values of the vector 1:100 48 times 
# and we save it in the column simulation of the data frame summary.table
summary.table$simulation = c(rep(1:100,each=48))

# we repeat 24 times the value 0.5 and 5 and later on
# we repeat the resulting vector 100 times 

summary.table$sigma2 = rep(c(rep(0.5, 24), rep(5, 24)), 100) 

# we repeat 12 times the value 0.5 and 5 and later on
# we repeat the resulting vector 2 times. The final vector
# is repeated 100 times 

summary.table$var = rep(c(rep(c(rep(0.5, 12), rep(5, 12)), 2)), 100) 

# we repeat 6 times the value 0.1 and 0.5 and later on
# we repeat the resulting vector 2 times. The resulting vector 
# is repeated 2 times and the final vector
# is repeated 100 times 

summary.table$corr = rep(c(rep(c(rep(c(rep(0.1, 6), rep(0.5, 6)), 2)), 2)), 100) 

# we repeat 3 times the value HC and HNC and later on
# we repeat the resulting vector 2 times. The resulting vector 
# is repeated 2 times and the final vector
# is repeated 100 times 

summary.table$optionHC = rep(c(rep(c(rep(c(rep("HC", 3), rep("HNC", 3)), 2)), 2)), 100) 

# we repeat the vector c("standardized","centred","Non") 2 times and later 
#on we repeat the results 2 times 
# The resulting vector is repeated 2 times and later on 100 times 
# we repeat the resulting vector 2 times. The resulting vector 

summary.table$scaling = rep(c(rep(c(rep(c(rep(c("standardize","centred","raw"), 2)), 2)), 2)), 100) 



# we start the for loop. We will create for each combination
# of the parameters 100 data sets 

for (d in 1:100){
  for (j in 1:length(sigma2.v)){
    for (k in 1:length(var.v)){
      for (m in 1:length(corr.v)){
        for (n_res in 1:length(optionHC.v)){
          for (l in 1:length(scaling.v)){
            
            ########################
            # Simulation Scenarios #
            ########################
            # we set as sigma2 the j-th element of the sigma2.v
            sigma2 = sigma2.v [j]
            # we set as var the k-th element of the var.v
            var = var.v [k]
            # we set as corr the m-th element of the corr.v
            corr = corr.v [m]
            # we set as optionHC  the n_res element of the oprionHC.v
            optionHC = optionHC.v [n_res]
            # we set as scaling the l-th element of the scaling.v
            scaling = scaling.v [l]
            
            
            ####################
            # Simulate dataset #
            ####################
            
            set.seed(1234 + d)
            
            N = N  # number of subjects 
            K = K  # number of planned repeated measurements per subject 
            t.max = t.max # maximum follow-up time
            
            # Parameters for the mixed effects model 
            betas = c("(Intercept)" = beta1, "year" = beta2, "sex" = beta3, "drug" = beta4)
            sigma.y = sqrt(sigma2)  # measurement error standard deviation
            # covariance of the random effects  
            cov = corr * var
            # variance-covariance matrix of the random effects
            D = matrix(c(var,  cov,  
                         cov,  var), nrow = 2, ncol = 2) 
            
            year = c(replicate(N, c(0, sort(runif(K-1, 0, t.max)))))  # we set elements of the year  
            sex = sample(c(rep(0, N/2), rep(1, N/2))) # group indicator, i.e., '0' male, '1' female
            drug<-rep(sample(c(0,1),N,replace=TRUE),each=K) # we set the variable drug to take values 0/1
            
            # we define the data frame DF. This data frame contains 
            # the variables sex, year and drug
            DF = data.frame(sex = factor(rep(sex, each = K)), year, drug)  
            # we define the design matrix X of the fixed part: year + sex + drug
            X = model.matrix(~ year + sex + drug, data = DF)
            # we define the design matrix Z of the random part: year
            Z = model.matrix(~ year, data = DF)
            # we take N random elements from a multivariate normal distribution
            # with mean zero and variance covariance matrix D 
            # these are the random terms 
            b = mvrnorm(N, rep(0, 2), Sigma = D)
            b = as.matrix(b)
            # we repeat each id of the patients K times each 
            id = rep(1:N, each = K)
            # we calculate the expected value of the model
            
            eta.y = as.vector(X %*% betas + rowSums(Z * b[id, ]))
           # we take N*K random values from a normal distribution 
            # with mean eta.y and variance sigma.y 
             y = rnorm(N * K, eta.y, sigma.y)
            
            # Final data set
            data = data.frame(y = y, year = DF$year, sex = DF$sex, drug = DF$drug, id = id)
            
            #############
            # fit model #
            #############
            
            formulas = list(y ~  year + sex + drug + ( year | id))
            data = data
            families = list(gaussian)
            engine =  "STAN"
            overdispersion = FALSE
            priors = NULL
            init = NULL
            control = NULL
            
            # we paste the path of the stan codes and we save it with 
            # the bame file.stan
            file.stan = paste0("V://Users//055609(A_Lika)//HC+NHC_models_simulations//stan codes new proposed 20-10-2020//stan_", optionHC, "_", scaling,".stan")
           # we name as sim.name the result of the paste command
            # this will be used for the definition of the 
            # columns of the data frame summary.table  
             sim.name = paste0("sim_sigma", sigma2, "_var", var, "_corr", corr,"_",optionHC,"_",scaling)
            
    ## We fit the HC model with standardized variables 
             #if the optionHC="HC" and scaling="standardize"
            
            if(((optionHC=="HC")&(scaling=="standardized"))){
              
              # we create the matrices for the model 
              #y ~  year + sex + drug + ( year | id) with the extractFrames_new 
              # function 
             
               necessary<-extractFrames_new(y ~  year + sex + drug + ( year | id),data=data)
              
              rstan_options(auto_write = TRUE)
              n<-necessary$n # total nymber of subjects 
              n_RE<-2 # total number of the random terms 
              N1<-necessary$N # total nymber of observations 
              ncx1<-necessary$ncx # total nymber of variables in the fixed part 
              id1<-necessary$id # the ids of subjects
              RE_ind1<-1:2 # vector from 1 to the number of the error terms 
              # which here is 2 
              y1<-necessary$y # the dependent variable 
              Z_1<-necessary$Z_ # the design matrix of the random part when we have raw data 
              Zs1<-necessary$Zs # the design matrix of the random part when we have standardized data 
              Zv1<-necessary$Zv # the transpose design matrix of the random part when we have raw data 
              Xhc1<-necessary$Xhc # the design matrix of the fixed part when we have raw data 
              ZrowsStart1<-necessary$ZrowsStart # defines where is the first observation of each patient
              ZrowsEnd1<-necessary$ZrowsEnd # defines where is the last observation of each patient
              
              
              Ztinv1<-necessary$Ztinv #  the transpose inverse design matrix of the random part when we have raw data
              Zinv1<-necessary$Zinv # the inverse design matrix of the random part when we have raw data
              XhcS1<-necessary$XhcS # the design matrix of the fixed part when we have standardized data 
              
              SDs_Xhc1<-necessary$SDs_Xhc # the sd of the variables in the fixed part 
              mean_sd_Xhc1<-necessary$mean_sd_Xhc # the mean/sd of the variables in the fixed part 
              
              
              SDs_Z1<-necessary$SDs_Z # the sd of the variables in the random part 
              mean_sd_Z1<-necessary$mean_sd_Z # the mean/sd of the variables in the random part 
              
              # we create the list of the elements required for the stan code 
              # of the HC model with standardized variables 
              pulpdat <- list(n=n,n_RE=n_RE,N1=N1,ncx1=ncx1,id1=id1,RE_ind1=RE_ind1,y1=y1,Zs1=Zs1,
                              scale_betas1=100,scale_sigmas=10,scale_diag_D=10,lkj_shape=2.0, Ztinv1= Ztinv1, Zinv1= Zinv1,Zv1=Zv1,Z_1=Z_1,Xhc1=Xhc1,
                              ZrowsStart1=ZrowsStart1,ZrowsEnd1=ZrowsEnd1,XhcS1=XhcS1,SDs_Xhc1=SDs_Xhc1,mean_sd_Xhc1=mean_sd_Xhc1,SDs_Z1=SDs_Z1,mean_sd_Z1=mean_sd_Z1)
              
              
              
              ## Fitting the HC model with standardized variables in STAN 
              # we run 2 chains with 10000 iterations (the first 4000 iterations are the 
              #warm ups), max_treedepth=20
              # and adapt_delta=0.95
              
              fit <- stan(file=file.stan,data=pulpdat,chains = 2, warmup = 4000, iter = 10000, thin = 1,
                          control = list(max_treedepth = 20, adapt_delta = 0.95))
             
              # we require the summary table of the parameters of the model 
    
              summary1 = as.data.frame(summary(fit, pars = c("betas1","sigma1","D","b"))$summary)
              # we save the summary table 
              save(summary1, file = paste0("V://Users//055609(A_Lika)//HC+NHC_models_simulations//summary_tables_10000//summary_", sim.name, "_d", d, ".RData"))
            }
            
            
             ## We fit the HC model with centred variables 
             #if the optionHC="HC" and scaling="centred"
             
             
            if(((optionHC=="HC")&(scaling=="centred"))){
              
              # we create the matrices for the model 
              #y ~  year + sex + drug + ( year | id) with the extractFrames_new 
              # function 
              
              necessary<-extractFrames_new(y ~  year + sex + drug + ( year | id),data=data)
              
              rstan_options(auto_write = TRUE)
              
              
              
              n<-necessary$n # total number of subjects 
              n_RE<-2 # total number of the random terms        
              N1<-necessary$N # total number of observations
              ncx1<-necessary$ncx # number of variables in the matrix of the fixed part 
              id1<-necessary$id # the vector with the ids of the subjects
              RE_ind1<-1:2 # vector from 1 to the number of the error terms 
              # which here is 2
  
              y1<-necessary$y # the dependent variable
              Z_1<-necessary$Z_ # the design matrix of the random part when we have raw variables
              Zc1<-necessary$Zc # the design matrix of the random part when we have centred variables
              Zv1<-necessary$Zv # the transpose design matrix of the random part when we have raw variables
              Xhc1<-necessary$Xhc # the design matrix of the fixed part when we have raw variables
              ZrowsStart1<-necessary$ZrowsStart# defines where is the first observation of each patient
              ZrowsEnd1<-necessary$ZrowsEnd# defines where is the last observation of each patient
              
              
              Ztinv1<-necessary$Ztinv # the transpose inverse design matrix of the random part when we have raw variables
              Zinv1<-necessary$Zinv # # the inverse design matrix of the random part when we have raw variables
              XhcC1<-necessary$XhcC # the design matrix of the random part when we have centred variables
              
              means_Xhc1<-necessary$means_Xhc # vector with the means of the variables in the fixed part 
              means_Z1<-necessary$means_Z # vector with the means of the variables in the random part 
              
              
              # we create the list of the elements required for the stan code 
              # of the HC model with centred variables
              
              
              pulpdat <- list(n=n,n_RE=n_RE,N1=N1,ncx1=ncx1,id1=id1,RE_ind1=RE_ind1,y1=y1,Zc1=Zc1,
                              scale_betas1=100,scale_sigmas=10,scale_diag_D=10,lkj_shape=2.0, Ztinv1= Ztinv1, Zinv1= Zinv1,Zv1=Zv1,Z_1=Z_1,Xhc1=Xhc1,
                              ZrowsStart1=ZrowsStart1,ZrowsEnd1=ZrowsEnd1,XhcC1=XhcC1,means_Z1=means_Z1,means_Xhc1=means_Xhc1)
              
              
              
             
              ## Fitting the HC model with centred variables in STAN 
              # we run 2 chains with 10000 iterations (the first 4000 iterations are the 
              #warm ups), max_treedepth=20
              # and adapt_delta=0.95
              
              fit <- stan(file=file.stan,data=pulpdat,chains = 2, warmup = 4000, iter = 10000, thin = 1,
                          control = list(max_treedepth = 20, adapt_delta = 0.95))
              # we require the summary table of the parameters of the model
              summary1 = as.data.frame(summary(fit, pars = c("betas1","sigma1","D","b"))$summary)
              # we save the summary1 
              save(summary1, file = paste0("V://Users//055609(A_Lika)//HC+NHC_models_simulations//summary_tables_10000//summary_", sim.name, "_d", d, ".RData"))
            }
            
             ## We fit the HC model with raw variables 
             #if the optionHC="HC" and scaling="raw"
             # and We follow the same process as was described previously for the other models 
            
            if(((optionHC=="HC")&(scaling=="raw"))){
              
              
              necessary<-extractFrames_new(y ~  year + sex + drug + ( year | id),data=data)
              
              rstan_options(auto_write = TRUE)
              
              n<-necessary$n # total number of subjects 
              n_RE<-2 # total number of random effects 
              N1<-necessary$N # total number of observations 
              ncx1<-necessary$ncx # number of variables in the matrix of the fixed part 
              id1<-necessary$id  # the vector with the ids of the subjects  
              RE_ind1<-1:2 # vector from 1 to the total number of the random terms 
              y1<-necessary$y # the dependent variable 
              Z_1<-necessary$Z_ # the design matrix of the random part when we have raw variables 
              
              Xhc1<-necessary$Xhc # the design matrix of the fixed part when we have raw variables 
              
              
              # list of data for the stan code 
              
              pulpdat <- list(n=n,n_RE=n_RE,N1=N1,ncx1=ncx1,id1=id1,RE_ind1=RE_ind1,y1=y1,Z_1=Z_1,
                              scale_betas1=100,scale_sigmas=10,scale_diag_D=10,lkj_shape=2.0,Xhc1=Xhc1)
              
              
              ## Fitting the model
              
              fit <- stan(file=file.stan,data=pulpdat,chains = 2, warmup = 4000, iter = 10000, thin = 1,
                          control = list(max_treedepth = 20, adapt_delta = 0.95))
              
              summary1 = as.data.frame(summary(fit, pars = c("betas1","sigma1","D","b"))$summary)
              save(summary1, file = paste0("V://Users//055609(A_Lika)//HC+NHC_models_simulations//summary_tables_10000//summary_", sim.name, "_d", d, ".RData"))
            }
            
            
            
            
             ## We fit the HNC model with standardized variables 
             #if the optionHC="HNC" and scaling="standardized"
             
            
            
            if(((optionHC=="HNC")&(scaling=="standardized"))){
              
              necessary<-extractFrames_new(y ~  year + sex + drug + ( year | id),data=data)
              
              rstan_options(auto_write = TRUE)
              
              n<-necessary$n # total number of subjects
              n_RE<-2 # total number of random terms
              N1<-necessary$N # total number of observations 
              ncx1<-necessary$ncx # total number of variables in the design matrix of the fixed part
              id1<-necessary$id # the ids of the subjects 
              RE_ind1<-1:2 # vector from 1 to the total number of the random terms 
              y1<-necessary$y # dependent variable
              Z_1<-necessary$Z_ # the design matrix of the random part when we have raw data
              Zs1<-necessary$Zs # the design matrix of the random part when we have standardized data
              Zv1<-necessary$Zv # the transposed design matrix of the random part when we have raw data
              Xhc1<-necessary$Xhc # # the design matrix of the fixed part when we have raw data
              ZrowsStart1<-necessary$ZrowsStart # defines where is the first observation of each patient
              ZrowsEnd1<-necessary$ZrowsEnd # defined where is the last observation of each patient 
              
              
              Ztinv1<-necessary$Ztinv # the transposed inverse design matrix of the random part when we have raw data
              Zinv1<-necessary$Zinv # inverse design matrix of the random part when we have raw data
              Xs1<-necessary$Xs # design matrix of the fixed part when we have standardized data
              
              SDs_X1<-necessary$SDs_X # vector with the sds of the variables in the fixed part 
              mean_sd_X1<-necessary$mean_sd_X # vector with the mean/sd of the variables in the fixed part 
              
              
              
              SDs_Z1<-necessary$SDs_Z # vector with the sds of the variables in the random part
              mean_sd_Z1<-necessary$mean_sd_Z # vector with the mean/sd of the variables in the random part 
              
              # list of data for the stan code 
              pulpdat <- list(n=n,n_RE=n_RE,N1=N1,ncx1=ncx1,id1=id1,RE_ind1=RE_ind1,y1=y1,Zs1=Zs1,
                              scale_betas1=100,scale_sigmas=10,scale_diag_D=10,lkj_shape=2.0, Ztinv1= Ztinv1, Zinv1= Zinv1,Zv1=Zv1,Z_1=Z_1,
                              ZrowsStart1=ZrowsStart1,ZrowsEnd1=ZrowsEnd1,Xs1=Xs1,SDs_X1=SDs_X1,mean_sd_X1=mean_sd_X1,mean_sd_Z1=mean_sd_Z1,SDs_Z1=SDs_Z1)
              
              
              
              # fitting of the model 
              
              fit = stan(file = file.stan, data = pulpdat, chains = 2, warmup = 4000, iter = 10000, thin = 1,
                         control = list(max_treedepth = 20, adapt_delta = 0.95))
              # summary table of the parameters 
              summary1 = as.data.frame(summary(fit, pars = c("betas1","sigma1","D","b"))$summary)
             # we save the summary table 
               save(summary1, file = paste0("V://Users//055609(A_Lika)//HC+NHC_models_simulations//summary_tables_10000//summary_", sim.name, "_d", d, ".RData"))
            }
            
             
             ## We fit the HNC model with centred variables 
             #if the optionHC="HNC" and scaling="centred"
             
            if(((optionHC=="HNC")&(scaling=="centred"))){
              
              
              necessary<-extractFrames_new(y ~  year + sex + drug + ( year | id),data=data)
              
              rstan_options(auto_write = TRUE)
              
              
              n<-necessary$n # total number of subjects 
              n_RE<-2 # total number of random terms
              N1<-necessary$N # total number of observations  
              ncx1<-necessary$ncx # total number of variables in the design matrix of the fixed part
              id1<-necessary$id # the ids of the subjects 
              RE_ind1<-1:2 # vector from 1 to the total number of the random terms 
              y1<-necessary$y # dependent variable
              Z_1<-necessary$Z_ # the design matrix of the random part when we have raw data
              Zc1<-necessary$Zc # the design matrix of the random part when we have centred data
              Zv1<-necessary$Zv # the transposed design matrix of the random part when we have raw data
              Xc1<-necessary$Xc # the design matrix of the fixed part when we have centred data
              
              
              
              Ztinv1<-necessary$Ztinv # the transposed inverse design matrix of the random part of all the subjects when we have raw data
              Zinv1<-necessary$Zinv # the inverse design matrix of the random part of all the subjects when we have raw data
        
              
              means_X1<-necessary$means_X # the means of the variables in the fixed part 
              means_Z1<-necessary$means_Z # the means of the variables in the random part 
              
              
              # list of data required for the stan code 
              
              pulpdat <- list(n=n,n_RE=n_RE,N1=N1,ncx1=ncx1,id1=id1,RE_ind1=RE_ind1,y1=y1,Zc1=Zc1,
                              scale_betas1=100,scale_sigmas=10,scale_diag_D=10,lkj_shape=2.0, Ztinv1= Ztinv1, Zinv1= Zinv1,Zv1=Zv1,Z_1=Z_1,Xc1=Xc1,
                              means_Z1=means_Z1,means_X1=means_X1)
              
              
              ## Fitting the model
              
              fit <- stan(file=file.stan,data=pulpdat,chains = 2, warmup = 4000, iter = 10000, thin = 1,
                          control = list(max_treedepth = 20, adapt_delta = 0.95))
              
              summary1 = as.data.frame(summary(fit, pars = c("betas1","sigma1","D","b"))$summary)
              save(summary1, file = paste0("V://Users//055609(A_Lika)//HC+NHC_models_simulations//summary_tables_10000//summary_", sim.name, "_d", d, ".RData"))
            }
            
            
             ## We fit the HNC model with raw variables 
             #if the optionHC="HNC" and scaling="raw"
             
            if(((optionHC=="HNC")&(scaling=="raw"))){
              
              
              necessary<-extractFrames_new(y ~  year + sex + drug + ( year | id),data=data)
              
              rstan_options(auto_write = TRUE)
              
              
              n<-necessary$n # total number of subjects  
              n_RE<-2 # total number of random terms 
              N1<-necessary$N # total number of observations 
              ncx1<-necessary$ncx # total number of variables in the design matrix of the fixed part 
              id1<-necessary$id # the ids of the subjects 
              RE_ind1<-1:2 # vector from 1 to the total number of the random terms 
              y1<-necessary$y # dependent variable
              Z_1<-necessary$Z_ # the design matrix of the random part when we have raw data 
              
              X1<-necessary$X # the design matrix of the fixed part when we have raw data 
              
              
              
              # the list of data required for the stan code 
              
              pulpdat <- list(n=n,n_RE=n_RE,N1=N1,ncx1=ncx1,id1=id1,RE_ind1=RE_ind1,y1=y1,Z_1=Z_1,
                              scale_betas1=100,scale_sigmas=10,scale_diag_D=10,lkj_shape=2.0,X1=X1)
              
              
              ## Fitting the model
              
              fit <- stan(file=file.stan,data=pulpdat,chains = 2, warmup = 4000, iter = 10000, thin = 1,
                          control = list(max_treedepth = 20, adapt_delta = 0.95))
              
              summary1 = as.data.frame(summary(fit, pars = c("betas1","sigma1","D","b"))$summary)
              save(summary1, file = paste0("V://Users//055609(A_Lika)//HC+NHC_models_simulations//summary_tables_10000//summary_", sim.name, "_d", d, ".RData"))
            }
            
            
            ####################
            #  summary.table   # 
            ####################
       
            # we save the elements elapsed time of the stan code 
            
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "time"] = sum(get_elapsed_time(fit))
            
             # we save the results of every parameter of the stan model   
             
             
             #beta1: for the coefficient of the intercept 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_beta1"] = summary1["betas1[1]",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_beta1"] = summary1["betas1[1]",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_beta1"] = summary1["betas1[1]",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_beta1"] = summary1["betas1[1]",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_beta1"] = summary1["betas1[1]",]$Rhat
            
            #beta2: for the coefficient of the year 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_beta2"] = summary1["betas1[2]",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_beta2"] = summary1["betas1[2]",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_beta2"] = summary1["betas1[2]",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_beta2"] = summary1["betas1[2]",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_beta2"] = summary1["betas1[2]",]$Rhat
            
            #beta3: for the coefficient of the sex 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_beta3"] = summary1["betas1[3]",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_beta3"] = summary1["betas1[3]",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_beta3"] = summary1["betas1[3]",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_beta3"] = summary1["betas1[3]",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_beta3"] = summary1["betas1[3]",]$Rhat
            
            #beta4: for the coefficient of the drug 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_beta4"] = summary1["betas1[4]",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_beta4"] = summary1["betas1[4]",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_beta4"] = summary1["betas1[4]",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_beta4"] = summary1["betas1[4]",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_beta4"] = summary1["betas1[4]",]$Rhat
            
            # sigma1: for the variance of the error terms 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_sigma2"] = summary1["sigma1",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_sigma2"] = summary1["sigma1",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_sigma2"] = summary1["sigma1",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_sigma2"] = summary1["sigma1",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_sigma2"] = summary1["sigma1",]$Rhat
            
            # D11: For the variance of the random intercept 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_D11"] = summary1["D[1,1]",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_D11"] = summary1["D[1,1]",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_D11"] = summary1["D[1,1]",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_D11"] = summary1["D[1,1]",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_D11"] = summary1["D[1,1]",]$Rhat
            
            # D12: for the covariance of the random effects 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_D12"] = summary1["D[1,2]",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_D12"] = summary1["D[1,2]",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_D12"] = summary1["D[1,2]",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_D12"] = summary1["D[1,2]",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_D12"] = summary1["D[1,2]",]$Rhat
            
            # D22: for the variance of the random slope 
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "mean_D22"] = summary1["D[2,2]",]$mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "n_eff_D22"] = summary1["D[2,2]",]$n_eff
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "se_mean_D22"] = summary1["D[2,2]",]$se_mean
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "sd_D22"] = summary1["D[2,2]",]$sd
            summary.table[( summary.table$sigma2 == sigma2 & summary.table$var == var 
                            & summary.table$corr == corr & summary.table$optionHC == optionHC 
                            & summary.table$scaling == scaling & summary.table$simulation == d), "Rhat_D22"] = summary1["D[2,2]",]$Rhat

          }
        }
      }
    }
  }
}

summary.table 

