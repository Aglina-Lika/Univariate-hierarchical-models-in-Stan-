# function to load and get the data of any RData file 

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# we set the working directory that we have saved the RData file and we load it 

setwd("E://stat paper 1")
data = loadRData("sim_total_1_100.RData")


# the 7th column has the running time of the Hamiltonian Monte Carlo in seconds. We transform it in minutes 

data[,7]<-data[,7]/60

# how many times are equal or greater that 60 minutes. These are outliers.  

sum(data$time >=60)



# we remove the records that have convergence time greater or equal to 60 minutes. 

data = data[data$time < 60, ]



# we check how many records remained

length(data[,7]) #4790

# we calculate the ESS(t)/Time for every parameter. The ESS(t) is the effective sample size and 
# the time is the convergence time.  

data$timeEff_intercept = data$n_eff_beta1/data$time
data$timeEff_time = data$n_eff_beta2/data$time
data$timeEff_sex = data$n_eff_beta3/data$time
data$timeEff_drug = data$n_eff_beta4/data$time

data$timeEff_sigma2 = data$n_eff_sigma2/data$time
data$timeEff_D11 = data$n_eff_D11/data$time 
data$timeEff_D12 = data$n_eff_D12/data$time
data$timeEff_D22 = data$n_eff_D22 /data$time

data$condition = paste0("sigma2 = ", data$sigma2, ", varD = ", data$var, ", corr = ", data$corr)





# we will create an new column which will define every time the labels of each box plot based on the 
# combination of the parameters that are used i.e., what was the sigma2, varD and corr.

data$condition2<-NULL
#######################################################################################################################

# sigma2=0.5, varD=0.5, corr=0.1, 
data$condition2[data$condition==tabs_matrix[1]]<-paste0("A:",data$condition[data$condition==tabs_matrix[1]])

# sigma2=0.5, varD=0.5, corr=0.5,
data$condition2[data$condition==tabs_matrix[2]]<-paste0("A:",data$condition[data$condition==tabs_matrix[2]])


# sigma2=0.5, varD=5, corr=0.1, 
data$condition2[data$condition==tabs_matrix[3]]<-paste0("B:",data$condition[data$condition==tabs_matrix[3]])

# sigma2=0.5, varD=5, corr=0.5,
data$condition2[data$condition==tabs_matrix[4]]<-paste0("B:",data$condition[data$condition==tabs_matrix[4]])


# sigma2=5, varD=0.5, corr=0.1, 
data$condition2[data$condition==tabs_matrix[5]]<-paste0("C:",data$condition[data$condition==tabs_matrix[5]])

# sigma2=5, varD=0.5, corr=0.5,
data$condition2[data$condition==tabs_matrix[6]]<-paste0("C:",data$condition[data$condition==tabs_matrix[6]])


# sigma2=5, varD=5, corr=0.1, 
data$condition2[data$condition==tabs_matrix[7]]<-paste0("D:",data$condition[data$condition==tabs_matrix[7]])

# sigma2=5, varD=5, corr=0.5,
data$condition2[data$condition==tabs_matrix[8]]<-paste0("D:",data$condition[data$condition==tabs_matrix[8]])





## We make expressions for the y-axis of the box plots. The expressions will combine the ESS_ with the naming of 
# the coefficients that we havegave in the models 


intercept<-expression(paste("ESS(",beta[0],")"," ","/"," ","Time"))

time<-expression(paste("ESS(",beta[1], ")", " ","/"," ","Time"," ", "and", " ", "ESS(", beta[Z[1]],")", " ","/"," ","Time"))

sex<- expression(paste("ESS(",beta[2], " )", " ","/", " ","Time"," ", "and", " ", "ESS(", beta[1]^-t,")"," ","/"," ","Time"))


drug<-expression(paste("ESS(",beta[3], " )", " ","/"," ", "Time"," ", "and", " ", "ESS(", beta[2]^-t,")"," ","/"," ","Time"))

s_var<-expression(paste("ESS(",sigma^2,")"," ","/"," ","Time"))

d0<-expression(paste("ESS(",sigma[b[0]]^2,")"," ","/"," ","Time"))

d12<-expression(paste("ESS(",sigma[b[12]],")"," ","/"," ","Time"))

d1<-expression(paste("ESS(",sigma[b[1]]^2,")"," ","/"," ","Time"))





# We require the library plyr. 

library(plyr)

# we revalue the variable scaling as raw=1, centre=2 and standardize=3 

data$scaling = revalue(data$scaling, c("raw"= 1, "centre"= 2, "standardize"= 3))



scaling<-data$scaling

# We set with 6 the scaling that has value 3 
scaling[scaling==3]<-6

# we save the new scaling to the data set. The scaling column will be used in the box plots to define the 
# line type of the median 

data$scaling<-scaling 

# we seperate the dataset in 2 datasets based on the value of the parameter corr. 

data_corr0.1 = data[data$corr == 0.1,]
data_corr0.5 = data[data$corr == 0.5,]





# We do the box plots of the ESS(t)/time of every parameter for every parametrization, 
# variable transformation (raw, centre, standardize) and combination of parameters.


# For the sub data sets with cor=0.1 we have:  


library(ggplot2)

# corr = 0.1
ggplot(data_corr0.1, aes(x = optionHC, y = time, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_intercept, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=intercept)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_time, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=time)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_sex, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=sex)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_drug, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=drug)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_sigma2, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=s_var)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_D11, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=d0)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_D12, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=d12)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centered", "standardized"))

ggplot(data_corr0.1, aes(x = optionHC, y = timeEff_D22, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=d1)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))




# For the sub-datasets with cor=0.5 we have:  


ggplot(data_corr0.5, aes(x = optionHC, y = time, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_intercept, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=intercept)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_time, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=time)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_sex, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=sex)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_drug, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=drug)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_sigma2, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=s_var)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_D11, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=d0)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_D12, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=d12)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

ggplot(data_corr0.5, aes(x = optionHC, y = timeEff_D22, linetype = scaling)) +
  geom_boxplot() +
  facet_wrap(~ condition2) +
  labs(y=d1)+
  scale_linetype_discrete(name = "Scaling", labels = c("raw", "centred", "standardized"))

