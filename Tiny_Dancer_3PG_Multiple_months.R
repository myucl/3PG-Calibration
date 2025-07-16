####################################
###Tiny dancer calibration of 3PG###
####################################

source("R/InitialiseRunFunctions.R")
library(SoNWaL)
library(lubridate)
library(coda)
library(BayesianTools)
library(miscTools)
library(ggpubr)
library(matrixStats)
library(future)
library(furrr)
library(parallel)
# Time step to run SonWal with
timeStep<-"weekly"


#get climate data
clm_df_full<-read.csv(paste0("data/clm_df_",timeStep,"_Harwood_CHESS-MET_1973-2024.csv"))

# get calibration data from Harwood
flxdata_daily <- read.csv("data/harwood_daily_2015_2020.csv")


# load in default parameters for Sitka 
sitka<-getParmsSitka(weather=clm_df_full,
                     waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                     timeStp = timeStep)

# run SoNWal with Sitka parameters
output<-do.call(SoNWaL,sitka)

#Check the variables we are working with
colnames(output)


# Name of parameters to calibrate for
nm <- c( "sigma_zR","E_S1", "E_S2", "shared_area", "maxRootDepth", "K_drain", 
         "pFS2", "pFS20", "aS", "nS", "pRx", "pRn", "gammaFx","gammaF0", "tgammaF", 
         "Rttover", "mF","mR","mS","SLA0", "SLA1", "tSLA","k", "alpha", "Y", "totNx",
         "kYmax", "kOmax", "hc", "qi", "qh", "qb", "eY", 
         "fYC", "fYN")


#Lower bound for variables
f.decrease <- c(
  0.1, #sigma_zR
  0.001,#E_S1
  0.01, #E_S2
  1, #shared_area
  0.5, #maxRootDepth
  0.01, #K_drain
  0.3, #pFS2
  0.1, #pFS20
  0.01, #aS
  2, #nS
  0.3, #pRx
  0.2, #PRn
  0.002, #gammaFx
  0.00005, #gammaF0
  20, #tgammaF
  0.003, #Rttover
  0.001, #mF
  0.1, #mR
  0.1, #mS
  4, #SLA0
  2, #SLA1
  2, #tSLA
  0.4, #k
  0.03, #alpha
  0.43, #Y
  2, #totNx
  0.001, #kYmax
  0.00001, #kOmax
  0.01, #hc
  90, #qi
  3, #qh
  0.6, #qb
  0.01, #eY
  0.1,#fYC
  0.1 #fYN
)



#Upper bound for variables
f.increase <- c(
  3, #sigma_zR
  1, #E_S1
  1, #E_S2
  6, #shared_area
  2, #maxRootDepth
  1, #K_drain
  1.6, #pFS2
  1.5, #pFS20
  0.3, #aS
  3, #nS
  1, #pRx
  0.4, #pRn
  0.1, #gammaFx
  0.05, #gammaF0
  110, #tgammaF
  0.16, #Rttover
  0.8, #mF
  0.8, #mR
  0.5, #mS
  9, #SLA0
  6, #SLA1
  10, #tSLA
  0.7, # k
  0.06, #alpha
  0.49, #Y
  20, #totNx
  0.1, #kYmax
  0.01, #kOmax
  0.8, #hc
  600, #qi
  50, #qh
  20, #qb
  0.9, #eY
  0.9, #fYC
  0.9 #fYN
)

#Create a lower and upper bound dataframe
para <- data.frame(
  rbind(f.decrease, f.increase)
)
colnames(para) <- nm

#Set modifier for NPP
modif<- 7.142857



#Set up dgpsi
library(lhs)
library(dgpsi)
init_py(verb = T)
sklearn <- reticulate::import('sklearn')
dgpsi_py <- reticulate::import('dgpsi')



#The subset of variables that are important for NPP calibration
para_sub<-para[,c("shared_area","K_drain","pFS2", "pFS20", "gammaFx", "SLA0", "SLA1", "tSLA", "k", "alpha", "Y", "fYC", "fYN","totNx", "kYmax", "kOmax", "hc", "qi", "qh", "qb", "eY")]

# "pRx", "pRn" are now deleted
#Suzzie's list
#nm <- c( "pFS2", "pFS20", "gammaFx",
#"SLA0", "SLA1", "tSLA",
#"k", "alpha", "Y", "totNx",
#"kYmax", "kOmax", "hc", "qi", "qh", "qb", "eY",
#"fYC", "fYN")
#Initial design of the subsetted variables
library(lhs)
lhd <- lhs::maximinLHS(2000, ncol(para_sub))
for (i in 1:ncol(lhd)){
  min_val <- para_sub[1,i]
  max_val <- para_sub[2,i]
  lhd[,i] <- lhd[,i]*(max_val-min_val) + min_val
}

lhd <- data.frame(lhd)

names(lhd) <- names(para_sub)


#Variables that are set to default values
fixed_names<-setdiff(colnames(para),colnames(para_sub))
para_rest<-para[,c(fixed_names)]
library(data.table)
library(dplyr)


#Getting the initial design 3PG results
all_output<-data.frame()


n <- nrow(lhd)
results <- vector("list", n) # preallocate list for speed

for (i in seq_len(n)) {
  
  sitka[colnames(para_sub)] <- lhd[i, ]
  
  output <- do.call(SoNWaL, sitka)
  output$id <- i
  
  results[[i]] <- output[, c(1,2,3,46,103)]
  
  if (i %% 10 == 0) cat("Iteration:", i, "\n") # reduce frequency of print statements
}


all_output <- rbindlist(results)

summary(all_output$NPP)


all_output<-all_output%>%
  mutate(NPP=NPP*modif)%>%
  group_by(Year, Month,id)%>%
  summarise(
    NPP = mean(NPP),     
  )

summary(all_output$NPP)


#Calibrating 2018 July
all_output_sum<-all_output[all_output$Year%in%c(2016,2018),]
all_output_sum<-na.omit(all_output_sum)

all_output_sum_1<-na.omit(all_output_sum)


plot(all_output_sum_1$NPP)
hist(all_output_sum_1$NPP)
library(dplyr)
all_output_sum_1 <- all_output_sum_1 %>% mutate(YM = paste(Year, Month,sep = ""))

library(ggplot2)
library(tidyr)


all_output_sum_1_20187<-all_output_sum_1[all_output_sum_1$YM=="20187",]
Train_ydfmm_20187<-all_output_sum_1_20187$NPP
Y_20187<-unname(as.matrix(Train_ydfmm_20187))
summary(Y_20187)

all_output_sum_1_20167<-all_output_sum_1[all_output_sum_1$YM=="20167",]
Train_ydfmm_20167<-all_output_sum_1_20167$NPP
Y_20167<-unname(as.matrix(Train_ydfmm_20167))
summary(Y_20167)



lhd2<-lhd
for (i in 1:ncol(lhd2)) {
  min_val <- para_sub[1, i]
  max_val <- para_sub[2, i]
  lhd2[, i] <- (lhd2[, i] - min_val) / (max_val - min_val)
}
Train_xdfmm<-lhd2
X<-Train_xdfmm

X<-unname(as.matrix(X))



#########################
#Building emulator

m_dgp_sub_20187_try <- dgp(X, Y_20187, name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(m_dgp_sub_20187_try, style=1, dim=4)
plot(m_dgp_sub_20187_try, style=1)
summary(m_dgp_sub_20187_try)
write(m_dgp_sub_20187_try, 'm_dgp_sub_20187_V0715_try')


m_dgp_sub_20167_try <- dgp(X, Y_20167, name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(m_dgp_sub_20167_try, style=1, dim=4)
plot(m_dgp_sub_20167_try, style=1)
summary(m_dgp_sub_20167_try)
write(m_dgp_sub_20167_try, 'm_dgp_sub_20167_V0715_try')
##############################################
library(dgpsi)
init_py()
library(tinydancer)
library(lhs)
library(tictoc)



#implausibility calculation
general_imp <- function(NewData, Emulator, Obs, Disc, ObsErr){
  tem <- predict(Emulator,NewData)$results
  abs(Obs-tem$mean)/sqrt(tem$var + Disc + ObsErr)
}

#implausibiility that can handel multiple months
shubert_imp <- function(x, targetLevel, levels, waves, 
                        EmulatorList = NULL,
                        ObsList = NULL, 
                        DiscList = NULL, 
                        ObsErrList = NULL){
  
  # Set defaults if not provided
  if(is.null(ObsList)) ObsList <- rep(5.240009995207778, waves)
  if(is.null(DiscList)) DiscList <- rep(0.01, waves)
  if(is.null(ObsErrList)) ObsErrList <- rep(0, waves)
  
  ans <- rep(Inf, levels)
  waveFail <- FALSE
  this.level <- 1
  wave.num <- 1
  Timp <- NA
  
  while ((this.level <= levels) & !waveFail) {
    Timp <- general_imp(
      NewData = x, 
      Emulator = EmulatorList[[wave.num]], 
      Obs = ObsList[wave.num],
      Disc = DiscList[wave.num],
      ObsErr = ObsErrList[wave.num]
    )
    wave.num <- wave.num + 1
    if ((Timp > targetLevel[this.level]) & (wave.num <= waves)) {
      waveFail <- TRUE
    }
    if ((!waveFail) & (wave.num > waves)) {
      ans[this.level:levels] <- Timp
      this.level <- levels + 1
    } else {
      ans[this.level] <- Timp
      this.level <- this.level + 1
    }
  }
  return(ans)
}


# Create a wrapper function of multiple months implausibility calculation
implausibility_function <- function(x, targetLevel, levels = 1, waves = 1) {
  return(shubert_imp(x, targetLevel, levels, waves, 
                     EmulatorList = your_emulator_list,
                     ObsList = my_obs_values,
                     DiscList = my_disc_values,
                     ObsErrList = my_obserr_values))
}
#Define emulator
gp_em2<-m_dgp_sub_20187_try
gp_em2_1<-m_dgp_sub_20167_try

# Define your parameter lists
EmulatorList <- list(gp_em2,gp_em2_1)
my_obs_values <- c(5.240009995207778, 6.578773935533921)
my_disc_values <- c(0.01, 0.01)
my_obserr_values <- c(0, 0)
#boundaries for variables
lim <- matrix(rep(c(0, 1), 21), nrow = 21, byrow = TRUE)

#Getting design for second wave
tic()
EmulatorList <- list(gp_em2,gp_em2_1)
control_list <- list(
  num_mutations = 10,
  num_iterations = 100,
  box_limits = matrix(rep(c(0, 1), 21), nrow = 21, byrow = TRUE),
  fork_chains=FALSE,
  estimate_volume=TRUE
)
control_list$one_per_level <- TRUE
new_ladder_w1 <- construct_temperature_ladder(
  implausibility = implausibility_function,
  dims = 21,
  target_levels = 0.5,
  control_list = control_list
)

# Define your parameter lists
my_obs_values <- c(5.240009995207778, 6.578773935533921)
my_disc_values <- c(0.01, 0.01)
my_obserr_values <- c(0, 0)

#Using tinydancer to find the design for second wave
result_w1 <- sample_slice_ptmcmc(new_ladder_w1, n_iter = 2000,
                                 implausibility = function(x, target_level) {
                                   shubert_imp(x, targetLevel = target_level, 
                                               levels = 1, waves = 1,
                                               EmulatorList = EmulatorList,
                                               ObsList = my_obs_values,
                                               DiscList = my_disc_values,
                                               ObsErrList = my_obserr_values)
                                 },
                                 control_list = control_list)

wave2_des <- result_w1$uniform_sample_list[[2]][,1:length(result_w1$restart[[1]]$Xt)]
wave2_des<-as.data.frame(wave2_des)
colnames(wave2_des)<-names(para_sub)

#Function to fun 3pg with parameter set
shubert <- function(x) {
  x<-as.data.frame(x)
  for (i in 1:ncol(x)){
    min_val <- para_sub[1,i]
    max_val <- para_sub[2,i]
    x[i] <- x[i]*(max_val-min_val) + min_val
  }
  sitka[c(colnames(para_sub))]<-x
  #sitka[c(fixed_names)]<-para_rest[1,]
  output<-do.call(SoNWaL,sitka)
  output<-output%>%
    mutate(NPP=NPP*modif)%>%
    group_by(Year, Month) %>%
    summarise(NPP = mean(NPP))
  output <- output %>% mutate(YM = paste(Year, Month,sep = ""))
  output<-output[output$YM%in%c("20187","20167"),]
  return(output$NPP)
}



#Running 3pg for second wave training data. 
wave2<-data.frame()

for (j in 1:nrow(wave2_des)){
  output<-shubert(wave2_des[j,])
  #output$ID<-j
  wave2<-rbind(wave2,output)
  print(j)
}

#Build emulator for second wave
gp_wave2_mat <- dgp(as.matrix(wave2_des), as.matrix(wave2[,1]), name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(gp_wave2_mat)

gp_wave2_mat_1 <- dgp(as.matrix(wave2_des), as.matrix(wave2[,2]), name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(gp_wave2_mat_1)


EmulatorList <- list(gp_em2, gp_em2_1, gp_wave2_mat, gp_wave2_mat_1, gp_wave2_mat, gp_wave2_mat_1)
my_obs_values <- c(5.240009995207778, 6.578773935533921,5.240009995207778, 6.578773935533921,5.240009995207778, 6.578773935533921)
my_disc_values <- c(0.01, 0.01,0.01, 0.01,0.01, 0.01)
my_obserr_values <- c(0, 0,0, 0,0, 0)
control_list$one_per_level <- TRUE
implausibility_function <- function(x, targetLevel, levels = 6, waves = 6) {
  return(shubert_imp(x, targetLevel, levels, waves, 
                     EmulatorList = EmulatorList,
                     ObsList = my_obs_values,
                     DiscList = my_disc_values,
                     ObsErrList = my_obserr_values))
}
new_ladder_w2 <- construct_temperature_ladder(
  implausibility = implausibility_function,
  dims = 21,
  target_levels = c(0.5,0.5,0.5,0.5,0.25,0.25),
  control_list = control_list
)


#new_ladder_w2 <- construct_temperature_ladder(
  #implausibility = function(x, target_level){shubert_imp(x, target_level, 3,3)},
  #dims = 23,
  #target_levels = c(0.5,0.5,0.25),
  #control_list = control_list
#)
result_w2 <- sample_slice_ptmcmc(new_ladder_w2, n_iter = 2000,
                                 implausibility = function(x, target_level) {
                                   shubert_imp(x, targetLevel = target_level, 
                                               levels = 6, waves = 6,
                                               EmulatorList = EmulatorList,
                                               ObsList = my_obs_values,
                                               DiscList = my_disc_values,
                                               ObsErrList = my_obserr_values)
                                 },
                                 control_list = control_list)
#NROY wave 3
w3_chain_no <- dplyr::first(which(lapply(new_ladder_w2$imp_levels, function(j) j[2]) <= 3)) + 1
nroyw2_des <- result_w2$uniform_sample_list[[w3_chain_no]][,1:length(result_w2$restart[[1]]$Xt)]
#values very close to
nchains <- length(result_w2$uniform_sample_list)
wave3_des_imp0.5 <- result_w2$uniform_sample_list[[nchains]][,1:length(result_w2$restart[[1]]$Xt)]
#100 NROY and 100 close
wave3_des <- rbind(nroyw2_des[101:2000,], wave3_des_imp0.5[1001:2000,])


wave3_des<-as.data.frame(wave3_des)
colnames(wave3_des)<-names(para_sub)
wave3<-data.frame()

for (j in 1:nrow(wave3_des)){
  output<-shubert(wave3_des[j,])
  #output$ID<-j
  wave3<-rbind(wave3,output)
  print(j)
}


wave5<-data.frame()

f_w2 <- function(x) {
  x<-as.data.frame(x)
  for (i in 1:ncol(x)){
    min_val <- para_sub[1,i]
    max_val <- para_sub[2,i]
    x[i] <- x[i]*(max_val-min_val) + min_val
  }
  sitka[c(colnames(para_sub))]<-x
  #sitka[c(fixed_names)]<-para_rest[1,]
  output<-do.call(SoNWaL,sitka)
  output<-output%>%
    mutate(NPP=NPP*modif)%>%
    group_by(Year, Month) %>%
    summarise(NPP = mean(NPP))
  output <- output %>% mutate(YM = paste(Year, Month,sep = ""))
  #output<-output[output$YM=="20183",]
  return(as.data.frame(output[,c("NPP","YM")]))
}

for (j in 2400:nrow(wave3_des)){
  output<-f_w2(wave3_des[j,])
  #output$ID<-j
  wave5<-rbind(wave5,output)
  print(j)
}
wave5$ID<-rep(1:501, each = 625)

filtered_output_w5 <- wave5 %>%
  filter(grepl("^2015|^2016|^2017|^2018", YM))
filtered_output_w5 <- filtered_output_w5 %>%
  mutate(
    YM = as.character(YM),  # Convert to character
    year = substr(YM, 1, 4),    # Extract year
    month = substr(YM, 5, nchar(YM)), # Extract month
    month = sprintf("%02d", as.numeric(month)), # Ensure two-digit month
    date = ymd(paste0(year, "-", month, "-01")) # Convert to date format
  ) %>%
  select(-year, -month)



startYear = 2015
endYear = 2018
observedVals_npp<-function(timeStep,data,sY=2015,eY=2018){
  
  data<-filter(data,year>=sY&year<=eY)
  data<-data[-1,]
  if(timeStep =="weekly"){
    data<-data%>%
      mutate(grp=week(as.Date(data$yday, origin = paste0(data$year,"-01-01"))))
  }
  
  if(timeStep =="monthly"){
    data<-data%>%
      mutate(grp=month(as.Date(data$yday, origin = paste0(data$year,"-01-01"))))
    
  }
  
  sdMin<-data%>% group_by(year,grp) %>%
    dplyr::summarise(sdnpp=mean(npp))
  
  observed <- c(pull(data%>% 
                       dplyr::group_by(year,grp) %>%
                       dplyr::summarise(npp=mean(npp))%>%
                       dplyr::select(npp))
  )
  
  coefVar1=0.2
  dev <- c(
    sapply(1:length(sdMin$sdnpp), function(i) max( coefVar1* abs(sdMin$sdnpp[i]),0.01)))
  return(list(observed,dev))
}
observed<-observedVals_npp(timeStep = "monthly",data=flxdata_daily, # we are calibrating against monthly means 
                           sY=startYear,eY=endYear)[[1]]
observed_NPP<-data.frame(observed)
observed_NPP$NPP<-observed_NPP$observed
observed_NPP$observed<-NULL

observed_NPP$date<-filtered_output_w5$date[1:48]

observed_NPP$ID<-"observed"

observed_NPP$YM<-1
observed_NPP <- observed_NPP[, colnames(filtered_output_w5)]

library(ggplot2)


ggplot(filtered_output_w5, aes(x = date, y = NPP, colour = ID, group = ID)) + 
  geom_line(alpha=0.3) +
  geom_line(data = observed_NPP, aes(x = date, y = NPP), color = "red", linewidth = 2) +
  theme(legend.position = "none")


ggplot(filtered_output_w5[filtered_output_w5$ID==6,], aes(x = date, y = NPP, colour = ID, group = ID)) + 
  geom_line(alpha=0.3) +
  geom_line(data = observed_NPP, aes(x = date, y = NPP), color = "red", linewidth = 2) +
  theme(legend.position = "none")


# Step 1: Merge the dataframes by Year
merged_df <- merge(filtered_output_w5, observed_NPP, by = "date")

# Step 2: Calculate MSE for each lai_bal0_nt group
library(dplyr)

mse_by_id <- merged_df %>%
  group_by(ID.x) %>%
  summarise(
    mse = mean((NPP.x - NPP.y)^2, na.rm = TRUE)
  ) %>%
  arrange(mse)

# Step 3: Get the ID with the lowest MSE
best_id <- mse_by_id %>% slice(1)







