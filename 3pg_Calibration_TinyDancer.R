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
para_sub<-para[,c("shared_area","K_drain","pFS2", "pFS20", "gammaFx", "SLA0", "SLA1", "tSLA", "k", "alpha", "Y", "pRx", "pRn", "fYC", "fYN","totNx", "kYmax", "kOmax", "hc", "qi", "qh", "qb", "eY")]


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

# Bind once at the end (very fast)
all_output <- rbindlist(results)
#plot(all_output$NPP)


summary(all_output$NPP)


all_output<-all_output%>%
  mutate(NPP=NPP*modif)%>%
  group_by(Year, Month,id)%>%
  summarise(
    NPP = mean(NPP),     
  )

summary(all_output$NPP)

all_output_sum<-all_output[all_output$Year==2018,]
all_output_sum<-na.omit(all_output_sum)
#all_output_sum_1<-all_output_sum%>%
#mutate(NPP=NPP*modif)%>%
#group_by(Year, Month,id)%>%
#summarise(
#NPP = mean(NPP),     
#)
all_output_sum_1<-na.omit(all_output_sum)


plot(all_output_sum_1$NPP)
hist(all_output_sum_1$NPP)
library(dplyr)
all_output_sum_1 <- all_output_sum_1 %>% mutate(YM = paste(Year, Month,sep = ""))

library(ggplot2)
library(tidyr)

# Assume your data is in `df` and the binary outcome is `Y`
df <- as.data.frame(lhd)     # ensure X is a data frame
df$Y <- as.factor(new_Y_20187)

# Pivot longer for ggplot
long_df <- pivot_longer(df, cols = -Y, names_to = "Variable", values_to = "Value")

# Density plots
ggplot(long_df, aes(x = Value, fill = Y)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal()
lhd2<-lhd
for (i in 1:ncol(lhd2)) {
  min_val <- para_sub[1, i]
  max_val <- para_sub[2, i]
  lhd2[, i] <- (lhd2[, i] - min_val) / (max_val - min_val)
}
Train_xdfmm<-lhd2
X<-Train_xdfmm
#X<-Train_xdfmm
X<-unname(as.matrix(X))
#########################
all_output_sum_1_20187<-all_output_sum_1[all_output_sum_1$YM=="20187",]
Train_ydfmm_20187<-all_output_sum_1_20187$NPP
Y_20187<-unname(as.matrix(Train_ydfmm_20187))
summary(Y_20187)
new_Y_20187 <- as.integer(Y_20187 >= 2.5)
m_gp_class_20187 <- dgp(X, new_Y_20187, depth=2, likelihood='Categorical', training=T, vecchia=T, share=F, name='matern2.5')
m_gp_class_20187 <- validate(m_gp_class_20187)
plot(m_gp_class_20187, style=2)
write(m_gp_class_20187, 'm_gp_class_20187_V0708')



##########################

sub_Y_20187 <- Y_20187[Y_20187 >= 2.5,]
sub_X_20187 <- X[Y_20187 >= 2.5,]
m_gp_sub_20187 <- gp(sub_X_20187, sub_Y_20187, name='matern2.5', vecchia=T, verb=F, nugget_est=T, prior='ref')
plot(m_gp_sub_20187, style=1)
m_gp_sub_20187$specs$lengthscales

m_dgp_sub_20187 <- dgp(sub_X_20187, sub_Y_20187, name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(m_dgp_sub_20187, style=1, dim=4)
plot(m_dgp_sub_20187, style=1)
summary(m_dgp_sub_20187)
write(m_dgp_sub_20187, 'm_dgp_sub_20187_V0708')

m_dgp_sub_20187_try <- dgp(X, Y_20187, name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(m_dgp_sub_20187_try, style=1, dim=4)
plot(m_dgp_sub_20187_try, style=1)
summary(m_dgp_sub_20187_try)
write(m_dgp_sub_20187_try, 'm_dgp_sub_20187_V0708_try')



##############################################

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
  output<-output[output$YM=="20187",]
  return(output$NPP)
}

library(dgpsi)
init_py()
library(tinydancer)
library(lhs)
library(tictoc)


general_imp <- function(NewData, Emulator, Obs, Disc, ObsErr){
  tem <- predict(Emulator,NewData)$results
  abs(Obs-tem$mean)/sqrt(tem$var + Disc + ObsErr)
}
shubert_imp <- function(x, targetLevel, levels, waves){
  ans <- rep(Inf, levels)
  waveFail <- FALSE
  this.level <- 1
  wave.num <- 1
  Timp <- NA
  while ((this.level <= levels) & !waveFail) {
    Timp <- general_imp(NewData = x, Emulator=EmulatorList[[wave.num]], Obs=5.240009995207778, Disc=0.01, ObsErr=0)
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

gp_em2<-m_dgp_sub_20187_try

library(tictoc)
tic()
source("ei.R")
lim <- matrix(rep(c(0, 1), 23), nrow = 23, byrow = TRUE)
gp_em2_opt <- design(gp_em2, N = 5, limits = lim, f = shubert, method = ei, eval = opt_monitor)
library(patchwork)
library(ggplot2)
draw(gp_em2_opt) +
  plot_annotation(title = 'Bayesian Optimization Tracker') +
  labs(y = "Minimum Value") +
  # Add a horizontal line to represent the global minimum for benchmarking
  geom_hline(
    aes(yintercept = y, linetype = "Global Minimum"), # Global minimum
    data = data.frame(y = 5.240009995207778),
    color = "#E31A1C",
    size = 0.5
  ) +
  scale_linetype_manual(
    values = c("Global Minimum" = "dashed"),
    name = "" # Remove the legend title
  )
toc()
#min location
gp_em2_opt$data$X[which.min(gp_em2_opt$data$Y),]
#minimum found
gp_em2_opt$data$Y[which.min(gp_em2_opt$data$Y)]


tic()
EmulatorList <- list(gp_em2)
control_list <- list(
  num_mutations = 10,
  num_iterations = 100,
  box_limits = matrix(rep(c(0, 1), 23), nrow = 23, byrow = TRUE),
  fork_chains=FALSE,
  estimate_volume=TRUE
)
new_ladder_w1 <- construct_temperature_ladder(
  implausibility = function(x, target_level){shubert_imp(x, target_level, 1,1)},
  dims = 23,
  target_levels = 0.5,
  control_list = control_list
)

result_w1 <- sample_slice_ptmcmc(new_ladder_w1, n_iter = 2000,
                                 implausibility=function(x, target_level){shubert_imp(x, target_level, 1,1)},
                                 control_list = control_list)

wave2_des <- result_w1$uniform_sample_list[[2]][,1:length(result_w1$restart[[1]]$Xt)]
#wave2 <- shubert(wave2_des)

wave2_des<-as.data.frame(wave2_des)
colnames(wave2_des)<-names(para_sub)
wave2<-data.frame()

for (j in 1:nrow(wave2_des)){
  output<-shubert(wave2_des[j,])
  #output$ID<-j
  wave2<-rbind(wave2,output)
  print(j)
}



gp_wave2_mat <- dgp(as.matrix(wave2_des), wave2$X2.3969934896736e.06, name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(gp_wave2_mat)



EmulatorList <- list(gp_em2, gp_wave2_mat, gp_wave2_mat)

control_list$one_per_level <- TRUE
new_ladder_w2 <- construct_temperature_ladder(
  implausibility = function(x, target_level){shubert_imp(x, target_level, 3,3)},
  dims = 23,
  target_levels = c(0.5,0.5,0.25),
  control_list = control_list
)
result_w2 <- sample_slice_ptmcmc(new_ladder_w2, n_iter = 1900,
                                 implausibility=function(x, target_level){shubert_imp(x, target_level, 3,3)},
                                 control_list = control_list)
#NROY wave 3
w3_chain_no <- dplyr::first(which(lapply(new_ladder_w2$imp_levels, function(j) j[2]) <= 3)) + 1
nroyw2_des <- result_w2$uniform_sample_list[[w3_chain_no]][,1:length(result_w2$restart[[1]]$Xt)]
#values very close to
nchains <- length(result_w2$uniform_sample_list)
wave3_des_imp0.5 <- result_w2$uniform_sample_list[[nchains]][,1:length(result_w2$restart[[1]]$Xt)]
#100 NROY and 100 close
wave3_des <- rbind(nroyw2_des[901:1900,], wave3_des_imp0.5[901:1900,])


wave3_des<-as.data.frame(wave3_des)
colnames(wave3_des)<-names(para_sub)
wave3<-data.frame()

for (j in 1:nrow(wave3_des)){
  output<-shubert(wave3_des[j,])
  #output$ID<-j
  wave3<-rbind(wave3,output)
  print(j)
}


gp_wave3_mat <- dgp(as.matrix(wave3_des), wave3$X5.18941129190778, name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(gp_wave3_mat)

EmulatorList <- list(gp_em2, gp_wave2_mat, gp_wave3_mat, gp_wave3_mat)
new_ladder_w3 <- construct_temperature_ladder(
  implausibility = function(x, target_level){shubert_imp(x, target_level, 4,4)},
  dims = 23,
  target_levels = c(0.5,0.5,0.25,0.15),
  control_list = control_list
)
result_w3 <- sample_slice_ptmcmc(new_ladder_w3, n_iter = 2000,
                                 implausibility=function(x, target_level){shubert_imp(x, target_level, 4,4)},
                                 control_list = control_list)
#NROY wave 3
w4_chain_no <- dplyr::first(which(lapply(new_ladder_w3$imp_levels, function(j) j[2]) <= 3)) + 1
nroyw3_des <- result_w3$uniform_sample_list[[w4_chain_no]][,1:length(result_w3$restart[[1]]$Xt)]
#values very close to
nchains <- length(result_w3$uniform_sample_list)
wave4_des_imp0.5 <- result_w3$uniform_sample_list[[nchains]][,1:length(result_w3$restart[[1]]$Xt)]
#100 NROY and 100 close
wave4_des <- rbind(nroyw3_des[1001:2000,], wave4_des_imp0.5[1001:2000,])
wave4_des<-as.data.frame(wave4_des)
colnames(wave4_des)<-names(para_sub)
wave4<-data.frame()

for (j in 1:nrow(wave4_des)){
  output<-shubert(wave4_des[j,])
  #output$ID<-j
  wave4<-rbind(wave4,output)
  print(j)
}
hist(wave4$X8.52936324640963)


nchains <- length(result_w3$uniform_sample_list)
wave4_des_imp0.5 <- result_w3$uniform_sample_list[[nchains]][,1:length(result_w3$restart[[1]]$Xt)]
#100 NROY and 100 close
wave4_des_1 <- wave4_des_imp0.5
wave4_des_1<-as.data.frame(wave4_des_1)
colnames(wave4_des_1)<-names(para_sub)
wave4_1<-data.frame()

for (j in 1:nrow(wave4_des_1)){
  output<-shubert(wave4_des_1[j,])
  #output$ID<-j
  wave4_1<-rbind(wave4_1,output)
  print(j)
}
hist(wave4_1$X5.15100744494931)


gp_wave4_mat <- dgp(as.matrix(wave4_des), wave4$X8.52936324640963, name=c('matern2.5','sexp'), vecchia=T, verb=T,node=1,share=F)
plot(gp_wave4_mat)

EmulatorList <- list(gp_em2, gp_wave2_mat, gp_wave3_mat, gp_wave4_mat,gp_wave4_mat)
new_ladder_w4 <- construct_temperature_ladder(
  implausibility = function(x, target_level){shubert_imp(x, target_level, 5,5)},
  dims = 23,
  target_levels = c(0.5,0.5,0.25,0.15,0.1),
  control_list = control_list
)
result_w4 <- sample_slice_ptmcmc(new_ladder_w4, n_iter = 2000,
                                 implausibility=function(x, target_level){shubert_imp(x, target_level, 5,5)},
                                 control_list = control_list)
#NROY wave 3
w5_chain_no <- dplyr::first(which(lapply(new_ladder_w4$imp_levels, function(j) j[2]) <= 3)) + 1
nroyw4_des <- result_w4$uniform_sample_list[[w5_chain_no]][,1:length(result_w4$restart[[1]]$Xt)]
#values very close to
nchains <- length(result_w4$uniform_sample_list)
wave5_des_imp0.5 <- result_w4$uniform_sample_list[[nchains]][,1:length(result_w4$restart[[1]]$Xt)]
#100 NROY and 100 close
wave5_des <- rbind(nroyw4_des[1001:2000,], wave5_des_imp0.5[1001:2000,])
wave5_des<-as.data.frame(wave5_des)
colnames(wave5_des)<-names(para_sub)
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

for (j in 1:nrow(wave5_des)){
  output<-f_w2(wave5_des[j,])
  #output$ID<-j
  wave5<-rbind(wave5,output)
  print(j)
}
wave5$ID<-rep(1:2000, each = 625)

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
observed_NPP <- observed_NPP[, colnames(filtered_output_w2)]

library(ggplot2)


ggplot(filtered_output_w5, aes(x = date, y = NPP, colour = ID, group = ID)) + 
  geom_line(alpha=0.3) +
  geom_line(data = observed_NPP, aes(x = date, y = NPP), color = "red", linewidth = 2) +
  theme(legend.position = "none")


ggplot(filtered_output_w5[filtered_output_w5$ID==1637,], aes(x = date, y = NPP, colour = ID, group = ID)) + 
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
