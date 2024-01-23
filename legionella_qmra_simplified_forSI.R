####-----Legionella Water Age - Premise Plumbing - Shower QMRA

set.seed(123) #random number generator
n<-10000 #number of mc samplings
library(growthrates)
library(truncnorm)
library(fitdistrplus)


#---Water Ages
#Water ages are input in a csv with the headers "Scenario" with numbers 1 through 8, "Time" referring to total hours and "WaterAge" referring to water age at shower use. Water Ages were generated in Clements et al. 2023 (Submitted)
wa<-read.csv("scenarios_52w.csv",header=TRUE) # shower

#Remove half of the scenario 2 water ages to limit it to one person living in a two person household. Use this code if you have water ages generated for multi-family households to examine risk to a single occupant.

# Subset the dataframe to include only rows with Scenario 2
subset_df <- wa[wa$Scenario == 2, ]
# Create a sequence of row indices to keep every other row (because two occupants. If three, it would be every third for instance)
rows_to_keep <- seq(from = 1, to = nrow(subset_df), by = 2)
# Subset the data using the row indices to keep
final_df <- subset_df[rows_to_keep, ]
# Combine the filtered rows with rows from other scenarios
wa <- rbind(final_df, wa[wa$Scenario != 2, ])

### WATER AGE SUMMARY

result_table <- wa %>%
  group_by(Scenario) %>%
  summarise(
    Average = mean(WaterAge),
    Median = median(WaterAge),
    Percentile_99 = quantile(WaterAge, 0.99),
    Maximum = max(WaterAge),
    number = length(WaterAge)
  )

# Display the result table
print(result_table)


#Scenario names and definitions
#Scenario 1, Scheduled purging, one person, no purging <-CONTROL
#Scenario 2, Scheduled purging, two people, no purging <-CONTROL
#Scenario 3, Scheduled purging, one person, purge every 12 hours @ 4am/pm
#Scenario 4, Scheduled purging, one person, purge every 24 hours @ 4am
#Scenario 5, Smart purging, one person, purge every 12 hours using 1 liter - both faucets
#Scenario 6, Smart purging, one person, purge every 12 hours using 0.25 liters - both faucets
#Scenario 7, Smart purging, one person, purge every 12 hours using 4 liters - both faucets
#Scenario 8, Smart purging, one person, purge every 24 hours using 1 liter - both faucets


#---Legionella starting concentration
C0<-1 #GC/L 

#---Legionella growth rates
#https://www.frontiersin.org/journals/water/articles/10.3389/frwa.2023.1114795/full

#stagnant, phase 2, PEX, hot water
mumax_hot<-0.88/24 #per day convert to hours
Nmax_hot<-10^5.94 #log10 copies/L convert to copies/L

#stagnant, phase 2, PEX, cold water
mumax_cold<-1.26/24 #per day convert to hours
Nmax_cold<-10^6.12 #log10 copies/L convert to copies/L

#ratio of hot/cold during shower 
ratio<-0.5 #so the hot and cold water is mixed at a 50/50 ratio in the shower

##-------Exposure model-------------------------

#- Kerry Hamilton - aerosol size profile approach----
# thanks https://github.com/mlee267 for conversion to R code help


B <- runif(n, min = 0.013, max = 0.017)# breathing rate, m^3 air/min 

tsh_start<-1 ## We use 1 minute shower length times, because the growth rate/water age is only valid for the time it takes for water to be drawn from the distribution system, which is <1 min. Using 1 min for conservatism. This can be adjusted depending on size of the house.

tsh_final <- rtruncnorm(n, a = 1, mean = 7.8, sd = 0.02) # shower time exposure duration, truncated normal distribution, min. This is a normal shower duration.


# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- rlnorm(n, meanlog = 17.5, sdlog = 0.30)
Caer2 <- rlnorm(n, meanlog = 17.5, sdlog = 0.17)
Caer3 <- rlnorm(n, meanlog = 19.4, sdlog = 0.35)
Caer4 <- rlnorm(n, meanlog = 19.4, sdlog = 0.35)
Caer5 <- rlnorm(n, meanlog = 19.4, sdlog = 0.35)
Caer6 <- rlnorm(n, meanlog = 20.0, sdlog = 0.31)
Caer7 <- rlnorm(n, meanlog = 20.0, sdlog = 0.31)
Caer8 <- rlnorm(n, meanlog = 20.0, sdlog = 0.31)
Caer9 <- rlnorm(n, meanlog = 20.0, sdlog = 0.31)
Caer10 <- rlnorm(n, meanlog = 20.0, sdlog = 0.31)


# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)

# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100

# Deposition efficiency of aerosols by size category
D1 <- runif(n, min = 0.23, max = 0.25)
D2 <- runif(n, min = 0.40, max = 0.53)
D3 <- runif(n, min = 0.36, max = 0.62)
D4 <- runif(n, min = 0.29, max = 0.61)
D5 <- runif(n, min = 0.19, max = 0.52)
D6 <- runif(n, min = 0.10, max = 0.40)
D7 <- runif(n, min = 0.06, max = 0.29)
D8 <- runif(n, min = 0.03, max = 0.19)
D9 <- runif(n, min = 0.01, max = 0.12)
D10 <- runif(n, min = 0.01, max = 0.06)

# Total aerosol concentration, m^3 water/m^3 air (conventional fixture)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)

# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)

# Equivalent inhaled aerosol volume, L (conventional fixture)
wat_conv_start <- B * tsh_start * (CVaer_conv * 1e3) * FD #CVaer_conv is in m3 so multiply by 1e3 to convert to L #BEGINNING OF SHOWER
wat_conv_final <- B * (tsh_final-tsh_start) * (CVaer_conv * 1e3) * FD #CVaer_conv is in m3 so multiply by 1e3 to convert to L #REST OF SHOWER



#log10GC -> infectious units currently 1 GC = 1 IU. Conservative assumption. 
gciu<-function(gc){
  iu<-gc/1
}


##--Dose response models

#-infection endpoint
kinf<-rlnorm(n=n,meanlog = -2.93, sdlog = 0.49)

#-clinically severe illness endpoint
kcsi<-rlnorm(n=n, meanlog = -9.69,sdlog=0.3)

#exponential dose response function
P<-function(k,d){
  p<-1-exp(-k*d)
  p
}

###--Use growth model to calculate a concentration for each water age, mixing concentrations from hot/cold water at a fixed ratio
wa_mix<-wa%>%
  group_by(Scenario) %>%
  mutate(Conc = ratio*(grow_baranyi(WaterAge, c(y0=C0, mumax=mumax_hot, K=Nmax_hot,h0=0))[,"y"])+
           (1-ratio)*(grow_baranyi(WaterAge, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,"y"])) #Concentration in GC/L different growth rates for different temperatures



##-------QMRA-----------------------------------

results_df <- data.frame() #initialize a results dataframe

# Loop through each row of the original dataframe. This code is computationally intensive.

for (i in 1:nrow(wa_mix)) {
  # Multiply Conc by each value in wat_conv
  dose <- as.numeric(wa_mix[i, "Conc"]) * wat_conv_start+C0*wat_conv_final #assuming that for the rest of the shower the concentration is constant at C0
  pinf <- P(kinf,gciu(dose))
  pcsi<-P(kcsi,gciu(dose))
  
  # Create a temporary dataframe with the results
  temp_df <- data.frame(
    Scenario = rep(as.numeric(wa_mix[i, "Scenario"]), n),
    WaterAge = rep(as.numeric(wa_mix[i, "WaterAge"]), n),
    # Time = rep(as.numeric(wa_mix[i, "Time"]), n),
    #Conc = rep(as.numeric(wa_mix[i, "Conc"]), n),
    # Dose = dose,
    Pinf = pinf, 
    Pcsi = pcsi,
    i=i
  )
  
  
  # Append the temporary dataframe to the results dataframe
  results_df <- rbind(results_df, temp_df)
  
  #print(round(i/nrow(wa_mix)*100))
}

results_df$id<-rep(1:n,times=length(unique(results_df$i)))



####ALL RESULTS, Pinf
pinf_sc1 <- dcast(data = filter(results_df, Scenario == 1), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")
pinf_sc2 <- dcast(data = filter(results_df, Scenario == 2), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")
pinf_sc3 <- dcast(data = filter(results_df, Scenario == 3), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")
pinf_sc4 <- dcast(data = filter(results_df, Scenario == 4), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")
pinf_sc5 <- dcast(data = filter(results_df, Scenario == 5), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")
pinf_sc6 <- dcast(data = filter(results_df, Scenario == 6), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")
pinf_sc7 <- dcast(data = filter(results_df, Scenario == 7), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")
pinf_sc8 <- dcast(data = filter(results_df, Scenario == 8), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pinf")

####ALL RESULTS, Pcsi
pcsi_sc1 <- dcast(data = filter(results_df, Scenario == 1), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")
pcsi_sc2 <- dcast(data = filter(results_df, Scenario == 2), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")
pcsi_sc3 <- dcast(data = filter(results_df, Scenario == 3), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")
pcsi_sc4 <- dcast(data = filter(results_df, Scenario == 4), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")
pcsi_sc5 <- dcast(data = filter(results_df, Scenario == 5), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")
pcsi_sc6 <- dcast(data = filter(results_df, Scenario == 6), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")
pcsi_sc7 <- dcast(data = filter(results_df, Scenario == 7), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")
pcsi_sc8 <- dcast(data = filter(results_df, Scenario == 8), formula = id ~ WaterAge, fun.aggregate = sum, value.var = "Pcsi")


#annual infection

Pinf_an <- data.frame(
  Scenario1 = 1 - apply(pinf_sc1[, -which(names(pinf_sc1) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario2 = 1 - apply(pinf_sc2[, -which(names(pinf_sc2) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario3 = 1 - apply(pinf_sc3[, -which(names(pinf_sc3) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario4 = 1 - apply(pinf_sc4[, -which(names(pinf_sc4) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario5 = 1 - apply(pinf_sc5[, -which(names(pinf_sc5) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario6 = 1 - apply(pinf_sc6[, -which(names(pinf_sc6) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario7 = 1 - apply(pinf_sc7[, -which(names(pinf_sc7) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario8 = 1 - apply(pinf_sc8[, -which(names(pinf_sc8) == 'id')], 1, function(row) { prod(1 - row)})
  
)

#annual csi
Pcsi_an <- data.frame(
  Scenario1 = 1 - apply(pcsi_sc1[, -which(names(pcsi_sc1) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario2 = 1 - apply(pcsi_sc2[, -which(names(pcsi_sc2) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario3 = 1 - apply(pcsi_sc3[, -which(names(pcsi_sc3) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario4 = 1 - apply(pcsi_sc4[, -which(names(pcsi_sc4) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario5 = 1 - apply(pcsi_sc5[, -which(names(pcsi_sc5) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario6 = 1 - apply(pcsi_sc6[, -which(names(pcsi_sc6) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario7 = 1 - apply(pcsi_sc7[, -which(names(pcsi_sc7) == 'id')], 1, function(row) { prod(1 - row)}),
  Scenario8 = 1 - apply(pcsi_sc8[, -which(names(pcsi_sc8) == 'id')], 1, function(row) { prod(1 - row)})
  
)

#DALY infecetion
DALY_inf <- data.frame(
  Scenario1 = 0.97*Pinf_an$Scenario1,
  Scenario2 = 0.97*Pinf_an$Scenario2,
  Scenario3 = 0.97*Pinf_an$Scenario3,
  Scenario4 = 0.97*Pinf_an$Scenario4,
  Scenario5 = 0.97*Pinf_an$Scenario5,
  Scenario6 =0.97*Pinf_an$Scenario6 ,
  Scenario7 = 0.97*Pinf_an$Scenario7,
  Scenario8 = 0.97*Pinf_an$Scenario8
  
)

#DALY csi
DALY_csi <- data.frame(
  Scenario1 = 0.97*Pcsi_an$Scenario1,
  Scenario2 = 0.97*Pcsi_an$Scenario2,
  Scenario3 = 0.97*Pcsi_an$Scenario3,
  Scenario4 = 0.97*Pcsi_an$Scenario4,
  Scenario5 = 0.97*Pcsi_an$Scenario5,
  Scenario6 =0.97*Pcsi_an$Scenario6 ,
  Scenario7 = 0.97*Pcsi_an$Scenario7,
  Scenario8 = 0.97*Pcsi_an$Scenario8
  
)

##-------Summary of Results--------------------
summary(Pinf_an) ## annual probability of infection
summary(Pcsi_an) ## annual probability of clinically severe illness
summary(DALY_inf)## DALY for infection endpoint
summary(DALY_csi)## DALY for csi endpoint

##-------Water ages vs Concentration, fig 3---
wa_sim<-runif(n,0,200)#130

pcsi_sim<-P(kcsi,gciu(((ratio*(grow_baranyi(wa_sim, c(y0=C0, mumax=mumax_hot, K=Nmax_hot,h0=0))[,"y"])+
                          (1-ratio)*(grow_baranyi(wa_sim, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,"y"]))
                      ))*wat_conv_start+C0*wat_conv_final)

pinf_sim<-P(kinf,gciu(((ratio*(grow_baranyi(wa_sim, c(y0=C0, mumax=mumax_hot, K=Nmax_hot,h0=0))[,"y"])+
                          (1-ratio)*(grow_baranyi(wa_sim, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,"y"]))
                      ))*wat_conv_start+C0*wat_conv_final)

wa_v_pinfcsi<-data.frame(wa_sim,pinf_sim,pcsi_sim)

summary(wa_v_pinfcsi)

##-------Sensitivity analysis-------------------

###Set all water ages to lognormal distributions
# Scenario 1
scenario_1_data <- wa %>% filter(Scenario == 1)
fit_1 <- fitdist(scenario_1_data$WaterAge, "lnorm", method = "mle")
waterage_1 <- rlnorm(n, meanlog = fit_1$estimate[1], sdlog = fit_1$estimate[2])

# Scenario 2
scenario_2_data <- wa %>% filter(Scenario == 2)
fit_2 <- fitdist(scenario_2_data$WaterAge, "lnorm", method = "mle")
waterage_2 <- rlnorm(n, meanlog = fit_2$estimate[1], sdlog = fit_2$estimate[2])

# Scenario 3
scenario_3_data <- wa %>% filter(Scenario == 3)
fit_3 <- fitdist(scenario_3_data$WaterAge, "lnorm", method = "mle")
waterage_3 <- rlnorm(n, meanlog = fit_3$estimate[1], sdlog = fit_3$estimate[2])

# Scenario 4
scenario_4_data <- wa %>% filter(Scenario == 4)
fit_4 <- fitdist(scenario_4_data$WaterAge, "lnorm", method = "mle")
waterage_4 <- rlnorm(n, meanlog = fit_4$estimate[1], sdlog = fit_4$estimate[2])

# Scenario 5
scenario_5_data <- wa %>% filter(Scenario == 5)
fit_5 <- fitdist(scenario_5_data$WaterAge, "lnorm", method = "mle")
waterage_5 <- rlnorm(n, meanlog = fit_5$estimate[1], sdlog = fit_5$estimate[2])

# Scenario 6
scenario_6_data <- wa %>% filter(Scenario == 6)
fit_6 <- fitdist(scenario_6_data$WaterAge, "lnorm", method = "mle")
waterage_6 <- rlnorm(n, meanlog = fit_6$estimate[1], sdlog = fit_6$estimate[2])

# Scenario 7
scenario_7_data <- wa %>% filter(Scenario == 7)
fit_7 <- fitdist(scenario_7_data$WaterAge, "lnorm", method = "mle")
waterage_7 <- rlnorm(n, meanlog = fit_7$estimate[1], sdlog = fit_7$estimate[2])

# Scenario 8
scenario_8_data <- wa %>% filter(Scenario == 8)
fit_8 <- fitdist(scenario_8_data$WaterAge, "lnorm", method = "mle")
waterage_8 <- rlnorm(n, meanlog = fit_8$estimate[1], sdlog = fit_8$estimate[2])

#Starting concentrations set to 0.01 to 100 gc/L
C0_s<-runif(n,0.01,100)

## INFECTION SCENARIOS------------
pinf_sens1<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_1, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))

cor(C0_s,pinf_sens1,method="spearman")
pinf_sens2<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_2, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pinf_sens3<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_3, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pinf_sens4<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_4, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pinf_sens5<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_5, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pinf_sens6<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_6, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pinf_sens7<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_7, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pinf_sens8<-P(kinf,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_8, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))

##CSI scenarios----------
pcsi_sens1<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_1, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))

pcsi_sens2<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_2, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pcsi_sens3<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_3, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pcsi_sens4<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_4, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pcsi_sens5<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_5, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pcsi_sens6<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_6, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pcsi_sens7<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_7, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))
pcsi_sens8<-P(kcsi,mapply(function(x, y) {
  gciu((as.numeric((ratio*(grow_baranyi(x, c(y0=y, mumax=mumax_hot, K=Nmax_hot,h0=0))[,2])+              
                      (1-ratio)*(grow_baranyi(x, c(y0=C0, mumax=mumax_cold, K=Nmax_cold,h0=0))[,2])))))
},
waterage_8, C0_s)*(B * tsh_start * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)))+
  C0_s*(B * (tsh_final-tsh_start) * (((Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)) * 1e3) * ((F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10))))




#INFECTION
sens1<-data.frame(waterage_1,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens1)
sens2<-data.frame(waterage_2,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens2)
sens3<-data.frame(waterage_3,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens3)
sens4<-data.frame(waterage_4,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens4)
sens5<-data.frame(waterage_5,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens5)
sens6<-data.frame(waterage_6,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens6)
sens7<-data.frame(waterage_7,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens7)
sens8<-data.frame(waterage_8,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kinf,pinf_sens8)


correlation_results_inf <- data.frame(
  Parameter = c(
    "waterage", "B", "tsh", "Caer1", "Caer2", "Caer3", "Caer4", "Caer5",
    "Caer6", "Caer7", "Caer8", "Caer9", "Caer10", "D1", "D2", "D3", "D4",
    "D5", "D6", "D7", "D8", "D9", "D10", "C0", "kinf"
  ),
  stringsAsFactors = FALSE
)

correlations_inf<- for (i in 1:8) {
  name<-paste0("sens",i)
  df<-get(name)
  correlations <- cor(df, method = "spearman")
  coefficients <- data.frame(Parameter = row.names(correlations), Spearman_rho = correlations[, 26])
  colnames(coefficients) <- c("Parameter", name)
  coefficients <- coefficients[1:25, ]
  correlation_results_inf[[name]] <- coefficients[2]
}

correlation_results_inf
correlation_results_inf<-as.data.frame(t(correlation_results_inf))
colnames(correlation_results_inf)<-correlation_results_inf[1,]
correlation_results_inf<-correlation_results_inf[2:9,]
correlation_results_inf$row<-row.names(correlation_results_inf)
correlation_results_inf<-reshape2::melt(correlation_results_inf,id.vars = 'row')
correlation_results_inf$value<-as.numeric(correlation_results_inf$value)

##Result
View(correlation_results_inf)


#ILLNESS

sens1<-data.frame(waterage_1,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens1)
sens2<-data.frame(waterage_2,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens2)
sens3<-data.frame(waterage_3,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens3)
sens4<-data.frame(waterage_4,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens4)
sens5<-data.frame(waterage_5,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens5)
sens6<-data.frame(waterage_6,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens6)
sens7<-data.frame(waterage_7,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens7)
sens8<-data.frame(waterage_8,B,tsh_final,Caer1,Caer2,Caer3,Caer4,Caer5,Caer6,Caer7,Caer8,Caer9,Caer10,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10, C0_s,kcsi,pcsi_sens8)



correlation_results_csi <- data.frame(
  Parameter = c(
    "waterage", "B", "tsh", "Caer1", "Caer2", "Caer3", "Caer4", "Caer5",
    "Caer6", "Caer7", "Caer8", "Caer9", "Caer10", "D1", "D2", "D3", "D4",
    "D5", "D6", "D7", "D8", "D9", "D10", "C0", "kcsi"
  ),
  stringsAsFactors = FALSE
)

correlations_csi<- for (i in 1:8) {
  name<-paste0("sens",i)
  df<-get(name)
  correlations <- cor(df, method = "spearman")
  coefficients <- data.frame(Parameter = row.names(correlations), Spearman_rho = correlations[, 26])
  colnames(coefficients) <- c("Parameter", name)
  coefficients <- coefficients[1:25, ]
  correlation_results_csi[[name]] <- coefficients[2]
}

correlation_results_csi
correlation_results_csi<-as.data.frame(t(correlation_results_csi))
colnames(correlation_results_csi)<-correlation_results_csi[1,]
correlation_results_csi<-correlation_results_csi[2:9,]
correlation_results_csi$row<-row.names(correlation_results_csi)
correlation_results_csi<-reshape2::melt(correlation_results_csi,id.vars = 'row')
correlation_results_csi$value<-as.numeric(correlation_results_csi$value)

##Result
View(correlation_results_csi)
