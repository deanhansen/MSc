# Author(s): Dean Hansen

source("Main/Scripts/MLE/Method of Moments.R")
source("Main/Scripts/MLE/Derivatives/DL1.R")
source("Main/Scripts/MLE/Derivatives/DL2.R")
source("Main/Scripts/MLE/Derivatives/DL3.R")
source("Main/Scripts/MLE/Derivatives/DRho.R")
source("Main/Scripts/MLE/Derivatives/DDL1.R")
source("Main/Scripts/MLE/Derivatives/DDL2.R")
source("Main/Scripts/MLE/Derivatives/DDL3.R")
source("Main/Scripts/MLE/Derivatives/DDRho.R")
source("Main/Scripts/MLE/Derivatives/DD_L12.R")
source("Main/Scripts/MLE/Derivatives/DD_L13.R")
source("Main/Scripts/MLE/Derivatives/DD_L1R.R")
source("Main/Scripts/MLE/Derivatives/DD_L23.R")
source("Main/Scripts/MLE/Derivatives/DD_L2R.R")
source("Main/Scripts/MLE/Derivatives/DD_L3R.R")
source("Main/Scripts/MLE/NR Algorithm.R")
source("Main/Scripts/MLE/Gradient.R")
source("Main/Scripts/MLE/Hessian.R")
source('Main/Scripts/Helpers - Probability Mass Functions.R')
source('Main/Scripts/Helpers - Maximum Likelihood Estimation Indices.R')
source("Main/Scripts/Helpers - Delta Method.R")
source('Main/Scripts/Helpers - Simulation Study.R')
set.seed(1)
N <- 500
B <- 10000

# Data ----------------------------------------------------------------
d1 <- bivpois::rbp(N, c(2,4,1))
d2 <- bzinb::rbnb(N,4,4,6,0.25,0.25)
d3 <- rbivpa(c(0.6,0.6,0.9,0.25),N)
statistics(d1)
statistics(d2)
statistics(d3)

# Case #1 -----------------------------------------------------------------
case_1_pois <- mle(d1,dist="POIS")
case_1_pois_mle <- as.vector(case_1_pois$lambda)
round(case_1_pois_mle,4)
round(case_1_pois$loglik[1],2)
mle_indices(case_1_pois_mle,dist="POIS")


d1_tab <- crosstabulate(d1)
case_1_PA <- mle(d1_tab,params=c(2,2.5,1,0.25),dist="PA")
case_1_PA_mle <- as.vector(case_1_PA[,1])
round(case_1_PA_mle,4)
round(evalLL(case_1_PA_mle,d1_tab),2)
mle_indices(case_1_PA_mle,dist="PA")


case_1_nbinom <- mle(d1,dist="NB")
case_1_nbinom_mle <- as.vector(coef(case_1_nbinom)[,1])
round(case_1_nbinom_mle,4)
round(case_1_nbinom$lik,2)
mle_indices(case_1_nbinom_mle,dist="NB")


# Case #2 -----------------------------------------------------------------
case_2_pois <- mle(d2,dist="POIS")
case_2_pois_mle <- as.vector(case_2_pois$lambda)
round(case_2_pois_mle,4)
round(case_2_pois$loglik[1],2)
mle_indices(case_2_pois_mle,dist="POIS")


d2_tab <- crosstabulate(d2)
case_2_PA <- mle(d2_tab,params=c(2,2.5,1,0.25),dist="PA")
case_2_PA_mle <- as.vector(case_2_PA[,1])
round(case_2_PA_mle,4)
round(evalLL(case_2_PA_mle,d2_tab),2)
mle_indices(case_2_PA_mle,dist="PA")


case_2_nbinom <- mle(d2,dist="NB")
case_2_nbinom_mle <- as.vector(coef(case_2_nbinom)[,1])
round(case_2_nbinom_mle,4)
round(case_2_nbinom$lik,2)
mle_indices(case_2_nbinom_mle,dist="NB")


# Case #3 -----------------------------------------------------------------
case_3_pois <- mle(d3,dist="POIS")
case_3_pois_mle <- as.vector(case_3_pois$lambda)
round(case_3_pois_mle,4)
round(case_3_pois$loglik[1],2)
mle_indices(case_3_pois_mle,dist="POIS")


d3_tab <- crosstabulate(d3)
case_3_PA <- mle(d3_tab,dist="PA")
case_3_PA_mle <- as.vector(case_3_PA[,1])
round(case_3_PA_mle,4)
round(evalLL(case_3_PA_mle,d3_tab),2)
mle_indices(case_3_PA_mle,dist="PA")


case_3_nbinom <- mle(d3,dist="NB")
case_3_nbinom_mle <- as.vector(coef(case_3_nbinom)[,1])
round(case_3_nbinom_mle,4)
round(case_3_nbinom$lik,2)
mle_indices(case_3_nbinom_mle,dist="NB")


# Simulation Study --------------------------------------------------------
simulation_study(B,50,d1)
simulation_study(B,50,d2)
simulation_study(B,50,d3)
