# Author(s): Dean Hansen

library('viridis')
source('Main/Scripts/Crosstabulate.R')
source('Main/Scripts/Sample - Indexes.R')
source('Main/Scripts/Indexes - Maximum Likelihood Estimation.R')
posterior_ratemaking_2017 <- read.csv('Main/Data/Bermudez and Karlis_2017/Posterior_Ratemaking_2017_Rows.csv')
posterior_ratemaking_2017_tab <- read.csv('Main/Data/Bermudez and Karlis_2017/Posterior_Ratemaking_2017_Table.csv')
N <- nrow(posterior_ratemaking_2017)

# MLEs for Simulated Data -------------------------------------------------------
pr_pois <- bivpois::bp.mle(posterior_ratemaking_2017)
pr_pois_mle <- as.vector(pr_pois$lambda)
round(pr_pois_mle,4)
mle_indices(pr_pois_mle, dist = "POIS")

# [1] 0.0670 0.0884 0.0140
# -logL -53271.05
# [1] 1.0000 1.0136 1.0000

data <- posterior_ratemaking_2017_tab
pr_PA <- NR.F(MofM(data),data)
pr_PA_mle <- as.vector(pr_PA[,1])
round(pr_PA_mle,4)
evalLL(pr_PA_mle,data)
mle_indices(pr_PA_mle, dist = "PA")

# [1] 0.0516 0.0658 0.0132 0.2167
# -logL -48087.48
# [1] 1.5533 1.7849 1.5533

data <- posterior_ratemaking_2017
pr_nbinom <- bzinb::bnb(data[,1],data[,2],em=T,maxiter=10000,tol=1e-5,showFlag=T)
pr_nbinom_mle <- as.vector(coef(pr_nbinom)[,1])
round(pr_nbinom_mle,4)
mle_indices(pr_nbinom_mle, dist = "NB")

# [1] 0.0822 0.0684 0.0710 0.5372 0.6695
# -logL -48050.41
# [1] 1.6034 1.7991 1.6044
