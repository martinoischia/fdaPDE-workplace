rm(list = ls())
set.seed(5847947)
library(latex2exp)

load("D:/connectome/all_patients/connectivityMaps.RData")
num_unit= ncol(observations)
library(fdaPDE)

GCVFLAG <- TRUE
GCVMETHODFLAG <- 'stochastic' 
# GCVMETHODFLAG <-"exact"
tol <- 1e-6



# (1): full mesh, full locs
wanted_unit = 59
new_obs = observations[,1:wanted_unit]
new_cov = covariates[1:(dim(observations)[1]*wanted_unit)]

length(new_cov)

lambda=seq(10^-5, 10^-3, length.out=10) #4 patients
lambda=seq(10^-4, 10^-2, length.out=20) #8 patients
lambda=seq(0.001, 0.005, length.out=20) #12 patients
lambda=seq(0.001, 0.005, length.out=8) #12 patients
output_CPP = smooth.FEM.mixed(observations = new_obs,
                              covariates = new_cov,
                              FEMbasis = FEMbasis,
                              lambda = lambda,
							  lambda.selection.lossfunction = 'GCV',
							  DOF.evaluation = GCVMETHODFLAG,
                              DOF.stochastic.realizations =60,
							  verbose=TRUE,
							  FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol, threshold_residual = tol)

#12 patients
output_CPP$bestlambda #8
lambda[output_CPP$bestlambda] # 0.002473684
output_CPP$beta[,output_CPP$bestlambda]
save.image("fixed_59_patients_imgenius")
# 0.006514516
#         b_11          b_21          b_31          b_41          b_51          b_61          b_71          b_81 
# 0.0066936700  0.0046181333 -0.0031220337  0.0005023062 -0.0018819550  0.0017680683  0.0002864598 -0.0057483129 
#         b_91         b_101         b_111         b_121 
# 0.0088133374 -0.0055524296 -0.0052526283 -0.0011246155 


# plot(FEM.mixed(output_CPP$fit.FEM.mixed$coeff[,output_CPP$bestlambda], 
#                num_units=wanted_unit,
#                FEMbasis))

# plot(FEM.mixed(as.vector(observations[,1:12]), 
#                num_units=wanted_unit,
#                FEMbasis)) #black and white...
# length(which(is.na(as.vector(observations[,1:12])))) #33552 so black and white...
# length(which(is.na(as.vector(output_CPP$fit.FEM.mixed$coeff[,output_CPP$bestlambda])))) #0



# #### plot test #####
# i=1
# which.min(output_CPP$fit.FEM.mixed$coeff[((i-1)*32492+1):(i*32492),output_CPP$bestlambda]) 
# #-0.8897433, 20216
# which.max(output_CPP$fit.FEM.mixed$coeff[((i-1)*32492+1):(i*32492),output_CPP$bestlambda]) 
# #1.002585, 24510
# #min should be -1, max should be 1

# plot(FEM(output_CPP$fit.FEM.mixed$coeff[((i-1)*32492+1):(i*32492),output_CPP$bestlambda], 
#                FEMbasis))
# pch3d(mesh$nodes[(20216-1):(20216+1),], col="pink", cex=10, pch=19) #located in blue area (min)
# pch3d(t(mesh$nodes[20216,]), col="pink", size=1000)

# plot(FEM(output_CPP$fit.FEM.mixed$coeff[((i-1)*32492+1):(i*32492),output_CPP$bestlambda], 
#          FEMbasis))
# pch3d(mesh$nodes[(24510-1):(24510+1),], col="pink", cex=10, pch=19) #located in red area (max)


# i=10
# min(output_CPP$fit.FEM.mixed$coeff[((i-1)*10000+1):(i*10000),output_CPP$bestlambda]) 
# #-0.2796454
# max(output_CPP$fit.FEM.mixed$coeff[((i-1)*10000+1):(10000),output_CPP$bestlambda]) 
# #1.033175
# #min should be -1, max should be 1
