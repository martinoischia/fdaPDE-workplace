rm(list = ls())
set.seed(5847947)

library(fdaPDE)

# memory profiling, temporarily disabled since there is other stuff to fix
# Rprof("Cshape.out", memory.profiling=TRUE, interval = 0.01)

data(horseshoe2D)
attach(horseshoe2D)
# Create mesh and FE basis
mesh <- create.mesh.2D(nodes = boundary_nodes, segments = boundary_segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

# lambda = 10^(-2)
lambda <- 10^seq(-2, 1, by = 0.1)
GCVFLAG <- TRUE
# GCVMETHODFLAG <- 'stochastic' 
GCVMETHODFLAG <-"exact"
tol <- 1e-6

# Almost the same as the one in Massardi Spaziani report !?
fs.test.time <- function(x, y, t = y) {
    K <- (y / 0.1 * as.double((abs(y) <= 0.1 & x > -0.5)) + as.double((abs(y) > 0.1 | x <= -0.5)))^2 # ???

    res <- numeric(length = length(x))

    for (i in seq_len(length(x))) {
        if (x[i] >= 0 && y[i] > 0) {
              res[i] <- cos(t[i]) * (0.25 * pi + x[i]) + (y[i] - 0.5)^2
          }

        if (x[i] >= 0 && y[i] <= 0) {
              res[i] <- cos(2 * t[i]) * (-0.25 * pi - x[i]) + (-y[i] - 0.5)^2
          } # ??

        if (x[i] < 0 && y[i] > 0) {
              res[i] <- cos(t[i]) * (-atan(y[i] / x[i]) * 0.5) + (sqrt(x[i]^2 + y[i]^2) - 0.5)^2 * K[i]
          }

        if (x[i] < 0 && y[i] <= 0) {
              res[i] <- cos(2 * t[i]) * (-atan(y[i] / x[i]) * 0.5) + (sqrt(x[i]^2 + y[i]^2) - 0.5)^2 * K[i]
          }
    }

    return(res)
}

# nodes number
nnodes <- dim(mesh$nodes)[1]
# Picking 100 random points on a ~ uniform grid on the C-shape
gridnum <- 20
nlocs <- 100
minV1 <- min(mesh$nodes[, 1])
maxV1 <- max(mesh$nodes[, 1])
minV2 <- min(mesh$nodes[, 2])
maxV2 <- max(mesh$nodes[, 2])
x <- seq(minV1 + 0.1, maxV1 - 0.1, length.out = gridnum)
y <- seq(minV2 + 0.1, maxV2 - 0.1, length.out = gridnum)
unif_grid <- expand.grid(x = x, y = y)
points1 <- eval.FEM(FEM(rep(1, nnodes), FEMbasis), locations = unif_grid)
loc <- unif_grid[-which(is.na(points1)), ]
ind <- sample.int(dim(loc)[1])[1:nlocs]
loc <- loc[ind, ]

# covariates
cov1 <- sin(2 * pi * loc[, 1]) * cos(2 * pi * loc[, 2])

# I think it would be better to use others in unit 2 and 3 but w/e
cov2 <- rnorm(nlocs, mean = 0, sd = 2)
W <- cbind(cov1, cov2)

# Fix betas
beta_exact1 <- c(3 - 5, 0.5)
beta_exact2 <- c(3, 0.5)
beta_exact3 <- c(3 + 5, 0.5)

## Fix f_i, for i=1,2,3
func_evaluation1 <- numeric(nlocs)
func_evaluation2 <- numeric(nlocs)
func_evaluation3 <- numeric(nlocs)
func_evaluation1 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(0.5, nlocs))
func_evaluation2 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(1, nlocs))
func_evaluation3 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(1.5, nlocs))

# by column
mixed_obs2 <- cbind(
    W %*% beta_exact1 + func_evaluation1,
    W %*% beta_exact2 + func_evaluation2
)

mixed_obs1 <- cbind(
    W %*% beta_exact1 + func_evaluation1,
    W %*% beta_exact2 + func_evaluation2,
    W %*% beta_exact3 + func_evaluation3
)
mixed_obs1 [3:nrow(mixed_obs1),3] = NA

covariates1 <- rbind(W, W, W)
covariates2 <- rbind(W,W)
# For postprocessing
# Set number of simulation trials
# N <- 3
# f_betamat <- matrix(data = NA, nrow = 2, ncol = N)
# # f_betamat = array(data=NA, dim = c (2, length(lambda), N))
# f_bimat <- matrix(data = NA, nrow = 3, ncol = N)

# RMSE_func <- function(v, w) {
#     sqrt(mean((v - w)^2))
# }

# f_rmse_beta1 <- NULL
# f_rmse_beta2 <- NULL
# f_rmse_b1_1 <- NULL
# f_rmse_b2_1 <- NULL
# f_rmse_b3_1 <- NULL
# f_rmse_f1 <- NULL # summed on the locations (just the f)
# f_rmse_f2 <- NULL # summed on the locations (just the f)
# f_rmse_f3 <- NULL # summed on the locations (just the f)
# f_global_rmse1 <- NULL # error on unit1
# f_global_rmse2 <- NULL # error on unit2
# f_global_rmse3 <- NULL # error on unit3

# Just for efficiency?

mod2_info <- smooth.FEM.mixed(
    locations = loc,
    observations = mixed_obs2, # obs doesn't effect tree, bary, dof
    covariates = covariates2,
    FEMbasis = FEMbasis,
    random_effect = c(1),
    lambda = lambda,
    lambda.selection.lossfunction = 'GCV',###
    DOF.stochastic.realizations = 100,###
    DOF.evaluation = GCVMETHODFLAG, ###
    FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol, threshold_residual = tol, verbose= TRUE)###
    # GCV = TRUE, ###
    # GCVmethod = 'Exact') ###

# mod1_info_na <- smooth.FEM.mixed(
#     locations = loc,
#     observations = mixed_obs1, # obs doesn't effect tree, bary, dof
#     covariates = covariates1,
#     FEMbasis = FEMbasis,
#     random_effect = c(1),
#     lambda = lambda,
#     lambda.selection.lossfunction = 'GCV',###
#     DOF.stochastic.realizations = 100,###
#     DOF.evaluation = GCVMETHODFLAG, ###
#     FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol, threshold_residual = tol, verbose = TRUE) ###
#     # GCV = TRUE, ###
#     # GCVmethod = 'Exact') ###

mod1_info <- smooth.FEM.mixed(
    locations = loc,
    observations = mixed_obs1, # obs doesn't effect tree, bary, dof
    covariates = covariates1,
    FEMbasis = FEMbasis,
    random_effect = c(1),
    lambda = lambda,
    lambda.selection.lossfunction = 'GCV',
    DOF.evaluation = GCVMETHODFLAG,
    FLAG_ITERATIVE = FALSE)

mod2_info$beta
# mod1_info_na$beta
mod1_info$beta
detach(horseshoe2D)