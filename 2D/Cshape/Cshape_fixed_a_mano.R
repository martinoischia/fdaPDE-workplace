rm(list = ls())
set.seed(5847947)
library(latex2exp)
library(fdaPDE)

# memory profiling, temporarily disabled since there is other stuff to fix
# Rprof("Cshape.out", memory.profiling=TRUE, interval = 0.01)

data(horseshoe2D)
attach(horseshoe2D)
# Create mesh and FE basis
mesh <- create.mesh.2D(nodes = boundary_nodes, segments = boundary_segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

lambda = 10^(-0.5)
GCVFLAG <- TRUE
GCVMETHODFLAG <- 'stochastic' 
# GCVMETHODFLAG <-"exact"
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

W <- cbind(cov1)

# Fix betas
beta_exact <- 3

## Fix f_i, for i=1,2,3
func_evaluation1 <- numeric(nlocs)
func_evaluation2 <- numeric(nlocs)
func_evaluation3 <- numeric(nlocs)
func_evaluation1 <- rep(1, nlocs) 
func_evaluation2 <- rep(0, nlocs)
func_evaluation3 <- rep(3, nlocs)

# by column
mixed_obs2 <- cbind(
    W %*% beta_exact + func_evaluation1,
    W %*% beta_exact + func_evaluation2,
    W %*% beta_exact + func_evaluation3
)

covariates <- rbind(W, W, W)

# For postprocessing
# Set number of simulation trials
N <- 50
f_betamat <- matrix(data = NA, nrow = 1, ncol = N)

RMSE_func <- function(v, w) {
    sqrt(mean((v - w)^2))
}

iterations <- NULL
residual <- NULL
f_rmse_beta1 <- NULL
f_rmse_f1 <- NULL # summed on the locations (just the f)
f_rmse_f2 <- NULL # summed on the locations (just the f)
f_rmse_f3 <- NULL # summed on the locations (just the f)
f_global_rmse1 <- NULL # error on unit1
f_global_rmse2 <- NULL # error on unit2
f_global_rmse3 <- NULL # error on unit3

mod2_info <- smooth.FEM.mixed(
    locations = loc,
    observations = mixed_obs2, # obs doesn't effect tree, bary, dof
    covariates = covariates,
    FEMbasis = FEMbasis,
    lambda = lambda,
    lambda.selection.lossfunction = 'GCV',
    DOF.evaluation = GCVMETHODFLAG,
    verbose=TRUE,
    DOF.stochastic.realizations =60,
    FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol, threshold_residual = tol)