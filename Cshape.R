rm(list = ls())
set.seed(5847947)
# Loading the output of Kim's fdaPDE for checking consistency of results
# (or you can install that library and run this same script modifying
# properly the lines where you find a triple pound [###])
kim= new.env() ###
load("D:/VM/Tesi/my/kim_full.RData",envir=kim) ###

library(fdaPDE)

# memory profiling, temporarily disabled since there is other stuff to fix
# Rprof("Cshape.out", memory.profiling=TRUE, interval = 0.01)

data(horseshoe2D)
attach(horseshoe2D)
# Create mesh and FE basis
mesh <- create.mesh.2D(nodes = boundary_nodes, segments = boundary_segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

lambda <- 10^seq(-2, 1, by = 0.1)
GCVFLAG <- TRUE
GCVMETHODFLAG <- "exact"

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

# Observations
# range
ran <- range(c(func_evaluation1, func_evaluation2, func_evaluation3))
# Set number of simulation trials
N <- 50
# by column
mixed_obs2 <- cbind(
    W %*% beta_exact1 + func_evaluation1,
    W %*% beta_exact2 + func_evaluation2,
    W %*% beta_exact3 + func_evaluation3
)

covariates <- rbind(W, W, W)

# For postprocessing
f_betamat <- matrix(data = NA, nrow = 2, ncol = N)
f_bimat <- matrix(data = NA, nrow = 3, ncol = N)

RMSE_func <- function(v, w) {
    sqrt(mean((v - w)^2))
}

f_rmse_beta1 <- NULL
f_rmse_beta2 <- NULL
f_rmse_b1_1 <- NULL
f_rmse_b2_1 <- NULL
f_rmse_b3_1 <- NULL
f_rmse_f1 <- NULL # summed on the locations (just the f)
f_rmse_f2 <- NULL # summed on the locations (just the f)
f_rmse_f3 <- NULL # summed on the locations (just the f)
f_global_rmse1 <- NULL # error on unit1
f_global_rmse2 <- NULL # error on unit2
f_global_rmse3 <- NULL # error on unit3

# Just for efficiency?

mod2_info <- smooth.FEM.mixed(
    locations = loc,
    observations = mixed_obs2, # obs doesn't effect tree, bary, dof
    covariates = covariates,
    FEMbasis = FEMbasis,
    random_effect = c(1),
    lambda = lambda,
    lambda.selection.lossfunction = 'GCV', ###
    #GCV = GCVFLAG, ###
    DOF.evaluation = GCVMETHODFLAG) ###
    #GCVmethod = 'Exact') ###

# plot(mod2_info$fit.FEM.time, 1) or mixed ?
# plot(mod2_info$fit.FEM.time, 2)
# plot(mod2_info$fit.FEM.time, 3)

all.equal(kim$mod2_info$fit.FEM.mixed$coeff, mod2_info$fit.FEM.mixed$coeff) ###
kim$mod2_info$beta - mod2_info$beta ###

time.taken <- 0
f_selected_lambda <- rep(0, N)
for (i in 1:N) {
    data2 <- mixed_obs2 + rnorm(nlocs * 3, mean = 0, sd = 0.05 * (ran[2] - ran[1]))

    # measure time
    start.time <- Sys.time()
    mod2 <- smooth.FEM.mixed(
        locations = loc,
        observations = data2,
        covariates = covariates,
        FEMbasis = mod2_info$fit.FEM$FEMbasis, # reuse tree info
        random_effect = c(1),
        lambda = lambda,
        lambda.selection.lossfunction = 'GCV', ###
        #GCV = GCVFLAG, ###
        bary.locations = mod2_info$bary.locations, # reuse bary info
        DOF.matrix = mod2_info$edf ###
        #DOF_matrix = mod2_info$edf ###
    ) # reuse dof info
    end.time <- Sys.time()
    time.taken <- (end.time - start.time) + time.taken

    f_selected_lambda[i] <- which.min(mod2$GCV) # returns index

    # beta estimates
    f_betamat[1, i] <- mod2$beta[, which.min(mod2$GCV)][1]
    f_betamat[2, i] <- mod2$beta[, which.min(mod2$GCV)][2]

    # beta MSE
    RMSE <- RMSE_func(beta_exact2[1], f_betamat[1, i])
    f_rmse_beta1 <- c(f_rmse_beta1, RMSE)

    RMSE <- RMSE_func(beta_exact2[2], f_betamat[2, i])
    f_rmse_beta2 <- c(f_rmse_beta2, RMSE)

    # bi estimates
    f_bimat[1, i] <- mod2$b_i[, which.min(mod2$GCV)][1]
    f_bimat[2, i] <- mod2$b_i[, which.min(mod2$GCV)][2]
    f_bimat[3, i] <- mod2$b_i[, which.min(mod2$GCV)][3]

    # bi MSE
    RMSE <- RMSE_func(-5, f_bimat[1, i])
    f_rmse_b1_1 <- c(f_rmse_b1_1, RMSE)

    RMSE <- RMSE_func(0, f_bimat[2, i])
    f_rmse_b2_1 <- c(f_rmse_b2_1, RMSE)

    RMSE <- RMSE_func(5, f_bimat[3, i])
    f_rmse_b3_1 <- c(f_rmse_b3_1, RMSE)

    # f RMSE
    fitted.function2 <- eval.FEM.mixed(mod2$fit.FEM,
        locations = loc
    )[, which.min(mod2$GCV)]

    RMSE <- RMSE_func(func_evaluation1, fitted.function2[1:nlocs])
    f_rmse_f1 <- c(f_rmse_f1, RMSE)

    RMSE <- RMSE_func(func_evaluation2, fitted.function2[(nlocs + 1):(2 * nlocs)])
    f_rmse_f2 <- c(f_rmse_f2, RMSE)

    RMSE <- RMSE_func(func_evaluation3, fitted.function2[(2 * nlocs + 1):(3 * nlocs)])
    f_rmse_f3 <- c(f_rmse_f3, RMSE)


    # global RMSE
    RMSE <- RMSE_func(
        W %*% beta_exact1 + func_evaluation1,
        W %*% f_betamat[, i] + cov1 * f_bimat[1, i] + fitted.function2[1:nlocs]
    )
    f_global_rmse1 <- c(f_global_rmse1, RMSE)

    RMSE <- RMSE_func(
        W %*% beta_exact2 + func_evaluation2,
        W %*% f_betamat[, i] + cov1 * f_bimat[2, i] + fitted.function2[(nlocs + 1):(2 * nlocs)]
    )
    f_global_rmse2 <- c(f_global_rmse2, RMSE)

    RMSE <- RMSE_func(
        W %*% beta_exact3 + func_evaluation3,
        W %*% f_betamat[, i] + cov1 * f_bimat[3, i] + fitted.function2[(2 * nlocs + 1):(3 * nlocs)]
    )
    f_global_rmse3 <- c(f_global_rmse3, RMSE)
}

# Rprof(NULL)
# summaryRprof("Cshape.out", memory = "tseries", diff = TRUE)
table(f_selected_lambda)
detach(horseshoe2D)
