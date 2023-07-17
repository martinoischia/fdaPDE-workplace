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
func_evaluation1 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(0.5, nlocs))
func_evaluation2 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(1, nlocs))
func_evaluation3 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(1.5, nlocs))

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
# mod1_info <- smooth.FEM.mixed(
#     locations = loc,
#     observations = mixed_obs2, # obs doesn't effect tree, bary, dof
#     covariates = covariates,
#     FEMbasis = FEMbasis,
#     random_effect = c(1),
#     lambda = lambda,
#     lambda.selection.lossfunction = 'GCV',
#     DOF.stochastic.realizations =30,
#     DOF.evaluation = GCVMETHODFLAG,
#     FLAG_ITERATIVE = FALSE)
	 

# load("D:/VM/Tesi/my/2D/Cshape/kim.RData")
# # print("vettore f")
# # mod2_info$fit.FEM.mixed$coeff [1:10]
# # print ("differenza vettore f tra i due solver")
# # (mod2_info$fit.FEM.mixed$coeff - mod1_info$fit.FEM.mixed$coeff)[1:10]
# # print ("differenza con kim")
# # (kim$mod2_info$fit.FEM.mixed$coeff - mod1_info$fit.FEM.mixed$coeff)[1:10]

# # print ("dof iterative: ")
# # mod2_info$edf
# # print ("kim ")
# # kim$mod2_info$edf

# # print("beta iterativo")
# # mod2_info$beta
# # print("beta normale")
# # mod1_info$beta
# # print("best lambda normale")
# # print(mod2_info$bestlambda)
# # print("best lambda iterativo")
# # print(mod1_info$bestlambda)

time.taken <- 0
f_selected_lambda <- rep(0, N)
ran <- range(c(func_evaluation1, func_evaluation2, func_evaluation3))
for (i in 1:N) {
    data2 <- mixed_obs2 + rnorm(nlocs * 3, mean = 0, sd = 0.05 * (ran[2] - ran[1]))

    # measure time
    start.time <- Sys.time()
    mod2 <- smooth.FEM.mixed(
        locations = loc,
        observations = data2,
        covariates = covariates,
        FEMbasis = mod2_info$fit.FEM$FEMbasis, # reuse tree info
        lambda = lambda,
        lambda.selection.lossfunction = 'GCV',
        bary.locations = mod2_info$bary.locations, # reuse bary info
        DOF.matrix = mod2_info$edf,
		FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol,
		threshold_residual = tol,
        verbose=TRUE
    ) # reuse dof info
    end.time <- Sys.time()
    time.taken <- (end.time - start.time) + time.taken

    iterations <- c(iterations, mod2$iterations[which.min(mod2$GCV)])	 
    residual <- c(residual, mod2$residuals[[which.min(mod2$GCV)]][length(mod2$residuals[[which.min(mod2$GCV)]])])	 
    f_selected_lambda[i] <- lambda[which.min(mod2$GCV)] # returns index

    # # beta estimates
    f_betamat[1, i] <- mod2$beta[, which.min(mod2$GCV)][1]
    # f_betamat[, ,i] <- mod2$beta
    # f_betamat[, ,i] <- mod2$beta

    # beta MSE
    RMSE <- RMSE_func(beta_exact[1], f_betamat[1, i])
    f_rmse_beta1 <- c(f_rmse_beta1, RMSE)

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
        W %*% beta_exact + func_evaluation1,
        W %*% f_betamat[, i] +  fitted.function2[1:nlocs]
    )
    f_global_rmse1 <- c(f_global_rmse1, RMSE)

    RMSE <- RMSE_func(
        W %*% beta_exact + func_evaluation2,
        W %*% f_betamat[, i] + fitted.function2[(nlocs + 1):(2 * nlocs)]
    )
    f_global_rmse2 <- c(f_global_rmse2, RMSE)

    RMSE <- RMSE_func(
        W %*% beta_exact + func_evaluation3,
        W %*% f_betamat[, i] +  fitted.function2[(2 * nlocs + 1):(3 * nlocs)]
    )
    f_global_rmse3 <- c(f_global_rmse3, RMSE)
}
half_range=0.15
x11()
boxplot(f_betamat[1,], ylim=c(beta_exact-half_range ,beta_exact+0.10 ), cex.main=1.25, cex.axis=1.25, cex.sub=1.25)
abline(h=beta_exact, col='gray60')

DF <- data.frame(f_global_rmse1,
                 f_global_rmse2,
                 f_global_rmse3)
x11()
boxplot(DF,  xaxt = "n", main='Global RMSE')
axis(side = 1,  labels = c("unit1","unit2","unit3"), at=1:3)