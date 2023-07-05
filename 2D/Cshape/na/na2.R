rm(list = ls())
set.seed(5847947)
library(latex2exp)
library(fdaPDE)
load("D:/VM/Tesi/my/2D/Cshape/kim.RData")

# memory profiling, temporarily disabled since there is other stuff to fix
# Rprof("Cshape.out", memory.profiling=TRUE, interval = 0.01)
setwd("d:\\VM\\Tesi\\my\\2D\\Cshape\\na")
data(horseshoe2D)
attach(horseshoe2D)
# Create mesh and FE basis
mesh <- create.mesh.2D(nodes = boundary_nodes, segments = boundary_segments)
mesh <- refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

# lambda = 10^(-2)
loglambda=seq(-2, 1, by = 0.1) 
lambda <- 10^loglambda
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
nlocs <- 60
nmore = 100
minV1 <- min(mesh$nodes[, 1])
maxV1 <- max(mesh$nodes[, 1])
minV2 <- min(mesh$nodes[, 2])
maxV2 <- max(mesh$nodes[, 2])
x <- seq(minV1 + 0.1, maxV1 - 0.1, length.out = gridnum)
y <- seq(minV2 + 0.1, maxV2 - 0.1, length.out = gridnum)
unif_grid <- expand.grid(x = x, y = y)
points1 <- eval.FEM(FEM(rep(1, nnodes), FEMbasis), locations = unif_grid)
loc <- unif_grid[-which(is.na(points1)), ]
ind <- sample.int(dim(loc)[1])[1:nmore]
loc <- loc[ind, ]
loc_less = loc[1:nlocs,]

# covariates
cov1 <- sin(2 * pi * loc[, 1]) * cos(2 * pi * loc[, 2])

# I think it would be better to use others in unit 2 and 3 but w/e
cov2 <- rnorm(nmore, mean = 0, sd = 2)
W <- cbind(cov1, cov2)

# Fix betas
beta_exact1 <- c(3 - 5, 0.5)
beta_exact2 <- c(3, 0.5)
beta_exact3 <- c(3 + 5, 0.5)

## Fix f_i, for i=1,2,3
func_evaluation1 <- numeric(nmore)
func_evaluation2 <- numeric(nmore)
func_evaluation3 <- numeric(nmore)
func_evaluation1 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(0.5, nmore))
func_evaluation2 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(1, nmore))
func_evaluation3 <- fs.test.time(x = loc[, 1], y = loc[, 2], t = rep(1.5, nmore))

# by column
mixed_obs <- cbind(
    W %*% beta_exact1 + func_evaluation1,
    W %*% beta_exact2 + func_evaluation2,
    W %*% beta_exact3 + func_evaluation3
)

mixed_obs1 <- mixed_obs

mixed_obs1 [(nlocs+1):nmore , ] = NA
mixed_obs2 = mixed_obs[1:nlocs,]

covariates <- rbind(W, W, W)
Wless = W[1:nlocs,]
covariates2 = rbind(Wless,Wless,Wless)
# For postprocessing
# Set number of simulation trials
N <- 50
f_betamat <- matrix(data = NA, nrow = 2, ncol = N)
# f_betamat = array(data=NA, dim = c (2, length(lambda), N))
f_bimat <- matrix(data = NA, nrow = 3, ncol = N)

RMSE_func <- function(v, w) {
    sqrt(mean((v - w)^2))
}

iterations <- NULL
residual <- NULL
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

# mod2_info <- smooth.FEM.mixed(
#     locations = loc_less,
#     observations = mixed_obs2, # obs doesn't effect tree, bary, dof
#     covariates = covariates2,
#     FEMbasis = FEMbasis,
#     random_effect = c(1),
#     lambda = lambda,
#                 lambda.selection.lossfunction = 'GCV',
#             DOF.evaluation = GCVMETHODFLAG,
#             verbose=TRUE,
#             DOF.stochastic.realizations =60,
#             FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol, threshold_residual = tol)

mod1_info_na <- smooth.FEM.mixed(
    locations = loc,
    observations = mixed_obs1, # obs doesn't effect tree, bary, dof
    covariates = covariates,
    FEMbasis = FEMbasis,
    random_effect = c(1),
    lambda = lambda,
    lambda.selection.lossfunction = 'GCV',###
    DOF.stochastic.realizations = 60,###
    DOF.evaluation = GCVMETHODFLAG, ###
    FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol, threshold_residual = tol, verbose = TRUE) ###


time.taken <- 0
f_selected_lambda <- rep(0, N)
ran <- range(c(func_evaluation1, func_evaluation2, func_evaluation3))
for (i in 1:N) {
    data2 <- mixed_obs1 + rnorm(nlocs * 3, mean = 0, sd = 0.05 * (ran[2] - ran[1]))

    # measure time
    start.time <- Sys.time()
    mod2 <- smooth.FEM.mixed(
        locations = loc,
        observations = data2,
        covariates = covariates,
        FEMbasis = mod1_info_na$fit.FEM$FEMbasis, # reuse tree info
        random_effect = c(1),
        lambda = lambda,
        lambda.selection.lossfunction = 'GCV',
        bary.locations = mod1_info_na$bary.locations, # reuse bary info
        DOF.matrix = mod1_info_na$edf,
		FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol,
		threshold_residual = tol,
        verbose=TRUE
    ) # reuse dof info
    end.time <- Sys.time()
    time.taken <- (end.time - start.time) + time.taken

    iterations <- c(iterations, mod2$iterations[which.min(mod2$GCV)])	 
    residual <- c(residual, mod2$residuals[[which.min(mod2$GCV)]][length(mod2$residuals[[which.min(mod2$GCV)]])])	 
    f_selected_lambda[i] <- loglambda[which.min(mod2$GCV)] # returns index

    # # beta estimates
    f_betamat[1, i] <- mod2$beta[, which.min(mod2$GCV)][1]
    f_betamat[2, i] <- mod2$beta[, which.min(mod2$GCV)][2]
    # f_betamat[, ,i] <- mod2$beta
    # f_betamat[, ,i] <- mod2$beta

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
kim$f_selected_lambda=loglambda[kim$f_selected_lambda]
m=matrix(ncol=length(loglambda),nrow=2)
rownames(m)=c("Monolithic", "Iterative")
colnames(m)=signif(loglambda, 2)
for (i in names(table(signif(kim$f_selected_lambda, 2)))){
    m[1,i]=table(signif(kim$f_selected_lambda,2))[i]
}
for (i in names(table(signif(f_selected_lambda, 2)))){
    m[2,i]=table(signif(f_selected_lambda,2))[i]
}
barplot(m, beside = TRUE, legend = rownames(m), xlab=TeX("$\\log(\\lambda)$"), ylab="counts")


png('beta.png')
par(mfrow = c(1,2))
half_range=0.15
boxplot(cbind(f_betamat[1,],kim$f_betamat[1,]), main=TeX("$\\beta_1$"), names=c("I", "M"), ylim=c(beta_exact2[1]-half_range ,beta_exact2[1]+0.10 ), cex.main=1.25, cex.axis=1.25, cex.sub=1.25)
abline(h=beta_exact2[1], col='gray60')
boxplot(cbind(f_betamat[2,],kim$f_betamat[2,]), main=TeX("$\\beta_2$"), names=c("I", "M"), ylim=c(beta_exact2[2]-half_range,beta_exact2[2]+0.10), cex.main=1.25, cex.axis=1.25, cex.sub=1.25)
abline(h=beta_exact2[2], col='gray60')
dev.off()

png('b.png')
par(mfrow = c(1,3))
boxplot(cbind(kim$f_bimat[1,],f_bimat[1,]), main=TeX("$b_1$"), names=c("I", "M"),ylim=c(-5 -half_range ,-5+half_range ), cex.main=2, cex.axis=2)
abline(h=-5, col='gray60')
boxplot(cbind(kim$f_bimat[2,],f_bimat[2,]),main=TeX("$b_2$"), names=c("I", "M"), ylim=c(-half_range, +half_range),cex.main=2, cex.axis=2)
abline(h=0, col='gray60')
# boxplot(f_global_rmse1)
boxplot(cbind(kim$f_bimat[3,],f_bimat[3,]), main=TeX("$b_3$"), names=c("I", "M"),ylim=c(5 -half_range ,5+half_range ),cex.main=2, cex.axis=2)
abline(h=5, col='gray60')
dev.off()

png('iterations.png')
plot(table(iterations), ylab="counts", xlab="number of iterations")
dev.off()

png('res.png')
boxplot(residual, main="Normalized residual", log="y")
dev.off()

DF <- data.frame(f_global_rmse1, kim$f_global_rmse1, 
                 f_global_rmse2, kim$f_global_rmse2, 
                 f_global_rmse3, kim$f_global_rmse3)
png('rmse.png')
boxplot(DF, col = c('gray60','gray90'), at = c(1:2,4:5,7:8), xaxt = "n", main='Global RMSE')
axis(side = 1, at = c(1.5,4.5,7.5), labels = c("unit1","unit2","unit3"))
legend("topleft", fill = c('gray60','gray90'), legend = c('Iterative','Monolithic'), horiz = T,
       pt.cex=1.5)
dev.off()

# load("D:/VM/Tesi/my/ite.RData")

# Rprof(NULL)
# summaryRprof("Cshape.out", memory = "tseries", diff = TRUE)
detach(horseshoe2D)
# print("beta iterativo di 3 diverse simulazioni - corrispondente  iterativo vecchio")
# ite$f_betamat[,1:3]-f_betamat
# print("beta iterativo di 3 diverse simulazioni - corrispondente  kim")
kim$f_betamat[,1:3]-f_betamat