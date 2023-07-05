library(fdaPDE)
rm(list=ls())
graphics.off()
set.seed(400000)

x = seq(0,1, length.out = 40)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations)
# plot(mesh)
nnodes = dim(mesh$nodes)[1]
numunits = 3
FEMbasis = create.FEM.basis(mesh)

f = function(x, y, z = 1)
{
    coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
    sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}




# cov1 <- sin(2 * pi * mesh$nodes[, 1]) * cos(2 * pi * mesh$nodes[, 2])
cov1 <- rnorm(nnodes, mean = 0, sd = 2)
cov2 <- rnorm(nnodes*numunits, mean = 0, sd = 2)
W <- cbind(rep(cov1, numunits ), cov2)
beta_exact <- cbind( rnorm( numunits , sd = 8 ), rep(0.5, numunits))

mixed_obs = matrix(nrow = nnodes, ncol = 0)
for (i in 1:numunits)
{
    mixed_obs = cbind( mixed_obs , W[(nnodes * (i-1) + 1) : (nnodes * i),] %*% beta_exact[i,] )
}

f_exact = matrix(nrow = nnodes, ncol = 0)
error = matrix(nrow = nnodes, ncol = 0)

for (i in runif(numunits, 0, 100))
{
    f_exact = cbind ( f_exact , f( mesh$nodes[,1], mesh$nodes[,2], z = i))
    #error = cbind( error, rnorm( nnodes , sd = 0.05 * 2)) # 2 is the range of f
}
mixed_obs = mixed_obs + f_exact #+ error


lambda <- 10^seq(-2, 1, by = 0.1)
GCVMETHODFLAG <- 'stochastic' 
# GCVMETHODFLAG <-"exact"
tol <- 1e-6

mixed_mod <- smooth.FEM.mixed(
    locations = mesh$nodes,
    observations = mixed_obs,
    covariates = W,
    FEMbasis = FEMbasis,
    random_effect = c(1),
    lambda = lambda,
    GCV = TRUE,
    GCVmethod = 'Exact'
)

mixed_mod$beta [1 , mixed_mod$bestlambda] - beta_exact [1, 2]
mixed_mod$b_i [ , mixed_mod$bestlambda] - beta_exact [, 1]