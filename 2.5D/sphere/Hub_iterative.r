#### Test 2: sphere domain ####
#            locations != nodes 
#            with covariates
#            no BC
#            order FE = 1
library(fdaPDE)
rm(list=ls())
graphics.off()
load("D:\\VM\\Tesi\\debugging\\fdaPDEnocran\\data\\sphere2.5D.RData")
# mesh = sphere2.5D

nodes = sphere2.5D$nodes
triangles = sphere2.5D$triangles

## Create mesh from nodes and connectivity matrix:
mesh = create.mesh.2.5D(nodes = nodes, triangles = triangles)
plot(mesh)
snapshot3d( 'D:/VM/Tesi/my/sphere.png', height = 450, width = 450)
FEMbasis=create.FEM.basis(mesh)

# Test function 
f = function(x, y, z, t){
  phi = (asin(y/sqrt(x^2+y^2)))
  theta = acos(z/sqrt(x^2+y^2+z^2))
  # rho = 1
  
  sin(4*(1/2*sin(t^2*theta)*exp(-sin(theta)^2)+1)*theta)*cos(2*(1/2*cos(t*phi)*exp(-cos(phi)^2)+1)*phi)*cos(t)+ cos(2*t)
}

# Exact solution (pointwise at nodes)
# sol_exact=f(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3])
# plot(FEM(sol_exact, FEMbasis))

# Generate data locations on the sphere
set.seed(598944)
ndata = 100
locations = matrix(rnorm(ndata*3), ncol = 3)
locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)

# Generate covariate and data
num_units = 5
cov1 = runif(ndata, min = -1, max = 1)
DatiEsatti = NULL;
ts = seq(-1,2,length.out = num_units)
bi = seq(-20,20,length.out = num_units)
for (i in 1:num_units)
{
  DatiEsatti = cbind(DatiEsatti, f(locations[,1],locations[,2], locations[,3], ts[i]) + bi[i]*cov1)
}
cov2= matrix(rnorm(ndata*num_units), ncol=num_units);
DatiEsatti = DatiEsatti+ cov2
cov = cbind(rep(cov1,num_units), as.vector(cov2))

# Add error to simulate data
set.seed(7893475)
ran=range(DatiEsatti)
data = DatiEsatti + matrix(rnorm(ndata*num_units, mean=0, sd=0.05*abs(ran[2]-ran[1])), ncol = num_units)

# Project locations on the mesh
projected_locations = projection.points.2.5D(mesh, locations)

# Set smoothing parameter
# lambda = 10^seq(-4,-2,by=0.25)
lambda = 0.01

#### Test 2.1: Without GCV
# output_CPP<-smooth.FEM(observations=data, 
#                        locations = projected_locations, 
#                        covariates = cov1,
#                        FEMbasis=FEMbasis, 
#                        lambda=lambda[1])
# # plot(output_CPP$fit.FEM)
# output_CPP$solution$beta

#### Test 2.2: grid with exact GCV

tol=10^(-6)
mod2_info <- smooth.FEM.mixed(
    locations = projected_locations,
    observations = data,
    covariates = cov,
    FEMbasis = FEMbasis,
    random_effect = c(1),
    lambda = lambda,
    lambda.selection.lossfunction = 'GCV',
    DOF.evaluation = "exact",
    FLAG_ITERATIVE = TRUE, max.steps = 1000, threshold = tol, threshold_residual = tol, verbose = TRUE)

plot(table(mod2_info$iterations), ylab="counts", xlab="number of iterations")
# plot(log10(lambda), mod2_info$optimization$GCV_vector)
# plot(mod2_info$fit.FEM.mixed)

# mod2_info$beta
# mod2_info$b_i


# #### Test 2.3: grid with stochastic GCV
# output_CPP<-smooth.FEM(observations=data, locations = projected_locations, 
#                        covariates = cov1,
#                        FEMbasis=FEMbasis, lambda=lambda,
#                        lambda.selection.criterion='grid', DOF.evaluation='stochastic', lambda.selection.lossfunction='GCV')
# plot(log10(lambda), output_CPP$optimization$GCV_vector)
# plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))

# output_CPP$solution$beta