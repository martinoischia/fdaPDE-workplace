library(R.matlab)
vertices.data<-readMat("D:/connectome/conte/matlab_GIfTI/vertices.mat")
faces.data<-readMat("D:/connectome/conte/matlab_GIfTI/faces.mat")

library(fdaPDE)
mesh<-create.mesh.2.5D(vertices.data$vertices, faces.data$faces)
#plot(mesh)
FEMbasis <- create.FEM.basis(mesh)

##########create FC maps ##########
BrainRegions = read.csv("D:/VM/Tesi/codici jkim/R codes/brain/BrainRegions.csv")
dim(BrainRegions) #64984    69

ROI_index = BrainRegions$L_precuneus[1:(nrow(mesh$nodes))] #nrow(mesh$nodes): 32492
ROI_index = as.logical(ROI_index)
#length(ROI_index)
#length(which(ROI_index == TRUE)) #1469

### FC map function ####
FCmap_func<-function(datamatrix) {
  
  obj_name = "patient.dat"
  time_series <-  get(obj_name)
  mesh <- get("mesh")
  # print(dim(time_series)) 32492  1200
  print(length(which(is.na(time_series))))  #0
  
  time_series=t(time_series)
  ROI_ts=time_series[,ROI_index]
  # print(dim(ROI_ts)) #1200 1469
  
  cross_sec_avg_ROI_row=as.vector(rowMeans(ROI_ts))
  # print(length(cross_sec_avg_ROI_row)) #1200
  cor_nodes<-NULL
  for(j in 1:nrow(mesh$nodes)){
    cor_nodes<-c(cor_nodes,cor(cross_sec_avg_ROI_row, time_series[,j])) #have 1200 length in each nodes
  }
  # print(length(cor_nodes)) #32492
  print(length(which(is.na(cor_nodes)))) #2796
  
  r=cor_nodes
  z = 0.5 * log((1+r)/(1-r)) #fisher transformation
  
  datamatrix<-cbind(datamatrix,z)
}

#cor(c(1,2,3), c(1,2,3)) #1
#cor(c(1,2,3), c(6,2,8)) #0.3273268

#### open 1D data of patients #####
# 101107: not usuable (2808 NA whereas all other is 2796

### fMRI ####
load("D:/connectome/all_patients/connectivityMaps_stopped_at_badguy.RData")
files = as.character(read.csv("D:/connectome//all_patients//func_files", header = FALSE)[,1])
setwd("D:/connectome/all_patients/")

for (file in files){
	patient.dat <-read.csv(file,  nrows=32492, header = FALSE, sep=' ')[,1:1200]
	observations = FCmap_func(observations)
	dim(observations)
	rm(list=ls(pattern="patient"))
}

#### cortical thickness ####

files_covariates = as.character(read.csv("D:/connectome//all_patients//thick_files", header = FALSE)[,1])
for (file in files_covariates){
	thick.dat <-read.csv(file,  nrows=32492, header = FALSE, sep=' ')[,1]
	covariates=c(covariates, thick.dat)
	length(covariates)
	rm(list=ls(pattern="thick"))
}

save.image("D:/connectome/all_patients/connectivityMaps.RData")
