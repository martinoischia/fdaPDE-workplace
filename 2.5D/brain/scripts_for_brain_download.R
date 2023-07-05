library ( neurohcp )
library ( dplyr )
knitr::opts_chunk$set( echo = TRUE, cache = FALSE, comment ="")
ACCESS_KEY ="AKIAXO65CT57AJSH5BMY"
SECRET_KEY ="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"
set_aws_api_key ( access_key = ACCESS_KEY, secret_key = SECRET_KEY
)
neurohcp::bucketlist ( verbose = FALSE)
ids_with_dwi = hcp_900_scanning_info %>%
filter ( scan_type %in% "dMRI" ) %>%
select ( id ) %>%
unique

patient_id = ids_with_dwi$id
length ( patient_id )

fmri=list()
thick=list()
folder='D:/connectome/all_patients/'
#'D:\connectome\all_patients\'
for (chosen_patient in patient_id[1:60]){
	print(chosen_patient)
	hcp_fmri_file =paste0 ( 'HCP/' , chosen_patient )
	hcp_fmri_file = paste0 ( hcp_fmri_file , '/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas.dtseries.nii')
	hcp_thick_file =paste0 ( 'HCP/' , chosen_patient )
	hcp_thick_file = paste0 ( hcp_thick_file , '/MNINonLinear/fsaverage_LR32k/' )
	hcp_thick_file = paste0 ( hcp_thick_file , chosen_patient )
	hcp_thick_file = paste0 ( hcp_thick_file , '.thickness.32k_fs_LR.dscalar.nii' )
	fmri = c(fmri, download_hcp_file ( hcp_fmri_file ,destfile=paste0(folder, chosen_patient), verbose = TRUE))
	thick = c(thick, download_hcp_file ( hcp_thick_file, verbose = TRUE))
}

# this is for CMD

wb_command -cifti-separate "D:\connectome\all_patients\   rfMRI_REST1_LR_Atlas.dtseries.nii" COLUMN -metric CORTEX_LEFT "D:\connectome\all_patients\    101410.L.func.gii" # correct numbers pls
wb_command -cifti-separate "C:\Users\ischi\AppData\Local\Temp\RtmpUf7WRI\100307.thickness.32k_fs_LR.dscalar.nii" COLUMN -metric CORTEX_LEFT "D:/connectome/101410.L.thickness.32k_fs_LR.func.gii"

# this for AFNI dont stick to this look ubuntu 
for file in /mnt/d/connectome/all_patients/*.L.func.gii; do
echo gifti_tool -infile $file -write_1D $(basename $file .gii).L.1D
done

gifti_tool -infile /mnt/d/connectome/101410.L.func.gii -write_1D /mnt/d/connectome/100307.L.1D
gifti_tool -infile /mnt/d/connectome/101410.L.thickness.32k_fs_LR.func.gii -write_1D /mnt/d/connectome/100307.L.thickness.32k_fs_LR.func.1D
#this for R 

fMRI1.inputFile <- "D:/connectome/100307.L.1D"
patient1.dat <- read.csv ( fMRI1.inputFile , nrows=32492 , header = FALSE, sep=' ') [ , 1 : 1200 ]
thick1.inputFile <- "D:/connectome/100307.L.thickness.32k_fs_LR.func.1D"
thick1.dat <- read.csv ( thick1.inputFile , nrows=32492 , header = FALSE, sep=' ') [ , 1 ]

