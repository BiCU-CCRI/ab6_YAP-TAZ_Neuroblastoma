# This script is designed to set up the environment for this R project.
# Before you begin, make sure that you are you activated the workspace.Rproj
# By openning the project ./workspace.Rproj first

# Set up the environment, installing all the required packages using renv
# This might take a long time ~1h
renv::activate(project = "/home/rstudio/workspace/")
renv::restore(project = "/home/rstudio/workspace/", prompt=FALSE)


# Pull the missing data from SRA/GEO
# TODO - write this module when GEO data is public

# Create output directories for results
list_of_dirs_to_create <- list(
  "~/workspace/neuroblastoma/results/",
  "~/workspace/neuroblastoma/results/RNA-seq",
  "~/workspace/neuroblastoma/results/ATAC-seq",
  "~/workspace/neuroblastoma/results/CnR"
)
for(dir_to_create in list_of_dirs_to_create){
  if(!dir.exists(dir_to_create)) {
    dir.create(dir_to_create)
  }
}
rm(list_of_dirs_to_create)

# TODO - Is this code necessary?
# files <- list.files(path="../Output/datasets/", pattern=".zip$")
# outDir <- "../Output/datasets/unzip"
# for (i in files) {
#   unzip(paste0("../Output/datasets/",i), exdir=outDir)
# }
