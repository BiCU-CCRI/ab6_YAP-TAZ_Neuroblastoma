# open the project ./workspace.Rproj first
# Set up the environment, installing all the required packages
# this might take a long time ~1h
renv::init()


if(!dir.exists("~/workspace/neuroblastoma/results/")) {
  dir.create("~/workspace/neuroblastoma/results/")
}

if(!dir.exists("~/workspace/neuroblastoma/results/RNA-seq")) {
  dir.create("~/workspace/neuroblastoma/results/RNA-seq")
}

if(!dir.exists("~/workspace/neuroblastoma/results/ATAC-seq")) {
  dir.create("~/workspace/neuroblastoma/results/ATAC-seq")
}


files <- list.files(path="../Output/datasets/", pattern=".zip$")
outDir <- "../Output/datasets/unzip"
for (i in files) {
  unzip(paste0("../Output/datasets/",i), exdir=outDir)
}
