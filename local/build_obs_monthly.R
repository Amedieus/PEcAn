# loading libraries
library(dplyr)
library(xts)
library(PEcAn.all)
library(purrr)
library(furrr)
library(lubridate)
library(nimble)
library(ncdf4)
library(PEcAnAssimSequential)
library(dplyr)
library(sp)
library(raster)
library(zoo)
library(ggplot2)
library(mnormt)
library(sjmisc)
library(stringr)
library(doParallel)
library(doSNOW)
library(data.table)
library(Kendall)
library(lgarch)
library(parallel)
library(foreach)
library(terra)
setwd("/projectnb/dietzelab/guYANG/pecan/")
## read settings xml file.
load("/projectnb/dietzelab/guYANG/pecan/pecan.Rdata")
# Change dir name and settings
settings$outdir      <- "/projectnb/dietzelab/guYANG/pecan/output/adjust/"
settings$rundir      <- file.path(settings$outdir, "run")
settings$modeloutdir <- file.path(settings$outdir, "out")
settings$host$rundir <- file.path(settings$outdir, "run")
settings$host$outdir <- file.path(settings$outdir, "out")
settings$host$folder <- file.path(settings$outdir, "out")
settings$ensemble$size <- 5
settings$state.data.assimilation$adjustment <- "FALSE"
settings$host$prerun <- "module load R/4.4.0"
###### Change Q type
settings$state.data.assimilation$q.type <- "wishart"
settings$state.data.assimilation$aqq.Init <- "1"
settings$state.data.assimilation$bqq.Init <- "1"
## Fix the multi output in one timestep bug
settings$model$jobtemplate <- "/projectnb/dietzelab/guYANG/pecan/runners/test7/sipnet_template.job"
# Load the selected sites
load("/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/sda_idx.Rdata")
all_ids  <- vapply(settings, \(s) as.character(s$run$site$id), "")
# Sub settings for this run
settings <- settings[all_ids %in% keep_ids]
settings <- PEcAn.settings::as.MultiSettings(settings)
# setup the batch job settings.
general.job <- list(cores = 28, folder.num = 80)
batch.settings = structure(list(
  general.job = general.job,
  qsub.cmd = "qsub -l h_rt=24:00:00 -l mem_per_core=4G -l buyin -pe omp @CORES@ -V -N @NAME@ -o @STDOUT@ -e @STDERR@ -S /bin/bash"
))
settings$state.data.assimilation$batch.settings <- batch.settings
# update settings with the actual PFTs.
settings <- PEcAn.settings::prepare.settings(settings)

## Modify OBS
settings$state.data.assimilation$Obs_Prep$outdir <- "/projectnb/dietzelab/guYANG/pecan/output/obs/"
settings$state.data.assimilation$Obs_Prep
# settings$state.data.assimilation$Obs_Prep$SMAP_SMP <- NULL

#### Change Landtrendr to GEDI
# buffer <- 0.005
# batch <- FALSE
# cores <- 22
# search_window <- "6 month"
# credential_path <- "~/.netrc"
# Obs_Prep <- settings$state.data.assimilation$Obs_Prep
# Obs_Prep$Landtrendr_AGB <- NULL
# Obs_Prep$GEDI_AGB <- list(
#   search_window   = "6 month",    
#   buffer          = "0.005",
#   batch           = "FALSE",
#   cores           = "22",
#   credential_path = "~/.netrc"
# )
# Obs_Prep$GEDI_AGB$timestep <- list(
#   unit = "year",
#   num  = "1"
# )
# settings$state.data.assimilation$Obs_Prep <- Obs_Prep

## Smoke test
Obs_Prep <- settings$state.data.assimilation$Obs_Prep

# Remove global timestep, forcing SDA_OBS_Assembler to use product-level timestep
Obs_Prep$timestep <- NULL

# Give each product its own start/end/timestep
for (nm in c("Landtrendr_AGB", "MODIS_LAI", "SMAP_SMP", "Soilgrids_SoilC")) {
  Obs_Prep[[nm]]$start.date <- "2012-07-15"
  Obs_Prep[[nm]]$end.date   <- "2024-07-15"
}

Obs_Prep$Landtrendr_AGB$timestep <- list(unit = "year",  num = "1")
Obs_Prep$MODIS_LAI$timestep      <- list(unit = "month", num = "1")
Obs_Prep$SMAP_SMP$timestep       <- list(unit = "month", num = "1")
Obs_Prep$Soilgrids_SoilC$timestep <- list(unit = "year", num = "1")
settings$state.data.assimilation$Obs_Prep <- Obs_Prep

OBS <- SDA_OBS_Assembler(settings)
obs.mean <- OBS$obs.mean
obs.cov  <- OBS$obs.cov
names(obs.mean)
