#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(PEcAn.all)
  library(PEcAnAssimSequential)
})

load_required <- function(path, expected = NULL) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  load(path, envir = .GlobalEnv)
  if (!is.null(expected) && !exists(expected, envir = .GlobalEnv)) {
    stop("Object `", expected, "` was not found after loading: ", path)
  }
}

subset_obs_sites <- function(obs_list, keep_ids) {
  lapply(obs_list, function(site_list) {
    if (is.null(site_list)) {
      return(site_list)
    }
    site_list[names(site_list) %in% keep_ids]
  })
}

apply_obs_floors <- function(obs_mean, obs_cov, mean_floor = 0.01, var_floor = 0.1) {
  for (ti in seq_along(obs_mean)) {
    if (is.null(obs_mean[[ti]]) || length(obs_mean[[ti]]) == 0) {
      next
    }

    for (sid in names(obs_mean[[ti]])) {
      mean_vec <- obs_mean[[ti]][[sid]]
      if (!is.null(mean_vec) && length(mean_vec) > 0) {
        mean_vec[mean_vec == 0] <- mean_floor
        obs_mean[[ti]][[sid]] <- mean_vec
      }

      cov_obj <- obs_cov[[ti]][[sid]]
      if (is.null(cov_obj) || length(cov_obj) == 0) {
        next
      }

      if (is.matrix(cov_obj) || is.data.frame(cov_obj)) {
        cov_mat <- as.matrix(cov_obj)
        d <- diag(cov_mat)
        d[d <= var_floor] <- var_floor
        diag(cov_mat) <- d
        obs_cov[[ti]][[sid]] <- cov_mat
      } else if (is.numeric(cov_obj)) {
        cov_obj[cov_obj <= var_floor] <- var_floor
        obs_cov[[ti]][[sid]] <- cov_obj
      }
    }
  }

  list(obs.mean = obs_mean, obs.cov = obs_cov)
}

project_dir <- "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/wishart_sda/"
setwd(project_dir)

paths <- list(
  settings = "/projectnb/dietzelab/guYANG/pecan/runners/test10/pecan_flux.RData",
  obs_mean = "/projectnb/dietzelab/guYANG/pecan/runners/test10/obs.mean.RData",
  obs_cov = "/projectnb/dietzelab/guYANG/pecan/runners/test10/obs.cov.RData",
  samples = "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/SDA_8k_site/samples.Rdata",
  met = "/projectnb/dietzelab/guYANG/Gap_fill/results/ec_3h.csv",
  matched = "/projectnb/dietzelab/guYANG/Validation/validation/matched_within_1km.csv"
)

job_template <- "/projectnb/dietzelab/guYANG/pecan/runners/test7/sipnet_template.job"
outdir <- "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_3/"
rundir <- file.path(outdir, "run")

load_required(paths$settings, expected = "settings")

settings$model$jobtemplate <- job_template
settings$outdir <- outdir
settings$host$prerun <- "module load R/4.4.0"
settings$host$rundir <- rundir

general.job <- list(cores = 28, folder.num = 80)
batch.settings <- list(
  general.job = general.job,
  qsub.cmd = "qsub -l h_rt=24:00:00 -l mem_per_core=4G -l buyin -pe omp @CORES@ -V -N @NAME@ -o @STDOUT@ -e @STDERR@ -S /bin/bash"
)
settings$state.data.assimilation$batch.settings <- batch.settings

settings$ensemble$size <- 25
settings <- PEcAn.settings::prepare.settings(settings)
settings <- PEcAn.settings::as.MultiSettings(settings)

load_required(paths$obs_mean, expected = "obs.mean")
load_required(paths$obs_cov, expected = "obs.cov")

resimet.df <- data.table::fread(paths$met)
close_points_df <- data.table::fread(paths$matched)

setDT(resimet.df)
setDT(close_points_df)

close_points_df <- close_points_df[order(min_dist_m), .SD[1], by = index]
resimet_2012_2017 <- resimet.df[
  !is.na(utc) &
    year(utc) >= 2012 &
    year(utc) <= 2017
]

site_ids_2012_2017 <- unique(resimet_2012_2017$Site_ID)
candidate_ids <- close_points_df[Site_ID %in% site_ids_2012_2017][order(min_dist_m)]

if (nrow(candidate_ids) < 31) {
  stop("Not enough candidate rows to select indices 31:66.")
}

pick_end <- min(66L, nrow(candidate_ids))
keep_ids <- as.character(candidate_ids[31:pick_end, index])

settings_list <- as.list(settings)
all_ids <- vapply(settings_list, function(s) as.character(s$run$site$id), character(1))
use_settings <- all_ids %in% keep_ids

if (!any(use_settings)) {
  stop("No settings site ids matched selected keep_ids.")
}

settings <- PEcAn.settings::as.MultiSettings(settings_list[use_settings])
obs.mean <- subset_obs_sites(obs.mean, keep_ids)
obs.cov <- subset_obs_sites(obs.cov, keep_ids)

floored_obs <- apply_obs_floors(obs.mean, obs.cov, mean_floor = 0.01, var_floor = 0.1)
obs.mean <- floored_obs$obs.mean
obs.cov <- floored_obs$obs.cov

load_required(paths$samples)
if (!exists("ensemble.samples")) {
  stop("Object `ensemble.samples` is missing from: ", paths$samples)
}

settings$state.data.assimilation$q.type <- "wishart"

Q <- NULL
restart <- NULL
pre_enkf_params <- NULL

control <- list(
  TimeseriesPlot = FALSE,
  OutlierDetection = FALSE,
  send_email = NULL,
  keepNC = FALSE,
  forceRun = TRUE,
  run_parallel = FALSE,
  MCMC.args = NULL,
  merge_nc = FALSE,
  execution = "qsub_parallel"
)

dir.create(settings$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(settings$host$rundir, recursive = TRUE, showWarnings = FALSE)

PEcAnAssimSequential:::sda.enkf.multisite(
  settings = settings,
  obs.mean = obs.mean,
  obs.cov = obs.cov,
  Q = Q,
  restart = restart,
  pre_enkf_params = pre_enkf_params,
  ensemble.samples = ensemble.samples,
  control = control
)
