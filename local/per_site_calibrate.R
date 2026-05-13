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

## Calibrated AGB 
agb_params <- fread("/projectnb/dietzelab/menglai/for_ppl/files/param_agb.csv")
site_ord <- fread("/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/matched_df.csv")

samples_file <- "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/SDA_8k_site/samples.Rdata" 
out_file     <- "/projectnb/dietzelab/guYANG/pecan/cali_samples.Rdata"

e <- new.env()
load(samples_file, envir = e)
stopifnot(exists("ensemble.samples", envir = e))
ensemble.samples <- get("ensemble.samples", envir = e)

target_pft <- "temperate.deciduous.HPDA"
stopifnot(target_pft %in% names(ensemble.samples))

old_df <- ensemble.samples[[target_pft]]
stopifnot(is.data.frame(old_df))

stopifnot(exists("agb_params"))
stopifnot(is.data.frame(agb_params))

# Delete site info
global_cols <- names(agb_params)[!grepl("^site_\\d+_", names(agb_params))]
agb_global <- agb_params[, ..global_cols]
common_cols <- intersect(names(old_df), names(agb_global))
if (length(common_cols) == 0) {
  stop("no param match found")
}

message("going to match：\n", paste(common_cols, collapse = ", "))

n_old <- nrow(old_df)
n_new <- nrow(agb_global)

## Deal with ensemble size
if (n_new == n_old) {
  idx <- seq_len(n_old)
} else if (n_new < n_old) {
  # ens > 50
  idx <- rep(seq_len(n_new), length.out = n_old)
  message("注意：agb_params 行数小于 ensemble 行数，已循环复用。")
} else {
  # ens < 50
  idx <- seq_len(n_old)
}

patched_df <- as.data.frame(old_df)
agb_global_df <- as.data.frame(agb_global)

patched_df[, common_cols] <- agb_global_df[idx, common_cols, drop = FALSE]
# 可选：确保 numeric
for (cc in common_cols) patched_df[[cc]] <- as.numeric(patched_df[[cc]])

ensemble.samples[[target_pft]] <- patched_df

# ==== 6) 写回环境并保存 ====
assign("ensemble.samples", ensemble.samples, envir = e)

# 尽量保留原 samples.Rdata 里的其他对象
obj_names <- ls(e, all.names = TRUE)
save(list = obj_names, file = out_file, envir = e)

message("完成。输出文件：", out_file)



