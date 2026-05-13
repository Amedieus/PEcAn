library(dplyr)

source("/projectnb/dietzelab/guYANG/pecan/pecan/local/diagnose_functions.R")

res_diag <- run_sda_diagnostic_pipeline(
  sda_dir = "/projectnb/dietzelab/guYANG/pecan/output/4obs_yearly",
  out_dir = "/projectnb/dietzelab/guYANG/pecan/output/4obs_yearly",
  file_index = 1:13,
  start_year = 2012,
  soilmoist_div10 = FALSE,
  prefix = "6obs_yearly",
  save_outputs = TRUE,
  make_site_nis_plots = TRUE
)

final_summary_all <- res_diag$final_summary_all
diag_summary_all  <- res_diag$diag_summary_all
obs_long_all      <- res_diag$obs_long_all
nis_dat           <- res_diag$nis_dat

res_diag$plots$scatter$plot
res_diag$plots$NIS_by_year$plot
res_diag$plots$R2_by_variable$plot
res_diag$plots$R2_improvement_by_variable$plot
res_diag$plots$R2_by_year_variable$plot
res_diag$plots$R2_improvement_by_year_variable$plot

###############################################################
#### Monthly ####

library(dplyr)

source("/projectnb/dietzelab/guYANG/pecan/pecan/local/diagnose_functions_monthly.R")

assim_dates <- seq(
  as.Date("2012-07-15"),
  as.Date("2021-07-15"),
  by = "1 month"
)

res_diag <- run_sda_diagnostic_pipeline_monthly(
  sda_dir = "/projectnb/dietzelab/guYANG/pecan/output/5obs_monthly",
  out_dir = "/projectnb/dietzelab/guYANG/pecan/output/5obs_monthly",
  file_index = seq_along(assim_dates),
  assim_dates = assim_dates,
  soilmoist_div10 = FALSE,
  prefix = "5obs_monthly",
  save_outputs = TRUE,
  make_site_nis_plots = TRUE
)

final_summary_all <- res_diag$final_summary_all
diag_summary_all  <- res_diag$diag_summary_all
obs_long_all      <- res_diag$obs_long_all
nis_dat           <- res_diag$nis_dat

res_diag$plots$scatter$plot
res_diag$plots$NIS_by_date$plot
res_diag$plots$R2_by_variable$plot
res_diag$plots$R2_improvement_by_variable$plot
res_diag$plots$R2_by_date_variable$plot
res_diag$plots$R2_improvement_by_date_variable$plot
res_diag$plots$R2_by_calendar_month_variable$plot