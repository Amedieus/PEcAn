### Random pick sites
# resimet.df <- fread("/projectnb/dietzelab/guYANG/Gap_fill/results/ec_3h.csv")
# close_points_df <- fread("/projectnb/dietzelab/guYANG/Validation/validation/matched_within_1km.csv")
# setDT(close_points_df)
# close_points_df <- close_points_df[order(min_dist_m), .SD[1], by = index]
# set.seed(123)
# keep_ids <- sample(as.character(close_points_df$index), 36)
# save(keep_ids, file = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/sda_idx.Rdata")

###
### All Ameriflux sites
resimet.df <- fread("/projectnb/dietzelab/guYANG/Gap_fill/results/ec_3h.csv")
close_points_df <- fread("/projectnb/dietzelab/guYANG/Validation/validation/matched_within_1km.csv")
setDT(close_points_df)
close_points_df <- close_points_df[order(min_dist_m), .SD[1], by = index]
gc()
keep_ids <- as.character(close_points_df$index)
save(keep_ids, file = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/sda_idx.Rdata")

### NEON Sites
resimet.df <- fread("/projectnb/dietzelab/guYANG/Gap_fill/results/ec_3h.csv")
close_points_df <- fread("/projectnb/dietzelab/guYANG/Validation/validation/matched_within_1km.csv")
setDT(close_points_df)
index_us_x <- close_points_df[
  grepl("^US-x", Site_ID),
  index
]
keep_ids <- sample(as.character(index_us_x))
save(keep_ids, file = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/sda_idx.Rdata")

index_td <- map[
  site_id %in% dcp[pft == "temperate.deciduous.HPDA", site],
  matched_index
]

index_td <- c(
  map[
    site_id %in% dcp[pft == "temperate.deciduous.HPDA", site],
    matched_index
  ],
  262, 705, 706, 707
)

keep_ids <- sample(as.character(index_td))
save(keep_ids, file = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/sda_idx.Rdata")
