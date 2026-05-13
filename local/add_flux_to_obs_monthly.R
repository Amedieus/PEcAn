library(data.table)
library(lubridate)

# ============================================================
# 0. User settings
# ============================================================

obs_rdata_dir <- "/projectnb/dietzelab/guYANG/pecan/output/obs/Rdata"

# Existing obs.mean / obs.cov files
# Change these if your current monthly obs files have different names.
in_obs_mean_file <- file.path(obs_rdata_dir, "obs.mean.Rdata")
in_obs_cov_file  <- file.path(obs_rdata_dir, "obs.cov.Rdata")

# EC gap-filled ensemble flux files
nee_file <- "/projectnb/dietzelab/guYANG/Gap_fill/results/ens_ec_3h.csv"
qle_file <- "/projectnb/dietzelab/guYANG/Gap_fill/results/le_ens_ec_3h.csv"

# Site_ID to SDA index matching file
match_file <- "/projectnb/dietzelab/guYANG/Validation/validation/matched_within_1km.csv"

# Output files
out_obs_mean_file <- file.path(obs_rdata_dir, "obs.mean.monthly.with_flux.Rdata")
out_obs_cov_file  <- file.path(obs_rdata_dir, "obs.cov.monthly.with_flux.Rdata")

# Monthly SDA observation window
start_date <- as.Date("2012-07-15")
end_date   <- as.Date("2024-07-15")

# SDA / XML state order
state_vars_xml <- c(
  "AbvGrndWood",
  "NEE",
  "Qle",
  "LAI",
  "TotSoilCarb"
)

# Inflation for NEE/Qle observation variances
flux_var_inflation <- 1.1


# ============================================================
# 1. Load existing obs.mean and obs.cov
# ============================================================

load(in_obs_mean_file)  # should load object: obs.mean
load(in_obs_cov_file)   # should load object: obs.cov

stopifnot(exists("obs.mean"))
stopifnot(exists("obs.cov"))


# ============================================================
# 2. Helper: normalize obs list dates
# ============================================================

get_obs_date_table <- function(obs.mean) {
  data.table(
    original_name = names(obs.mean),
    date = as.Date(substr(names(obs.mean), 1, 10))
  )
}


# ============================================================
# 3. Helper: aggregate 3-hourly ensemble EC data to monthly SDA time
# ============================================================

aggregate_flux_to_monthly <- function(flux_file,
                                      match_file,
                                      flux_name = c("NEE", "Qle"),
                                      start_date = as.Date("2012-07-15"),
                                      end_date   = as.Date("2024-07-15")) {
  
  flux_name <- match.arg(flux_name)
  
  flux_df <- fread(flux_file)
  close_points_df <- fread(match_file)
  
  setDT(flux_df)
  setDT(close_points_df)
  
  # Map Site_ID to the SDA obs.mean site index.
  idx_map <- close_points_df[, .(Site_ID, index)]
  idx_map[, Site_ID := as.character(Site_ID)]
  idx_map[, index := as.character(index)]
  setkey(idx_map, Site_ID)
  
  flux_df[, Site_ID := as.character(Site_ID)]
  flux_df <- idx_map[flux_df, on = "Site_ID"]
  
  dt <- as.data.table(flux_df)
  
  # Clean time
  dt <- dt[!is.na(utc)]
  dt[, utc := as.POSIXct(utc, tz = "UTC")]
  
  # Remove rows without matched SDA index
  dt <- dt[!is.na(index)]
  
  # ------------------------------------------------------------
  # Monthly assimilation time anchored on the 15th.
  #
  # Day 1-15    -> current month 15th
  # Day 16-end  -> next month 15th
  #
  # Examples:
  # 2012-07-01 to 2012-07-15 -> 2012-07-15
  # 2012-07-16 to 2012-08-15 -> 2012-08-15
  # ------------------------------------------------------------
  
  dt[, time := as.Date(floor_date(utc, unit = "month") + days(14))]
  
  dt[day(utc) > 15,
     time := as.Date(floor_date(utc %m+% months(1), unit = "month") + days(14))]
  
  # Keep only SDA monthly window
  dt <- dt[time >= start_date & time <= end_date]
  
  # Ensemble columns: ens01, ens02, ...
  ens_cols <- grep("^ens\\d{2}$", names(dt), value = TRUE)
  ens_cols <- sort(ens_cols)
  
  if (length(ens_cols) < 1) {
    stop("No ensemble columns found. Expected columns like ens01, ens02, ...")
  }
  
  if (!("ens_mean" %in% names(dt))) {
    stop("Column ens_mean is missing from flux file: ", flux_file)
  }
  
  # Monthly mean for each site index x monthly SDA time
  monthly_df <- dt[
    ,
    c(
      lapply(.SD, mean, na.rm = TRUE),
      list(ens_mean = mean(ens_mean, na.rm = TRUE))
    ),
    by = .(index, time),
    .SDcols = ens_cols
  ]
  
  # Remove rows where all ensemble values and ens_mean are NA
  cols_check <- c(ens_cols, "ens_mean")
  monthly_df <- monthly_df[
    rowSums(!is.na(monthly_df[, ..cols_check])) > 0
  ]
  
  # Rename ensemble columns
  old_names <- c(ens_cols, "ens_mean")
  new_names <- c(paste0(flux_name, "_", ens_cols), flux_name)
  setnames(monthly_df, old = old_names, new = new_names)
  
  flux_ens_cols <- paste0(flux_name, "_", ens_cols)
  
  # Ensemble spread based on monthly ensemble means
  monthly_df[
    ,
    paste0(flux_name, "_sd") := apply(.SD, 1, sd, na.rm = TRUE),
    .SDcols = flux_ens_cols
  ]
  
  monthly_df[
    ,
    paste0(flux_name, "_var") := apply(.SD, 1, var, na.rm = TRUE),
    .SDcols = flux_ens_cols
  ]
  
  # Keep only needed columns
  keep_cols <- c(
    "index",
    "time",
    flux_name,
    paste0(flux_name, "_sd"),
    paste0(flux_name, "_var")
  )
  
  monthly_df <- monthly_df[, ..keep_cols]
  monthly_df[, index := as.character(index)]
  monthly_df[, time := as.Date(time)]
  
  setorder(monthly_df, index, time)
  
  return(monthly_df)
}


# ============================================================
# 4. Build monthly flux observation table: NEE + Qle
# ============================================================

nee_monthly_df <- aggregate_flux_to_monthly(
  flux_file = nee_file,
  match_file = match_file,
  flux_name = "NEE",
  start_date = start_date,
  end_date = end_date
)

qle_monthly_df <- aggregate_flux_to_monthly(
  flux_file = qle_file,
  match_file = match_file,
  flux_name = "Qle",
  start_date = start_date,
  end_date = end_date
)

monthly_flux_df <- merge(
  nee_monthly_df,
  qle_monthly_df,
  by = c("index", "time"),
  all = TRUE
)

monthly_flux_df[, index := as.character(index)]
monthly_flux_df[, time := as.Date(time)]

setorder(monthly_flux_df, index, time)

message("Monthly flux table:")
print(dim(monthly_flux_df))
print(head(monthly_flux_df))


# ============================================================
# 5. Diagnostics before insertion
# ============================================================

obs_date_table <- get_obs_date_table(obs.mean)

message("Date overlap between obs.mean and monthly flux table:")
date_overlap <- intersect(obs_date_table$date, unique(monthly_flux_df$time))
print(length(date_overlap))
print(head(date_overlap, 20))

first_obs_name <- names(obs.mean)[1]
obs_site_names <- names(obs.mean[[first_obs_name]])

message("Example obs.mean date names:")
print(head(names(obs.mean)))

message("Example obs.mean site names:")
print(head(obs_site_names))

message("Example monthly_flux_df index values:")
print(head(unique(monthly_flux_df$index)))

message("Site-name overlap:")
print(length(intersect(obs_site_names, unique(monthly_flux_df$index))))


# ============================================================
# 6. Add NEE and Qle to obs.mean by monthly date and site name
# ============================================================

add_flux_to_obsmean_by_name <- function(obs.mean,
                                        monthly_flux_df,
                                        state_vars_xml) {
  
  obs.mean_new <- obs.mean
  
  obs_date_table <- get_obs_date_table(obs.mean_new)
  setkey(obs_date_table, date)
  
  flux <- as.data.table(copy(monthly_flux_df))
  flux[, time := as.Date(time)]
  flux[, site_chr := as.character(index)]
  
  n_add_nee <- 0L
  n_add_qle <- 0L
  n_skip_date <- 0L
  n_skip_site <- 0L
  n_skip_no_flux <- 0L
  
  for (r in seq_len(nrow(flux))) {
    
    d_date <- flux$time[r]
    sid <- flux$site_chr[r]
    
    obs_name <- obs_date_table[J(d_date), original_name]
    
    if (length(obs_name) == 0 || is.na(obs_name[1])) {
      n_skip_date <- n_skip_date + 1L
      next
    }
    
    # If duplicate obs names map to same Date, use the first one.
    d <- obs_name[1]
    
    if (!(sid %in% names(obs.mean_new[[d]]))) {
      n_skip_site <- n_skip_site + 1L
      next
    }
    
    nee <- flux$NEE[r]
    qle <- flux$Qle[r]
    
    if (is.na(nee) && is.na(qle)) {
      n_skip_no_flux <- n_skip_no_flux + 1L
      next
    }
    
    df0 <- obs.mean_new[[d]][[sid]]
    
    if (is.null(df0)) {
      df0 <- data.frame()
    }
    
    if (!is.data.frame(df0)) {
      df0 <- as.data.frame(df0)
    }
    
    if (!is.na(nee)) {
      df0$NEE <- nee
      n_add_nee <- n_add_nee + 1L
    }
    
    if (!is.na(qle)) {
      df0$Qle <- qle
      n_add_qle <- n_add_qle + 1L
    }
    
    # Reorder variables according to SDA/XML order
    vars_present <- names(df0)
    vars_ordered <- state_vars_xml[state_vars_xml %in% vars_present]
    extra_vars <- setdiff(vars_present, vars_ordered)
    
    df0 <- df0[, c(vars_ordered, extra_vars), drop = FALSE]
    
    obs.mean_new[[d]][[sid]] <- df0
  }
  
  message("Added NEE values: ", n_add_nee)
  message("Added Qle values: ", n_add_qle)
  message("Skipped rows because date not found: ", n_skip_date)
  message("Skipped rows because site not found: ", n_skip_site)
  message("Skipped rows because both NEE and Qle are NA: ", n_skip_no_flux)
  
  return(obs.mean_new)
}


obs.mean_new <- add_flux_to_obsmean_by_name(
  obs.mean = obs.mean,
  monthly_flux_df = monthly_flux_df,
  state_vars_xml = state_vars_xml
)


# ============================================================
# 7. Add NEE and Qle variances to obs.cov by monthly date and site name
# ============================================================

add_flux_to_obscov_by_name <- function(obs.mean_old,
                                       obs.cov_old,
                                       obs.mean_new,
                                       monthly_flux_df,
                                       state_vars_xml,
                                       inflation = 1.1) {
  
  obs.cov_new <- obs.cov_old
  
  obs_date_table <- get_obs_date_table(obs.mean_new)
  setkey(obs_date_table, date)
  
  flux <- as.data.table(copy(monthly_flux_df))
  flux[, time := as.Date(time)]
  flux[, site_chr := as.character(index)]
  setkey(flux, time, site_chr)
  
  n_add_nee_var <- 0L
  n_add_qle_var <- 0L
  n_drop_missing_var <- 0L
  
  for (d in names(obs.mean_new)) {
    
    d_date <- as.Date(substr(d, 1, 10))
    
    for (sid in names(obs.mean_new[[d]])) {
      
      y_new <- obs.mean_new[[d]][[sid]]
      
      if (is.null(y_new) || !is.data.frame(y_new) || ncol(y_new) == 0) {
        obs.cov_new[[d]][[sid]] <- matrix(0, 0, 0)
        next
      }
      
      vars_new <- names(y_new)
      
      # Old observation mean/cov for old variables
      y_old <- NULL
      R_old <- NULL
      
      if (d %in% names(obs.mean_old) &&
          sid %in% names(obs.mean_old[[d]])) {
        y_old <- obs.mean_old[[d]][[sid]]
      }
      
      if (d %in% names(obs.cov_old) &&
          sid %in% names(obs.cov_old[[d]])) {
        R_old <- obs.cov_old[[d]][[sid]]
      }
      
      old_var_diag <- numeric(0)
      
      if (!is.null(y_old) &&
          is.data.frame(y_old) &&
          ncol(y_old) > 0 &&
          !is.null(R_old)) {
        
        R_old <- as.matrix(R_old)
        old_vars <- names(y_old)
        
        if (nrow(R_old) == length(old_vars) &&
            ncol(R_old) == length(old_vars)) {
          old_var_diag <- diag(R_old)
          names(old_var_diag) <- old_vars
        }
      }
      
      # New variance vector in the same order as obs.mean_new
      vv <- rep(NA_real_, length(vars_new))
      names(vv) <- vars_new
      
      # Fill old variables from old obs.cov
      common_old <- intersect(names(old_var_diag), vars_new)
      vv[common_old] <- old_var_diag[common_old]
      
      # Fill NEE/Qle from monthly flux ensemble spread
      row_flux <- flux[J(d_date, as.character(sid)), nomatch = 0]
      
      if (nrow(row_flux) > 0) {
        row_flux <- row_flux[1]
        
        if ("NEE" %in% vars_new && !is.na(row_flux$NEE_var)) {
          vv["NEE"] <- row_flux$NEE_var * inflation
          n_add_nee_var <- n_add_nee_var + 1L
        }
        
        if ("Qle" %in% vars_new && !is.na(row_flux$Qle_var)) {
          vv["Qle"] <- row_flux$Qle_var * inflation
          n_add_qle_var <- n_add_qle_var + 1L
        }
      }
      
      # If a variable has mean but no variance, drop it to avoid dimension mismatch
      ok <- !is.na(vv)
      
      if (!all(ok)) {
        vars_keep <- names(vv)[ok]
        vars_drop <- names(vv)[!ok]
        
        if (length(vars_drop) > 0) {
          n_drop_missing_var <- n_drop_missing_var + length(vars_drop)
        }
        
        y_new <- y_new[, vars_keep, drop = FALSE]
        obs.mean_new[[d]][[sid]] <- y_new
        
        vv <- vv[ok]
        vars_new <- names(vv)
      }
      
      if (length(vv) == 0) {
        obs.cov_new[[d]][[sid]] <- matrix(0, 0, 0)
      } else {
        R_new <- diag(as.numeric(vv), nrow = length(vv), ncol = length(vv))
        dimnames(R_new) <- list(vars_new, vars_new)
        obs.cov_new[[d]][[sid]] <- R_new
      }
    }
  }
  
  message("Added NEE variances: ", n_add_nee_var)
  message("Added Qle variances: ", n_add_qle_var)
  message("Dropped variables because variance was missing: ", n_drop_missing_var)
  
  return(list(
    obs.mean = obs.mean_new,
    obs.cov = obs.cov_new
  ))
}


tmp <- add_flux_to_obscov_by_name(
  obs.mean_old = obs.mean,
  obs.cov_old = obs.cov,
  obs.mean_new = obs.mean_new,
  monthly_flux_df = monthly_flux_df,
  state_vars_xml = state_vars_xml,
  inflation = flux_var_inflation
)

obs.mean_new <- tmp$obs.mean
obs.cov_new  <- tmp$obs.cov


# ============================================================
# 8. Final format checks
# ============================================================

check_obs_pair <- function(obs.mean, obs.cov) {
  
  bad <- list()
  k <- 1L
  
  for (d in names(obs.mean)) {
    
    if (!(d %in% names(obs.cov))) {
      bad[[k]] <- data.frame(
        time = d,
        site = NA,
        issue = "date missing in obs.cov"
      )
      k <- k + 1L
      next
    }
    
    for (sid in names(obs.mean[[d]])) {
      
      y <- obs.mean[[d]][[sid]]
      
      if (!(sid %in% names(obs.cov[[d]]))) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = "site missing in obs.cov"
        )
        k <- k + 1L
        next
      }
      
      R <- obs.cov[[d]][[sid]]
      
      if (is.null(y) || !is.data.frame(y) || ncol(y) == 0) {
        next
      }
      
      if (is.null(R)) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = "R is NULL"
        )
        k <- k + 1L
        next
      }
      
      R <- as.matrix(R)
      
      if (ncol(y) != nrow(R) || ncol(y) != ncol(R)) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = paste0(
            "dim mismatch: y=",
            ncol(y),
            ", R=",
            paste(dim(R), collapse = "x")
          )
        )
        k <- k + 1L
        next
      }
      
      if (!identical(names(y), colnames(R))) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = paste0(
            "name mismatch: y=",
            paste(names(y), collapse = ","),
            " | R=",
            paste(colnames(R), collapse = ",")
          )
        )
        k <- k + 1L
      }
    }
  }
  
  if (length(bad) == 0) {
    message("obs.mean and obs.cov are aligned.")
    return(invisible(NULL))
  }
  
  return(rbindlist(bad, fill = TRUE))
}


check_result <- check_obs_pair(obs.mean_new, obs.cov_new)

if (!is.null(check_result)) {
  print(check_result)
  stop("obs.mean and obs.cov are not aligned. Fix before saving.")
}


# ============================================================
# 9. Confirm NEE and Qle are present
# ============================================================

vars_all <- unique(unlist(lapply(obs.mean_new, function(x) {
  unlist(lapply(x, names))
})))

message("Variables found in obs.mean_new:")
print(vars_all)

n_nee <- sum(unlist(lapply(obs.mean_new, function(date_list) {
  lapply(date_list, function(df) {
    is.data.frame(df) && "NEE" %in% names(df)
  })
})))

n_qle <- sum(unlist(lapply(obs.mean_new, function(date_list) {
  lapply(date_list, function(df) {
    is.data.frame(df) && "Qle" %in% names(df)
  })
})))

message("Number of site-time entries containing NEE:")
print(n_nee)

message("Number of site-time entries containing Qle:")
print(n_qle)


# ============================================================
# 10. Extra diagnostics: monthly counts
# ============================================================

nee_by_date <- rbindlist(lapply(names(obs.mean_new), function(d) {
  data.table(
    time = d,
    date = as.Date(substr(d, 1, 10)),
    n_NEE = sum(unlist(lapply(obs.mean_new[[d]], function(df) {
      is.data.frame(df) && "NEE" %in% names(df)
    }))),
    n_Qle = sum(unlist(lapply(obs.mean_new[[d]], function(df) {
      is.data.frame(df) && "Qle" %in% names(df)
    })))
  )
}), fill = TRUE)

nee_by_date <- nee_by_date[date >= start_date & date <= end_date]
setorder(nee_by_date, date)

message("Monthly NEE/Qle insertion counts:")
print(nee_by_date)


# ============================================================
# 11. Save using object names expected by SDA
# ============================================================

obs.mean <- obs.mean_new
obs.cov  <- obs.cov_new

save(obs.mean, file = out_obs_mean_file)
save(obs.cov,  file = out_obs_cov_file)

message("Saved obs.mean to: ", out_obs_mean_file)
message("Saved obs.cov to: ", out_obs_cov_file)