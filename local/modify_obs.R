drop_state_vars_from_obs <- function(vars_drop,
                                     obs_mean_in = NULL,
                                     obs_cov_in  = NULL,
                                     save_files = FALSE,
                                     out_mean_file = NULL,
                                     out_cov_file  = NULL,
                                     assign_to_global = FALSE,
                                     envir = parent.frame()) {
  
  vars_drop <- unique(as.character(vars_drop))
  
  if (length(vars_drop) == 0) {
    stop("vars_drop is empty.")
  }
  
  # If not provided explicitly, take obs.mean / obs.cov from current environment
  if (is.null(obs_mean_in)) {
    if (!exists("obs.mean", envir = envir, inherits = TRUE)) {
      stop("obs.mean is not found. Please load obs.mean first or provide obs_mean_in.")
    }
    obs_mean_in <- get("obs.mean", envir = envir, inherits = TRUE)
  }
  
  if (is.null(obs_cov_in)) {
    if (!exists("obs.cov", envir = envir, inherits = TRUE)) {
      stop("obs.cov is not found. Please load obs.cov first or provide obs_cov_in.")
    }
    obs_cov_in <- get("obs.cov", envir = envir, inherits = TRUE)
  }
  
  obs.mean_new <- obs_mean_in
  obs.cov_new  <- obs_cov_in
  
  n_removed_mean <- 0L
  n_removed_cov  <- 0L
  
  normalize_R <- function(R, vars_old) {
    if (is.null(R)) return(NULL)
    
    if (length(vars_old) == 1 && is.null(dim(R))) {
      R_mat <- matrix(as.numeric(R), nrow = 1, ncol = 1)
      dimnames(R_mat) <- list(vars_old, vars_old)
      return(R_mat)
    }
    
    R_mat <- as.matrix(R)
    
    if (nrow(R_mat) != length(vars_old) ||
        ncol(R_mat) != length(vars_old)) {
      return(NULL)
    }
    
    if (is.null(rownames(R_mat))) rownames(R_mat) <- vars_old
    if (is.null(colnames(R_mat))) colnames(R_mat) <- vars_old
    
    R_mat
  }
  
  for (d in names(obs.mean_new)) {
    
    if (!(d %in% names(obs.cov_new))) {
      warning("Date ", d, " exists in obs.mean but not obs.cov. Skipped.")
      next
    }
    
    site_names <- names(obs.mean_new[[d]])
    
    if (is.null(site_names) || length(site_names) == 0) {
      next
    }
    
    for (sid in site_names) {
      
      y <- obs.mean_new[[d]][[sid]]
      
      if (is.null(y)) next
      
      if (!is.data.frame(y)) {
        y <- as.data.frame(y)
      }
      
      vars_old <- names(y)
      
      if (length(vars_old) == 0) {
        obs.mean_new[[d]][[sid]] <- y
        obs.cov_new[[d]][[sid]]  <- matrix(0, 0, 0)
        next
      }
      
      vars_hit  <- intersect(vars_old, vars_drop)
      vars_keep <- setdiff(vars_old, vars_drop)
      
      if (length(vars_hit) == 0) {
        next
      }
      
      # ---- remove from obs.mean ----
      y_new <- y[, vars_keep, drop = FALSE]
      obs.mean_new[[d]][[sid]] <- y_new
      n_removed_mean <- n_removed_mean + length(vars_hit)
      
      # ---- remove matching rows / cols from obs.cov ----
      if (!(sid %in% names(obs.cov_new[[d]]))) {
        warning("Site ", sid, " at date ", d, " exists in obs.mean but not obs.cov. Skipped obs.cov update.")
        next
      }
      
      R_old <- obs.cov_new[[d]][[sid]]
      R_mat <- normalize_R(R_old, vars_old)
      
      if (is.null(R_mat)) {
        stop(
          "Cannot normalize obs.cov at date = ", d,
          ", site = ", sid,
          ". Check obs.mean / obs.cov dimensions."
        )
      }
      
      if (length(vars_keep) == 0) {
        R_new <- matrix(0, 0, 0)
      } else {
        R_new <- R_mat[vars_keep, vars_keep, drop = FALSE]
        dimnames(R_new) <- list(vars_keep, vars_keep)
      }
      
      obs.cov_new[[d]][[sid]] <- R_new
      n_removed_cov <- n_removed_cov + length(vars_hit)
    }
  }
  
  # ============================================================
  # Final alignment check
  # ============================================================
  
  bad <- list()
  k <- 1L
  
  for (d in names(obs.mean_new)) {
    
    if (!(d %in% names(obs.cov_new))) {
      bad[[k]] <- data.frame(
        time = d,
        site = NA_character_,
        issue = "date missing in obs.cov"
      )
      k <- k + 1L
      next
    }
    
    site_names <- names(obs.mean_new[[d]])
    if (is.null(site_names) || length(site_names) == 0) next
    
    for (sid in site_names) {
      
      y <- obs.mean_new[[d]][[sid]]
      
      if (is.null(y)) next
      if (!is.data.frame(y)) y <- as.data.frame(y)
      if (ncol(y) == 0) next
      
      if (!(sid %in% names(obs.cov_new[[d]]))) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = "site missing in obs.cov"
        )
        k <- k + 1L
        next
      }
      
      R <- as.matrix(obs.cov_new[[d]][[sid]])
      
      if (ncol(y) != nrow(R) || ncol(y) != ncol(R)) {
        bad[[k]] <- data.frame(
          time = d,
          site = sid,
          issue = paste0(
            "dim mismatch: y = ",
            ncol(y),
            ", R = ",
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
            "name mismatch: y = ",
            paste(names(y), collapse = ","),
            " | R = ",
            paste(colnames(R), collapse = ",")
          )
        )
        k <- k + 1L
      }
    }
  }
  
  if (length(bad) > 0) {
    bad_df <- do.call(rbind, bad)
    print(bad_df)
    stop("obs.mean_new and obs.cov_new are not aligned. Fix before saving.")
  }
  
  message("Removed variables from obs.mean entries: ", n_removed_mean)
  message("Removed corresponding rows/cols from obs.cov entries: ", n_removed_cov)
  message("Final obs.mean and obs.cov are aligned.")
  
  if (save_files) {
    if (is.null(out_mean_file) || is.null(out_cov_file)) {
      stop("If save_files = TRUE, provide both out_mean_file and out_cov_file.")
    }
    
    obs.mean <- obs.mean_new
    obs.cov  <- obs.cov_new
    
    save(obs.mean, file = out_mean_file)
    save(obs.cov,  file = out_cov_file)
    
    message("Saved obs.mean to: ", out_mean_file)
    message("Saved obs.cov to: ", out_cov_file)
  }
  
  if (assign_to_global) {
    assign("obs.mean", obs.mean_new, envir = .GlobalEnv)
    assign("obs.cov",  obs.cov_new,  envir = .GlobalEnv)
    message("Assigned obs.mean and obs.cov back to .GlobalEnv.")
  }
  
  list(
    obs.mean = obs.mean_new,
    obs.cov  = obs.cov_new
  )
}

# 假设你已经 load 了原始 obs.mean 和 obs.cov
# load("/projectnb/dietzelab/guYANG/pecan/output/obs/Rdata/obs.mean.monthly.with_flux.Rdata")
# load("/projectnb/dietzelab/guYANG/pecan/output/obs/Rdata/obs.cov.monthly.with_flux.Rdata")

res <- drop_state_vars_from_obs(
  vars_drop = c("AbvGrndWood"),
  save_files = TRUE,
  out_mean_file = "/projectnb/dietzelab/guYANG/pecan/output/obs/Rdata/obs.mean.monthly.noAGBSM.Rdata",
  out_cov_file  = "/projectnb/dietzelab/guYANG/pecan/output/obs/Rdata/obs.cov.monthly.noAGBSM.Rdata"
)
