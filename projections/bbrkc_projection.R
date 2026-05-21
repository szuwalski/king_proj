# ============================================================
# RKC forward projection with:
#   1) one No-GAM scenario
#   2) one GAM scenario per supplied projection temperature
#
# For now, each GAM scenario uses a single constant temperature
# over the entire projection period.
# ============================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------
# Helpers to read ADMB .rep and .par files
# ------------------------------------------------------------

read_rep_sections <- function(rep_file) {
  x <- readLines(rep_file, warn = FALSE)
  idx <- grep("^\\$", x)
  nm <- sub("^\\$", "", x[idx])
  
  out <- vector("list", length(idx))
  names(out) <- nm
  
  for (i in seq_along(idx)) {
    start <- idx[i] + 1
    end   <- if (i < length(idx)) idx[i + 1] - 1 else length(x)
    out[[i]] <- trimws(x[start:end])
  }
  
  out
}

extract_numeric <- function(lines) {
  if (length(lines) == 0) return(numeric(0))
  txt <- paste(lines, collapse = " ")
  as.numeric(scan(text = txt, quiet = TRUE))
}

extract_scalar <- function(lines) {
  vals <- extract_numeric(lines)
  if (length(vals) != 1) stop("Expected scalar, got length = ", length(vals))
  vals
}

extract_vector <- function(lines, expected_len = NULL) {
  vals <- extract_numeric(lines)
  if (!is.null(expected_len) && length(vals) != expected_len) {
    stop("Expected vector length ", expected_len, ", got ", length(vals))
  }
  vals
}

extract_matrix <- function(lines, nrow, ncol) {
  vals <- extract_numeric(lines)
  if (length(vals) != nrow * ncol) {
    stop("Expected matrix with ", nrow * ncol, " values, got ", length(vals))
  }
  matrix(vals, nrow = nrow, ncol = ncol, byrow = TRUE)
}

read_par_file <- function(par_file) {
  x <- readLines(par_file, warn = FALSE)
  
  hdr <- grep("^# ", x)
  section_names <- trimws(sub("^# ", "", x[hdr]))
  section_names <- sub(":$", "", section_names)
  
  out <- list()
  
  for (i in seq_along(hdr)) {
    nm <- section_names[i]
    
    if (grepl("^Number of parameters", nm)) next
    
    start <- hdr[i] + 1
    end   <- if (i < length(hdr)) hdr[i + 1] - 1 else length(x)
    vals  <- extract_numeric(x[start:end])
    
    out[[nm]] <- vals
  }
  
  out
}

# ------------------------------------------------------------
# Read ADMB .cor file
# ------------------------------------------------------------

read_cor_file <- function(cor_file) {
  x <- readLines(cor_file, warn = FALSE)
  
  keep <- grepl("^\\s*[0-9]+\\s+[A-Za-z0-9_]+\\s+[-+0-9.eE]+\\s+[-+0-9.eE]+", x)
  x <- x[keep]
  
  parse_one <- function(line) {
    toks <- strsplit(trimws(line), "\\s+")[[1]]
    data.frame(
      index = as.integer(toks[1]),
      name  = toks[2],
      value = as.numeric(toks[3]),
      sd    = as.numeric(toks[4]),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, lapply(x, parse_one))
}

# ------------------------------------------------------------
# Align a .cor series to a reference series from .rep
# ------------------------------------------------------------

align_cor_series <- function(cor_df, name, ref_values, ref_sd = NULL) {
  sub <- cor_df[cor_df$name == name, , drop = FALSE]
  
  if (nrow(sub) == 0) {
    stop("Series '", name, "' not found in .cor file.")
  }
  
  n_ref <- length(ref_values)
  n_cor <- nrow(sub)
  
  if (n_cor < n_ref) {
    stop("Series '", name, "' in .cor has fewer values than reference series.")
  }
  
  if (n_cor == n_ref) {
    out <- sub
    out$aligned_pos <- seq_len(n_ref)
    return(out)
  }
  
  starts <- seq_len(n_cor - n_ref + 1)
  rmse <- rep(NA_real_, length(starts))
  
  for (i in seq_along(starts)) {
    idx <- starts[i]:(starts[i] + n_ref - 1)
    rmse[i] <- sqrt(mean((sub$value[idx] - ref_values)^2, na.rm = TRUE))
  }
  
  best_start <- starts[which.min(rmse)]
  idx <- best_start:(best_start + n_ref - 1)
  
  out <- sub[idx, , drop = FALSE]
  out$aligned_pos <- seq_len(n_ref)
  out
}

# ------------------------------------------------------------
# Build model objects from rep + par + cor
# ------------------------------------------------------------

build_rkc_objects <- function(rep_file, par_file, cor_file = NULL) {
  rep_sec <- read_rep_sections(rep_file)
  par_obj <- read_par_file(par_file)
  
  styr   <- extract_scalar(rep_sec[["styr"]])
  endyr  <- extract_scalar(rep_sec[["endyr"]])
  years  <- styr:endyr
  n_year <- length(years)
  
  sizes  <- extract_vector(rep_sec[["sizes"]])
  n_size <- length(sizes)
  
  obj <- list(
    years                = years,
    styr                 = styr,
    endyr                = endyr,
    sizes                = sizes,
    n_year               = n_year,
    n_size               = n_size,
    n_obs                = extract_vector(rep_sec[["n_obs"]], n_year),
    surv_n_cv            = extract_vector(rep_sec[["surv_n_cv"]], n_year),
    n_obs_len            = extract_matrix(rep_sec[["n_obs_len"]], n_year, n_size),
    numbers_pred         = extract_vector(rep_sec[["numbers_pred"]], n_year),
    recruits_hist        = extract_vector(rep_sec[["recruits"]], n_year),
    M_mat_hist           = extract_matrix(rep_sec[["natural mortality"]], n_year, n_size),
    pop_num_hist         = extract_matrix(rep_sec[["pred_pop_num"]], n_year, n_size),
    pred_numbers_at_size = extract_matrix(rep_sec[["pred numbers at size"]], n_year, n_size),
    surv_sel_mat         = extract_matrix(rep_sec[["survey selectivity"]], n_year, n_size),
    fish_tot_sel         = extract_vector(rep_sec[["total_fish_sel"]], n_size),
    fish_ret_sel         = extract_vector(rep_sec[["ret_fish_sel"]], n_size),
    f_hist               = extract_vector(rep_sec[["est_fishing_mort"]], n_year),
    size_trans           = extract_matrix(rep_sec[["size_trans"]], n_size, n_size),
    in_prob_molt         = extract_vector(rep_sec[["in_prob_molt"]], n_size),
    prop_rec             = extract_vector(rep_sec[["temp_prop_rec"]], 3),
    par                  = par_obj
  )
  
  obj$total_abund_hist <- rowSums(obj$pop_num_hist)
  
  obs_n_at_size <- sweep(
    obj$n_obs_len,
    1,
    obj$n_obs,
    `*`
  )
  
  obj$mean_size_hist <- apply(
    obs_n_at_size,
    1,
    function(z) {
      if (all(is.na(z)) || sum(z, na.rm = TRUE) == 0) {
        NA_real_
      } else {
        weighted.mean(obj$sizes, w = z, na.rm = TRUE)
      }
    }
  )
  
  obj$M_hist <- obj$M_mat_hist[, 1]
  obj$p_mort_hist <- pmin(pmax(1 - exp(-obj$M_hist), 1e-6), 1 - 1e-6)
  
  obj$surv_sel_terminal <- obj$surv_sel_mat[n_year, ]
  
  obj$cor <- NULL
  if (!is.null(cor_file)) {
    cor_df <- read_cor_file(cor_file)
    
    cor_numbers_pred <- align_cor_series(
      cor_df = cor_df,
      name = "numbers_pred",
      ref_values = obj$numbers_pred
    )
    
    nat_sub <- cor_df[cor_df$name == "nat_m_dev", , drop = FALSE]
    if (nrow(nat_sub) < n_year) {
      stop("Series 'nat_m_dev' in .cor has fewer values than historical years.")
    }
    if (nrow(nat_sub) == n_year) {
      cor_nat_m <- nat_sub
      cor_nat_m$aligned_pos <- seq_len(n_year)
    } else {
      cor_nat_m <- nat_sub[seq_len(n_year), , drop = FALSE]
      cor_nat_m$aligned_pos <- seq_len(n_year)
    }
    
    obj$cor <- list(
      raw = cor_df,
      numbers_pred = cor_numbers_pred,
      nat_m_dev = cor_nat_m
    )
    
    obj$numbers_pred_sd_cor <- cor_numbers_pred$sd
    obj$M_hist_sd_cor <- cor_nat_m$sd
  }
  
  obj
}

# ------------------------------------------------------------
# Fit GAM to historical mortality
# ------------------------------------------------------------

fit_mortality_gam <- function(obj,
                              hist_temperature,
                              hist_ice,
                              k_abund = 4,
                              k_temp = 4,
                              k_size = 4,
                              k_ice = 4) {
  
  if (length(hist_temperature) != obj$n_year) {
    stop("hist_temperature must have length ", obj$n_year)
  }
  if (length(hist_ice) != obj$n_year) {
    stop("hist_ice must have length ", obj$n_year)
  }
  
  mort_sd <- if (!is.null(obj$M_hist_sd_cor)) obj$M_hist_sd_cor else rep(1, obj$n_year)
  
  gam_dat <- data.frame(
    year        = obj$years,
    p_mort      = obj$p_mort_hist,
    M           = obj$M_hist,
    M_sd_log    = mort_sd,
    Abundance   = obj$total_abund_hist / max(obj$total_abund_hist, na.rm = TRUE),
    Temperature = hist_temperature,
    Size        = obj$mean_size_hist,
    Ice         = hist_ice
  )
  
  gam_dat$p_mort <- pmin(pmax(gam_dat$p_mort, 1e-6), 1 - 1e-6)
  gam_dat$wts <- 1 / gam_dat$M_sd_log
  
  gam_fit <- mgcv::gam(
    p_mort ~
      s(Abundance,   k = k_abund) +
      s(Temperature, k = k_temp)  +
      s(Size,        k = k_size),
    data    = gam_dat,
    weights = wts,
    family  = betar(link = "logit"),
    method  = "REML"
  )
  
  mf <- model.frame(gam_fit)
  rows_used <- as.numeric(rownames(mf))
  
  gam_dat$fit_hist    <- NA_real_
  gam_dat$fit_hist_lo <- NA_real_
  gam_dat$fit_hist_hi <- NA_real_
  
  pred_link <- predict(gam_fit, type = "link", se.fit = TRUE)
  
  gam_dat$fit_hist[rows_used] <-
    plogis(pred_link$fit)
  
  gam_dat$fit_hist_lo[rows_used] <-
    plogis(pred_link$fit - 1.96 * pred_link$se.fit)
  
  gam_dat$fit_hist_hi[rows_used] <-
    plogis(pred_link$fit + 1.96 * pred_link$se.fit)
  
  list(data = gam_dat, fit = gam_fit)
}

# ------------------------------------------------------------
# One-step forward simulator from the ADMB template
# ------------------------------------------------------------

project_one_year <- function(state_n,
                             F_future,
                             M_future,
                             size_trans,
                             in_prob_molt,
                             fish_tot_sel,
                             fish_ret_sel,
                             prop_rec,
                             recruits_future,
                             discard_survival = 0) {
  
  temp_n <- state_n * exp(-0.17 * M_future)
  
  temp_catch_n <- temp_n * (1 - exp(-(F_future * fish_tot_sel)))
  temp_n_post  <- temp_n * exp(-(F_future * fish_tot_sel)) +
    temp_catch_n * (1 - fish_ret_sel) * discard_survival
  
  trans_n <- as.vector(size_trans %*% (temp_n_post * in_prob_molt))
  temp_n2 <- trans_n + temp_n_post * (1 - in_prob_molt)
  
  temp_n2[1:3] <- temp_n2[1:3] + recruits_future * prop_rec
  
  next_state <- temp_n2 * exp(-0.83 * M_future)
  
  list(
    next_state    = next_state,
    catch_numbers = temp_catch_n,
    postgrowth_n  = temp_n2
  )
}

# ------------------------------------------------------------
# Main projection engine
# ------------------------------------------------------------

run_rkc_projection <- function(rep_file,
                               par_file,
                               cor_file,
                               hist_temperature,
                               hist_ice,
                               gam_projection_temperatures,
                               future_ice,
                               future_F,
                               n_proj,
                               n_hist_ridges = 15,
                               n_sims = 1000,
                               recruitment_start_year = 1976,
                               discard_survival = 0,
                               seed = 123) {
  
  set.seed(seed)
  
  obj <- build_rkc_objects(rep_file, par_file, cor_file = cor_file)
  
  gam_obj <- fit_mortality_gam(
    obj = obj,
    hist_temperature = hist_temperature,
    hist_ice = hist_ice
  )
  
  if (length(future_ice) != n_proj) {
    stop("future_ice must have length = n_proj")
  }
  
  if (length(gam_projection_temperatures) < 1) {
    stop("gam_projection_temperatures must have at least one value.")
  }
  
  proj_years <- (obj$endyr + 1):(obj$endyr + n_proj)
  
  rec_pool_idx <- which(obj$years >= recruitment_start_year)
  if (length(rec_pool_idx) == 0) {
    stop("No recruitment years >= ", recruitment_start_year, " found in model years.")
  }
  
  rec_pool <- obj$recruits_hist[rec_pool_idx]
  m_pool <- obj$M_hist
  init_state <- obj$pop_num_hist[obj$n_year, ]
  
  gam_scenario_names <- paste0("GAM_", gam_projection_temperatures)
  scenario_names <- c("No GAM", gam_scenario_names)
  
  sim_out <- vector("list", length(scenario_names))
  names(sim_out) <- scenario_names
  
  for (scen in scenario_names) {
    
    survey_store <- matrix(NA_real_, nrow = n_sims, ncol = n_proj)
    M_store      <- matrix(NA_real_, nrow = n_sims, ncol = n_proj)
    R_store      <- matrix(NA_real_, nrow = n_sims, ncol = n_proj)
    
    nsize_store <- array(
      NA_real_,
      dim = c(n_sims, n_proj, obj$n_size),
      dimnames = list(NULL, proj_years, obj$sizes)
    )
    
    if (scen != "No GAM") {
      scen_temp <- as.numeric(sub("^GAM_", "", scen))
      scenario_temp_vec <- rep(scen_temp, n_proj)
    }
    
    for (s in seq_len(n_sims)) {
      state <- init_state
      
      for (tt in seq_len(n_proj)) {
        
        R_t <- sample(rec_pool, size = 1, replace = TRUE)
        
        survey_state <- state * obj$surv_sel_terminal
        abund_t <- sum(state)
        size_t  <- weighted.mean(obj$sizes, w = survey_state)
        
        if (scen == "No GAM") {
          M_t <- sample(m_pool, size = 1, replace = TRUE)
        } else {
          new_dat <- data.frame(
            Abundance   = abund_t / max(obj$total_abund_hist, na.rm = TRUE),
            Temperature = scenario_temp_vec[tt],
            Size        = size_t
          )
          
          pred_link <- predict(
            gam_obj$fit,
            newdata = new_dat,
            type = "link",
            se.fit = TRUE
          )
          
          eta_draw <- rnorm(
            n    = 1,
            mean = as.numeric(pred_link$fit),
            sd   = as.numeric(pred_link$se.fit)
          )
          
          p_draw <- plogis(eta_draw)
          p_draw <- pmin(pmax(p_draw, 1e-6), 1 - 1e-6)
          M_t <- -log(1 - p_draw)
        }
        
        step <- project_one_year(
          state_n          = state,
          F_future         = future_F,
          M_future         = M_t,
          size_trans       = obj$size_trans,
          in_prob_molt     = obj$in_prob_molt,
          fish_tot_sel     = obj$fish_tot_sel,
          fish_ret_sel     = obj$fish_ret_sel,
          prop_rec         = obj$prop_rec,
          recruits_future  = R_t,
          discard_survival = discard_survival
        )
        
        state <- step$next_state
        survey_t <- sum(state * obj$surv_sel_terminal)
        
        survey_store[s, tt]  <- survey_t
        M_store[s, tt]       <- M_t
        R_store[s, tt]       <- R_t
        nsize_store[s, tt, ] <- state
      }
    }
    
    sim_out[[scen]] <- list(
      survey = survey_store,
      M      = M_store,
      R      = R_store,
      nsize  = nsize_store,
      projection_temperature = if (scen == "No GAM") NA_real_ else scen_temp
    )
  }
  
  out <- list(
    obj        = obj,
    gam        = gam_obj,
    sims       = sim_out,
    proj_years = proj_years,
    inputs     = list(
      hist_temperature           = hist_temperature,
      hist_ice                   = hist_ice,
      gam_projection_temperatures = gam_projection_temperatures,
      future_ice                 = future_ice,
      future_F                   = future_F,
      n_proj                     = n_proj,
      n_hist_ridges              = n_hist_ridges,
      n_sims                     = n_sims,
      recruitment_start_year     = recruitment_start_year,
      discard_survival           = discard_survival
    )
  )
  
  class(out) <- "rkc_projection"
  out
}

# ------------------------------------------------------------
# Summaries for plotting
# ------------------------------------------------------------

summarize_projection <- function(proj_obj,
                                 probs = c(0.025, 0.5, 0.975),
                                 n_hist_ridges = 15) {
  
  obj <- proj_obj$obj
  proj_years <- proj_obj$proj_years
  
  y_max_cap <- 1.2 * max(obj$n_obs, na.rm = TRUE)
  
  hist_df <- data.frame(
    year        = obj$years,
    observed    = obj$n_obs,
    observed_cv = obj$surv_n_cv,
    estimate    = obj$numbers_pred
  ) %>%
    mutate(
      log_sd_obs = sqrt(log(1 + observed_cv^2)),
      obs_lo = observed * exp(-1.96 * log_sd_obs),
      obs_hi = observed * exp( 1.96 * log_sd_obs),
      obs_lo = pmax(obs_lo, 0),
      obs_hi = pmin(obs_hi, y_max_cap)
    )
  
  if (!is.null(obj$numbers_pred_sd_cor)) {
    hist_df$est_lo <- pmax(obj$numbers_pred - 1.96 * obj$numbers_pred_sd_cor, 0)
    hist_df$est_hi <- pmin(obj$numbers_pred + 1.96 * obj$numbers_pred_sd_cor, y_max_cap)
  } else {
    hist_df$est_lo <- NA_real_
    hist_df$est_hi <- NA_real_
  }
  
  proj_df <- bind_rows(lapply(names(proj_obj$sims), function(scen) {
    x <- proj_obj$sims[[scen]]$survey
    data.frame(
      year     = proj_years,
      scenario = scen,
      lo       = apply(x, 2, quantile, probs = probs[1], na.rm = TRUE),
      med      = apply(x, 2, quantile, probs = probs[2], na.rm = TRUE),
      hi       = apply(x, 2, quantile, probs = probs[3], na.rm = TRUE)
    )
  }))
  
  hist_mort_df <- data.frame(
    year   = obj$years,
    M      = obj$M_hist,
    p_mort = obj$p_mort_hist
  )
  
  if (!is.null(obj$M_hist_sd_cor)) {
    hist_mort_df$M_lo <- hist_mort_df$M * exp(-1.96 * obj$M_hist_sd_cor)
    hist_mort_df$M_hi <- hist_mort_df$M * exp( 1.96 * obj$M_hist_sd_cor)
    hist_mort_df$p_lo <- 1 - exp(-hist_mort_df$M_lo)
    hist_mort_df$p_hi <- 1 - exp(-hist_mort_df$M_hi)
  } else {
    hist_mort_df$p_lo <- NA_real_
    hist_mort_df$p_hi <- NA_real_
  }
  
  hist_gam_fit_df <- proj_obj$gam$data %>%
    transmute(
      year   = year,
      fit    = fit_hist,
      fit_lo = fit_hist_lo,
      fit_hi = fit_hist_hi
    )
  
  proj_mort_df <- bind_rows(lapply(names(proj_obj$sims), function(scen) {
    m_mat <- proj_obj$sims[[scen]]$M
    p_mat <- 1 - exp(-m_mat)
    data.frame(
      year     = proj_years,
      scenario = scen,
      lo       = apply(p_mat, 2, quantile, probs = probs[1], na.rm = TRUE),
      med      = apply(p_mat, 2, quantile, probs = probs[2], na.rm = TRUE),
      hi       = apply(p_mat, 2, quantile, probs = probs[3], na.rm = TRUE)
    )
  }))
  
  hist_keep <- tail(seq_len(obj$n_year), n_hist_ridges)
  
  hist_n_at_size <- sweep(
    obj$n_obs_len[hist_keep, , drop = FALSE],
    1,
    obj$n_obs[hist_keep],
    `*`
  )
  
  ridge_hist <- as.data.frame(hist_n_at_size)
  colnames(ridge_hist) <- paste0("size_", seq_len(obj$n_size))
  
  ridge_hist <- ridge_hist %>%
    mutate(year = obj$years[hist_keep], scenario = "Historical") %>%
    pivot_longer(
      cols      = starts_with("size_"),
      names_to  = "size_bin",
      values_to = "numbers"
    ) %>%
    mutate(
      size_index = as.integer(sub("size_", "", size_bin)),
      size = obj$sizes[size_index]
    ) %>%
    select(year, scenario, size, numbers)
  
  ridge_proj <- bind_rows(lapply(names(proj_obj$sims), function(scen) {
    
    arr    <- proj_obj$sims[[scen]]$nsize
    n_sims <- dim(arr)[1]
    n_proj <- dim(arr)[2]
    n_size <- dim(arr)[3]
    
    surv_sel <- obj$surv_sel_terminal
    
    survey_nsize_arr <- array(
      NA_real_,
      dim = c(n_sims, n_proj, n_size)
    )
    
    for (ss in seq_len(n_sims)) {
      for (tt in seq_len(n_proj)) {
        survey_nsize_arr[ss, tt, ] <- arr[ss, tt, ] * surv_sel
      }
    }
    
    med_mat <- apply(survey_nsize_arr, c(2, 3), median, na.rm = TRUE)
    
    out <- as.data.frame(med_mat)
    colnames(out) <- paste0("size_", seq_len(obj$n_size))
    
    out %>%
      mutate(year = proj_years, scenario = scen) %>%
      pivot_longer(
        cols      = starts_with("size_"),
        names_to  = "size_bin",
        values_to = "numbers"
      ) %>%
      mutate(
        size_index = as.integer(sub("size_", "", size_bin)),
        size = obj$sizes[size_index]
      ) %>%
      select(year, scenario, size, numbers)
  }))
  
  list(
    hist_df          = hist_df,
    proj_df          = proj_df,
    hist_mort_df     = hist_mort_df,
    hist_gam_fit_df  = hist_gam_fit_df,
    proj_mort_df     = proj_mort_df,
    ridge_hist       = ridge_hist,
    ridge_proj       = ridge_proj,
    y_max_cap        = y_max_cap
  )
}

# ------------------------------------------------------------
# Plot method
# ------------------------------------------------------------

plot_rkc_projection <- function(proj_obj,
                                n_hist_ridges = 15,
                                ridge_scale = 9,
                                smooth_height_ratio = 0.25,
                                mort_height_ratio = 1/3,
                                colors = NULL,
                                hist_color = "grey55",
                                legend_position = c(0.50, 0.86)) {
  
  s <- summarize_projection(proj_obj, n_hist_ridges = n_hist_ridges)
  proj_start_year <- min(proj_obj$proj_years)
  
  scen_names <- names(proj_obj$sims)
  
  if (is.null(colors)) {
    default_cols <- c(
      "No GAM" = "#D55E00",
      "GAM_2"  = "#0072B2",
      "GAM_4"  = "#009E73",
      "GAM_6"  = "#CC79A7",
      "GAM_8"  = "#56B4E9"
    )
    
    if (all(scen_names %in% names(default_cols))) {
      colors <- default_cols[scen_names]
    } else {
      extra_cols <- grDevices::hcl.colors(length(scen_names), palette = "Dark 3")
      colors <- stats::setNames(extra_cols, scen_names)
      if ("No GAM" %in% scen_names) colors["No GAM"] <- "#D55E00"
    }
  }
  
  build_gam_smooth_row <- function(gam_fit) {
    
    dev_expl_pct <- round(summary(gam_fit)$dev.expl * 100)
    
    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    
    plotted <- plot(
      gam_fit,
      pages = 1,
      se = TRUE,
      shade = FALSE,
      scale = 0
    )
    
    smooth_plots <- vector("list", length(plotted))
    
    for (i in seq_along(plotted)) {
      pd <- plotted[[i]]
      
      tmp <- data.frame(
        x    = pd$x,
        fit  = pd$fit,
        y_up = pd$fit + 2 * pd$se,
        y_dn = pd$fit - 2 * pd$se
      )
      
      p <- ggplot(tmp, aes(x = x, y = fit)) +
        geom_ribbon(
          aes(ymin = y_dn, ymax = y_up),
          fill = "grey80",
          alpha = 0.35
        ) +
        geom_line(linewidth = 0.75, color = "black") +
        geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
        theme_bw(base_size = 8) +
        labs(
          x = pd$xlab,
          y = NULL
        ) +
        theme(
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(face = "bold", size = 8),
          plot.margin = margin(2, 4, 2, 4)
        )
      
      if (grepl("Abundance", pd$xlab, fixed = TRUE)) {
        x_rng <- range(tmp$x, na.rm = TRUE)
        y_rng <- range(c(tmp$y_dn, tmp$y_up), na.rm = TRUE)
        
        p <- p +
          annotate(
            "text",
            x = x_rng[1] + 0.03 * diff(x_rng),
            y = y_rng[2] - 0.08 * diff(y_rng),
            label = paste0(dev_expl_pct, "%"),
            hjust = 0,
            vjust = 1,
            size = 2.8,
            color = "grey30"
          )
      }
      
      smooth_plots[[i]] <- p
    }
    
    patchwork::wrap_plots(smooth_plots, nrow = 1)
  }
  
  survey_y_top <- s$y_max_cap
  mort_y_top   <- 1
  
  p_mort <- ggplot() +
    geom_vline(xintercept = proj_start_year - 0.5, linetype = 3, color = "grey40") +
    geom_errorbar(
      data = s$hist_mort_df,
      aes(x = year, ymin = p_lo, ymax = p_hi),
      width = 0,
      color = "black",
      alpha = 0.55
    ) +
    geom_point(
      data = s$hist_mort_df,
      aes(x = year, y = p_mort),
      color = "black",
      size = 1.4
    ) +
    geom_ribbon(
      data = s$hist_gam_fit_df,
      aes(x = year, ymin = fit_lo, ymax = fit_hi),
      fill = "grey80",
      alpha = 0.35
    ) +
    geom_line(
      data = s$hist_gam_fit_df,
      aes(x = year, y = fit),
      color = "grey45",
      linewidth = 0.95
    ) +
    geom_ribbon(
      data = s$proj_mort_df,
      aes(x = year, ymin = lo, ymax = hi, fill = scenario),
      alpha = 0.18
    ) +
    geom_line(
      data = s$proj_mort_df,
      aes(x = year, y = med, color = scenario),
      linewidth = 1
    ) +
    annotate(
      "text",
      x = min(s$hist_mort_df$year, na.rm = TRUE) + 2,
      y = mort_y_top * 0.93,
      label = "Historical",
      hjust = 0,
      vjust = 1,
      size = 3,
      color = "grey25"
    ) +
    annotate(
      "text",
      x = proj_start_year + 1,
      y = mort_y_top * 0.93,
      label = "Projection",
      hjust = 0,
      vjust = 1,
      size = 3,
      color = "grey25"
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = NULL,
      y = "p(mortality)"
    )
  
  p_survey <- ggplot() +
    geom_vline(xintercept = proj_start_year - 0.5, linetype = 3, color = "grey40") +
    geom_errorbar(
      data = s$hist_df,
      aes(x = year, ymin = obs_lo, ymax = obs_hi),
      width = 0,
      color = "black",
      alpha = 0.55
    ) +
    geom_point(
      data = s$hist_df,
      aes(x = year, y = observed),
      color = "black",
      size = 1.7
    ) +
    geom_ribbon(
      data = s$hist_df,
      aes(x = year, ymin = est_lo, ymax = est_hi),
      fill = "grey82",
      alpha = 0.40
    ) +
    geom_line(
      data = s$hist_df,
      aes(x = year, y = estimate),
      color = "grey45",
      linewidth = 0.95
    ) +
    geom_ribbon(
      data = s$proj_df,
      aes(x = year, ymin = lo, ymax = hi, fill = scenario),
      alpha = 0.18
    ) +
    geom_line(
      data = s$proj_df,
      aes(x = year, y = med, color = scenario),
      linewidth = 1.05
    ) +
    annotate(
      "text",
      x = min(s$hist_df$year, na.rm = TRUE) + 2,
      y = survey_y_top * 0.95,
      label = "Historical",
      hjust = 0,
      vjust = 1,
      size = 3,
      color = "grey25"
    ) +
    annotate(
      "text",
      x = proj_start_year + 1,
      y = survey_y_top * 0.95,
      label = "Projection",
      hjust = 0,
      vjust = 1,
      size = 3,
      color = "grey25"
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    coord_cartesian(ylim = c(0, s$y_max_cap)) +
    theme_bw() +
    theme(
      legend.position = legend_position,
      legend.background = element_rect(
        fill = alpha("white", 0.75),
        color = NA
      ),
      legend.key = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Year",
      y = "Survey numbers",
      color = NULL,
      fill = NULL
    )
  
  all_numbers <- c(s$ridge_hist$numbers, s$ridge_proj$numbers)
  ridge_denom <- max(all_numbers, na.rm = TRUE)
  
  ridge_hist <- s$ridge_hist %>%
    mutate(height = numbers / ridge_denom)
  
  ridge_proj <- s$ridge_proj %>%
    mutate(height = numbers / ridge_denom)
  
  p_ridge <- ggplot() +
    geom_ridgeline(
      data = ridge_hist,
      aes(x = size, y = year, height = height, group = year),
      fill = "grey75",
      color = "grey45",
      alpha = 0.60,
      scale = ridge_scale,
      size = 0.25
    )
  
  for (scen in scen_names) {
    if (scen == "Historical") next
    
    p_ridge <- p_ridge +
      geom_ridgeline(
        data = subset(ridge_proj, scenario == scen),
        aes(x = size, y = year, height = height, group = interaction(year, scenario)),
        fill = alpha(colors[scen], 0.14),
        color = colors[scen],
        alpha = 0.30,
        scale = ridge_scale,
        size = 0.42
      )
  }
  
  p_ridge <- p_ridge +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Size",
      y = "Year"
    )
  
  p_top <- build_gam_smooth_row(proj_obj$gam$fit)
  
  survey_units <- 1 / mort_height_ratio
  p_left <- patchwork::wrap_plots(
    p_mort,
    p_survey,
    ncol = 1,
    heights = c(1, survey_units)
  )
  
  p_bottom <- patchwork::wrap_plots(
    p_left,
    p_ridge,
    ncol = 2,
    widths = c(1.05, 1)
  )
  
  bottom_units <- 1 / smooth_height_ratio
  
  patchwork::wrap_plots(
    p_top,
    p_bottom,
    ncol = 1,
    heights = c(1, bottom_units)
  )
}

# ------------------------------------------------------------
# Convenience wrapper (UPDATED)
# ------------------------------------------------------------

run_and_plot_rkc_projection <- function(rep_file,
                                        par_file,
                                        cor_file,
                                        hist_temperature,
                                        hist_ice,
                                        future_ice,
                                        future_F,
                                        n_proj,
                                        gam_temperatures,   # <- vector of temps (e.g., c(3,5,7))
                                        n_hist_ridges = 15,
                                        n_sims = 1000,
                                        recruitment_start_year = 1976,
                                        discard_survival = 0,
                                        seed = 123) {
  
  # ---- build projection ----
  proj <- run_rkc_projection(
    rep_file = rep_file,
    par_file = par_file,
    cor_file = cor_file,
    hist_temperature = hist_temperature,
    hist_ice = hist_ice,
    future_ice = future_ice,
    future_F = future_F,
    n_proj = n_proj,
    gam_projection_temperatures = gam_temperatures,  # <- key input
    n_hist_ridges = n_hist_ridges,
    n_sims = n_sims,
    recruitment_start_year = recruitment_start_year,
    discard_survival = discard_survival,
    seed = seed
  )
  
  # ---- plot ----
  fig <- plot_rkc_projection(
    proj_obj = proj,
    n_hist_ridges = n_hist_ridges
  )
  
  return(list(
    projection = proj,
    figure = fig
  ))
}

rep_file <- "ADMB/bbrkc/test/rkc.rep"
par_file <- "ADMB/bbrkc/test/rkc.par"
cor_file <- "ADMB/bbrkc/test/rkc.cor"

obj <- build_rkc_objects(rep_file, par_file, cor_file = cor_file)

in_ice  <- read.csv("data/ice_extent.csv")
in_temp <- read.csv("data/alt_metrics_calc.csv")

hist_temperature <- dplyr::filter(in_temp, stock == "BBRKC")$Temperature
hist_temperature <- c(hist_temperature[1:45], NA, hist_temperature[46:length(hist_temperature)])
hist_ice <- in_ice$Ice


n_proj <- 10
future_temperature <- rep(5.5, n_proj)
future_ice         <- rep(20, n_proj)


out <- run_and_plot_rkc_projection(
  rep_file = rep_file,
  par_file = par_file,
  cor_file = cor_file,
  hist_temperature = hist_temperature,
  hist_ice = hist_ice,
  future_ice = rep(20, 10),
  future_F = 0.25,
  n_proj = 10,
  gam_temperatures = c(1, 4, 5.5),  # Cold / Avg / Hot (or whatever you want)
  n_sims = 500
)

print(out$figure)

 ggsave("plots/rkc_projection.png", out$figure, width = 12, height = 9, dpi = 300)
 