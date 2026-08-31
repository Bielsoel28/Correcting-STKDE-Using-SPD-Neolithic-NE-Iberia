##########################################################################
# This script is Main for Radiocarbon SPD Analysis                      #
#                                                                       #   
#                                                                       #                                                   
#                                                                       #                                                   
#                                                                       #                                                   
# Author: Biel Soriano Elias & Andreu Monforte-Barberan                 #
# Affiliation : Autonomous University of Barcelona                      #
# Creation date : XX/XX/XXXX                                            #
# E-mail: biel.soriano@uab.cat & andreumonbar@gmail.com                 #
##########################################################################

# 0 ENVIRONMENT SETUP ##########################################################

## 0.1 Prepare environment ======================================================

# Clean workspace
rm(list = ls())

# Output folders (create if they do not exist)
dir.create("Results/1_SPD", showWarnings = FALSE, recursive = TRUE)

## 0.2 Install packages =========================================================

# Required packages
packages <- c(
  "readxl","terra","sf","dplyr","tidyr","foreach","doParallel",
  "ggplot2","spatstat","stars","viridis","rcarbon", "magick","corrplot", "mclust",
  "vegan", "Metrics", "scales","gdata","data.table","Bchron","BBmisc","sqldf", "reshape2"
)

# Install packages if missing
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

## 0.3 Show session information =================================================

sessionInfo()

# 1 SPD ANALYSIS ###############################################################

## 1.0 Load data ================================================================

# Raw radiocarbon dataset
c14_raw <- read_excel("Data/Raw_burials/c14_raw_burials.xlsx")

# Ensure numeric coordinates
c14_raw$LOG <- as.numeric(c14_raw$LOG)
c14_raw$LAT <- as.numeric(c14_raw$LAT)

# Convert to spatial object
coords_sf <- st_as_sf(
  c14_raw,
  coords = c("LOG","LAT"),
  crs = 4326
)

# Transform to UTM (ETRS89 / UTM zone 31N)
coords_utm <- st_transform(coords_sf, 25831)
c14_sf <- coords_utm #make a copy for latter analysis

utm <- st_coordinates(coords_utm)

c14_raw$X <- utm[,1]
c14_raw$Y <- utm[,2]

### 1.0.1 Standard partition ===================================================

### This block is the standard procedure of the paper, if holdout validation is aimed skip to 1.0.2 ###
### If date & site reconstruction is aime, run this part and subsequently parts 1.0.3.1 & 2 ###

# Radiocarbon calibration
C14_raw_calibration <- rcarbon::calibrate(
  x = c14_raw$DATE,
  errors = c14_raw$SD,
  ids = c14_raw$ID
)

# Raw no radiocarbon dataset
no_c14_raw <- read_excel("Data/Raw_burials/no_c14_raw_burials.xlsx")

#Erase NAs
no_c14_raw <- na.omit(no_c14_raw)

# Ensure numeric coordinates
no_c14_raw$LOG <- as.numeric(no_c14_raw$LOG)
no_c14_raw$LAT <- as.numeric(no_c14_raw$LAT)

# Convert to spatial object
coords_sf <- st_as_sf(
  no_c14_raw,
  coords = c("LOG","LAT"),
  crs = 4326
)

# Transform to UTM (ETRS89 / UTM zone 31N)
coords_utm <- st_transform(coords_sf, 25831)
no_c14_sf <- coords_utm #make a copy for latter analysis

utm <- st_coordinates(coords_utm)

no_c14_raw$X <- utm[,1]
no_c14_raw$Y <- utm[,2]

## Compute number of site for weight calculation
# Combine both datasets for counting
combined_data <- bind_rows(c14_raw, no_c14_raw)

# Correction for site overrepresentation using BOTH datasets
site_counts <- combined_data %>%
  group_by(SITE) %>%
  summarise(n_dates = n(), .groups = "drop")

# Join counts back only
c14_raw <- c14_raw %>%
  left_join(site_counts, by = "SITE")

no_c14_raw <- no_c14_raw %>%
  left_join(site_counts, by = "SITE")

### 1.0.2 (Optional) Holdout validation ========================================

### This section must only be run if the holdout validation is intended ###
### To do so, instead of running 1.0.1, run this part (1.0.2) ###

# HOLDOUT SPLIT: fix 25% of dated rows as "treated as undated"
HOLDOUT_FRACTION <- 0.25
set.seed(42)  # fixed seed guarantees the same subset every run

#Subset dates 
holdout_idx <- sample(
  seq_len(nrow(c14_raw)),
  size = floor(HOLDOUT_FRACTION * nrow(c14_raw)),
  replace = FALSE
)

# Rows kept in the dated pipeline
c14_dated   <- c14_raw[-holdout_idx, ]

# Rows moved to the undated pipeline (dates are NOT used)
c14_holdout <- c14_raw[holdout_idx, ]

# Raw no radiocarbon dataset
no_c14_raw <- read_excel("Data/Raw_burials/no_c14_raw_burials.xlsx")

#Erase NAs
no_c14_raw <- na.omit(no_c14_raw)

# Ensure numeric coordinates
no_c14_raw$LOG <- as.numeric(no_c14_raw$LOG)
no_c14_raw$LAT <- as.numeric(no_c14_raw$LAT)

# Convert to spatial object
coords_sf <- st_as_sf(
  no_c14_raw,
  coords = c("LOG","LAT"),
  crs = 4326
)

# Transform to UTM (ETRS89 / UTM zone 31N)
coords_utm <- st_transform(coords_sf, 25831)
no_c14_sf <- coords_utm #make a copy for latter analysis

utm <- st_coordinates(coords_utm)

no_c14_raw$X <- utm[,1]
no_c14_raw$Y <- utm[,2]

# Append holdout rows to the undated pool (dates ignored from here on)
no_c14_raw <- bind_rows(no_c14_raw, c14_holdout)

## Compute number of site for weight calculation
# Combine both datasets for counting
combined_data <- bind_rows(c14_dated, no_c14_raw)

# Correction for site overrepresentation using BOTH datasets
site_counts <- combined_data %>%
  group_by(SITE) %>%
  summarise(n_dates = n(), .groups = "drop")

# Join counts back only
c14_dated <- c14_dated %>%
  left_join(site_counts, by = "SITE")

no_c14_raw <- no_c14_raw %>%
  left_join(site_counts, by = "SITE")

# Keep c14_raw name for legacy compatibility 
c14_raw <- c14_dated

# Radiocarbon calibration (dated subset only)
C14_raw_calibration <- rcarbon::calibrate(
  x = c14_raw$DATE,
  errors = c14_raw$SD,
  ids = c14_raw$ID
)

### 1.0.3 (Optional) Date & site reconstruction validation =====================

### These section must only be run if the site & reconstruction validation is intended ###
### To do so, once you ran 1.0.1, run any of these two sections (1.0.3. 1 & 2) ###

#### 1.0.3.1 Date reconstruction ===============================================

material_cols <- c("VAR","OBS","SIL BEDU", "AXE EX", "BRACE", "MONT", "BQ", "CH")

time_range  <- c(6400, 5000)   # calBP range for spd(), adjust to your data span
n_random    <- 1000         # number of random SPD comparisons
n_resample  <- 5000         # sample size drawn from each density for Wilcoxon
weighting   <- "equal"      # "equal" | "count" -> how material SPDs are combined

set.seed(1)                 # for reproducibility of the random null model

##### Helper functions =========================================================

# row index in c14_raw / C14_raw_calibration for a given date ID
get_idx <- function(id) match(id, c14_raw$ID)

# extract a normalised (sums to 1) calBP/PrDens data.frame for one date
get_grid <- function(idx) {
  g <- C14_raw_calibration$grids[[as.character(idx)]]
  g$PrDens <- g$PrDens / sum(g$PrDens)
  g
}

# align two grids on the union of calBP, filling missing years with 0, and re-normalise both to sum to 1
align_grids <- function(g1, g2) {
  bp <- sort(union(g1$calBP, g2$calBP), decreasing = TRUE)
  p1 <- g1$PrDens[match(bp, g1$calBP)]; p1[is.na(p1)] <- 0
  p2 <- g2$PrDens[match(bp, g2$calBP)]; p2[is.na(p2)] <- 0
  list(calBP = bp, p1 = p1 / sum(p1), p2 = p2 / sum(p2))
}

hellinger_distance <- function(p, q) sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)

wasserstein_distance <- function(bp, p, q) {
  ord <- order(bp)
  bp <- bp[ord]; p <- p[ord]; q <- q[ord]
  cdf_diff <- abs(cumsum(p) - cumsum(q))
  w <- abs(diff(bp))
  sum(cdf_diff[-length(bp)] * w)
}

wilcoxon_compare <- function(bp, p, q, n = n_resample) {
  s1 <- sample(bp, n, replace = TRUE, prob = p)
  s2 <- sample(bp, n, replace = TRUE, prob = q)
  wt <- suppressWarnings(wilcox.test(s1, s2))
  c(W = unname(wt$statistic), p_value = wt$p.value)
}

compare_distributions <- function(g1, g2) {
  al <- align_grids(g1, g2)
  hd <- hellinger_distance(al$p1, al$p2)
  wd <- wasserstein_distance(al$calBP, al$p1, al$p2)
  wc <- wilcoxon_compare(al$calBP, al$p1, al$p2)
  c(hellinger   = hd,
    wasserstein = wd,
    wilcoxon_W  = as.numeric(wc["W"]),
    wilcoxon_p  = as.numeric(wc["p_value"]))
}

# build a normalised SPD (calBP/PrDens data.frame) from a set of row indices
build_spd <- function(idx, timeRange = time_range) {
  if (length(idx) == 0) return(NULL)
  sub <- C14_raw_calibration[idx]
  s <- spd(sub, timeRange = timeRange, spdnormalised = TRUE, verbose = FALSE)
  s$grid
}

# combine several material-specific SPDs into one reconstructed distribution
combine_spds <- function(spd_list, weights) {
  bp_all <- sort(unique(unlist(lapply(spd_list, `[[`, "calBP"))), decreasing = TRUE)
  mat <- sapply(seq_along(spd_list), function(i) {
    g <- spd_list[[i]]
    v <- g$PrDens[match(bp_all, g$calBP)]
    v[is.na(v)] <- 0
    v * weights[i]
  })
  combined <- rowSums(mat)
  combined <- combined / sum(combined)
  data.frame(calBP = bp_all, PrDens = combined)
}

# function to reconstruct date
reconstruct_date <- function(target_id, weighting = "equal") {
  
  target_idx <- get_idx(target_id)
  if (is.na(target_idx)) stop("target_id not found in c14_raw$ID")
  
  # materials present at the target date
  date_row <- c14_raw[c14_raw$ID == target_id, material_cols]
  present  <- material_cols[as.numeric(date_row) > 0]
  if (length(present) == 0) stop("target date has no material associations")
  
  # pool of candidate dates = everything except the target
  pool <- c14_raw[c14_raw$ID != target_id, ]
  
  spd_list <- list()
  weights  <- c()
  used_idx <- c()
  
  for (m in present) {
    ids_m <- pool$ID[pool[[m]] > 0]
    idx_m <- match(ids_m, c14_raw$ID)
    idx_m <- idx_m[!is.na(idx_m)]
    if (length(idx_m) == 0) next
    used_idx <- union(used_idx, idx_m)
    spd_m <- build_spd(idx_m)
    if (is.null(spd_m)) next
    spd_list[[m]] <- spd_m
    weights[m] <- if (weighting == "count") sum(date_row[[m]]) else 1
  }
  
  if (length(spd_list) == 0) stop("no material-linked dates found for this date")
  
  weights <- weights / sum(weights)
  reconstructed <- combine_spds(spd_list, weights)
  original      <- get_grid(target_idx)
  
  metrics <- compare_distributions(reconstructed, original)
  
  list(
    target_id     = target_id,
    materials     = present,
    n_dates_used  = length(used_idx),
    used_idx      = used_idx,
    reconstructed = reconstructed,
    original      = original,
    metrics       = metrics
  )
}

# function of random model
random_null <- function(target_id, n_used, n_iter = n_random) {
  
  target_idx <- get_idx(target_id)
  original   <- get_grid(target_idx)
  pool_idx   <- setdiff(seq_len(nrow(c14_raw)), target_idx)
  
  res <- matrix(NA_real_, nrow = n_iter, ncol = 4,
                dimnames = list(NULL, c("hellinger", "wasserstein",
                                        "wilcoxon_W", "wilcoxon_p")))
  
  for (i in seq_len(n_iter)) {
    samp_idx <- sample(pool_idx, size = min(n_used, length(pool_idx)), replace = FALSE)
    rand_spd <- build_spd(samp_idx)
    if (is.null(rand_spd)) next
    res[i, ] <- compare_distributions(rand_spd, original)
  }
  
  res
}

# function to plot: original vs reconstructed distribution 
plot_reconstruction <- function(recon, results_row = NULL) {
  
  orig  <- recon$original
  rec   <- recon$reconstructed
  
  xlim_range <- range(c(orig$calBP, rec$calBP))
  
  op <- par(mar = c(4.5, 4.5, 3, 1))
  on.exit(par(op))
  
  plot(orig$calBP, orig$PrDens, type = "n",
       xlim = rev(xlim_range),                 # calBP: older (larger) on the left
       ylim = c(0, max(orig$PrDens, rec$PrDens) * 1.1),
       xlab = "cal BP", ylab = "Probability density",
       main = paste0("date ", recon$target_id,
                     " - observed vs. reconstructed date"))
  
  # observed (true) calibrated date, filled
  polygon(c(orig$calBP, rev(orig$calBP)),
          c(orig$PrDens, rep(0, length(orig$PrDens))),
          col = adjustcolor("grey40", alpha.f = 0.5), border = NA)
  
  # reconstructed distribution, outline only
  lines(rec$calBP, rec$PrDens, col = "firebrick", lwd = 2)
  
  legend("topright", bty = "n",
         legend = c("Observed date (calibrated)", "Reconstructed (material SPD)"),
         fill   = c(adjustcolor("grey40", alpha.f = 0.5), NA),
         border = c(NA, NA),
         lty    = c(NA, 1),
         lwd    = c(NA, 2),
         col    = c(NA, "firebrick"))
  
  if (!is.null(results_row)) {
    mtext(sprintf("Hellinger = %.3f | Wasserstein = %.1f | Wilcoxon p = %.3f",
                  results_row$hellinger, results_row$wasserstein, results_row$wilcoxon_p),
          side = 3, line = 0.3, cex = 0.8)
  }
}

#function to validate date
validate_date <- function(target_id, make_plot = TRUE, plot_dir = "LOO_plots") {
  
  recon   <- reconstruct_date(target_id, weighting = weighting)
  nullmat <- random_null(target_id, n_used = recon$n_dates_used, n_iter = n_random)
  
  res_row <- data.frame(
    target_id               = target_id,
    n_dates_used            = recon$n_dates_used,
    hellinger               = recon$metrics["hellinger"],
    wasserstein             = recon$metrics["wasserstein"],
    wilcoxon_W              = recon$metrics["wilcoxon_W"],
    wilcoxon_p              = recon$metrics["wilcoxon_p"],
    random_mean_hellinger   = mean(nullmat[, "hellinger"], na.rm = TRUE),
    random_mean_wasserstein = mean(nullmat[, "wasserstein"], na.rm = TRUE),
    random_mean_wilcoxon_W  = mean(nullmat[, "wilcoxon_W"], na.rm = TRUE),
    random_mean_wilcoxon_p  = mean(nullmat[, "wilcoxon_p"], na.rm = TRUE),
    empirical_p_hellinger   = mean(nullmat[, "hellinger"]   <= recon$metrics["hellinger"],   na.rm = TRUE),
    empirical_p_wasserstein = mean(nullmat[, "wasserstein"] <= recon$metrics["wasserstein"], na.rm = TRUE),
    row.names = NULL
  )
  
  if (make_plot) {
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    png(file.path(plot_dir, paste0("LOO_date_", target_id, "_plot.png")),
        width = 1400, height = 900, res = 150)
    plot_reconstruction(recon, res_row)
    dev.off()
  }
  
  list(summary = res_row, reconstruction = recon, null_distribution = nullmat)
}

##### Loop over all dates ======================================================

#Set ID to loop over
loo_ids <- c(1:131)       

loo_runs <- setNames(vector("list", length(loo_ids)), loo_ids)

for (id in loo_ids) {
  message("Validating date ID: ", id)
  loo_runs[[as.character(id)]] <- tryCatch(
    validate_date(id, make_plot = TRUE, plot_dir = "LOO_plots"),
    error = function(e) {
      warning("date ", id, " failed: ", conditionMessage(e))
      NULL
    }
  )
}

# combine all summary rows into one table
all_results <- do.call(rbind, lapply(loo_runs, function(x) if (!is.null(x)) x$summary))
#print(all_results)

write.csv(all_results, "LOO_dates_validation_all_summary.csv", row.names = FALSE)
# saveRDS(loo_runs, "LOO_dates_validation_all_full.rds")


# Mean of each metric across all left-out dates, for observed reconstructions and for the random null model, so you can compare the two at a glance.
metric_cols <- c("hellinger", "wasserstein", "wilcoxon_W", "wilcoxon_p")

overall_summary <- data.frame(
  metric = metric_cols,
  mean_observed = sapply(metric_cols, function(m) mean(all_results[[m]], na.rm = TRUE)),
  mean_random   = sapply(metric_cols, function(m) mean(all_results[[paste0("random_mean_", m)]], na.rm = TRUE)),
  row.names = NULL
)

# mean empirical p-values (how often random reconstructions beat the
# material-informed one, averaged across all left-out dates)
overall_summary_p <- data.frame(
  metric = c("hellinger", "wasserstein"),
  mean_empirical_p = c(
    mean(all_results$empirical_p_hellinger, na.rm = TRUE),
    mean(all_results$empirical_p_wasserstein, na.rm = TRUE)
  )
)

# proportion of dates where the Wilcoxon test rejected H0 (p < 0.05),
# i.e. reconstructed and observed distributions were significantly different
prop_wilcoxon_significant <- mean(all_results$wilcoxon_p < 0.05, na.rm = TRUE)



# percentage of dates where the material-informed reconstruction actually
# beat its own random-null baseline, per metric.
metric_direction <- list(
  hellinger   = "lower",
  wasserstein = "lower",
  wilcoxon_p  = "higher"
)

pct_better_than_random <- sapply(names(metric_direction), function(m) {
  obs <- all_results[[m]]
  rnd <- all_results[[paste0("random_mean_", m)]]
  if (metric_direction[[m]] == "lower") {
    mean(obs < rnd, na.rm = TRUE) * 100
  } else {
    mean(obs > rnd, na.rm = TRUE) * 100
  }
})

pct_better_df <- data.frame(
  metric = names(pct_better_than_random),
  pct_reconstruction_better_than_random = round(unname(pct_better_than_random), 1)
)

write.csv(overall_summary, "LOO_dates_validation_overall_metric_means.csv", row.names = FALSE)
write.csv(pct_better_df, "LOO_dates_validation_pct_better_than_random.csv", row.names = FALSE)

#Clean objects created for this section
rm(material_cols, time_range, n_random, n_resample, weighting,
   get_idx, get_grid, align_grids, hellinger_distance, wasserstein_distance,
   wilcoxon_compare, compare_distributions, build_spd, combine_spds,
   reconstruct_date, random_null, plot_reconstruction, validate_date,
   loo_ids, loo_runs, id, all_results, metric_cols, overall_summary,
   overall_summary_p, prop_wilcoxon_significant, metric_direction,
   pct_better_than_random, pct_better_df)

gc()

#### 1.0.3.2 Site reconstruction ===============================================

material_cols <- c("VAR","OBS","SIL BEDU", "AXE EX", "BRACE", "MONT", "BQ", "CH")

time_range  <- c(6400, 5000)   # calBP range for spd(), adjust to your data span
n_random    <- 1000         # number of random SPD comparisons
n_resample  <- 5000         # sample size drawn from each density for Wilcoxon
weighting   <- "equal"      # "equal" | "count" -> how material SPDs are combined

set.seed(1)                 # for reproducibility of the random null model

##### Helper functions =========================================================

# all row indices in c14_raw / C14_raw_calibration belonging to a given site
get_site_idx <- function(site_name) {
  which(c14_raw$SITE == site_name)
}

# a filesystem-safe version of a site name, for filenames
safe_name <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)

# align two grids on the union of calBP, filling missing years with 0, and re-normalise both to sum to 1
align_grids <- function(g1, g2) {
  bp <- sort(union(g1$calBP, g2$calBP), decreasing = TRUE)
  p1 <- g1$PrDens[match(bp, g1$calBP)]; p1[is.na(p1)] <- 0
  p2 <- g2$PrDens[match(bp, g2$calBP)]; p2[is.na(p2)] <- 0
  list(calBP = bp, p1 = p1 / sum(p1), p2 = p2 / sum(p2))
}

hellinger_distance <- function(p, q) sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)

wasserstein_distance <- function(bp, p, q) {
  ord <- order(bp)
  bp <- bp[ord]; p <- p[ord]; q <- q[ord]
  cdf_diff <- abs(cumsum(p) - cumsum(q))
  w <- abs(diff(bp))
  sum(cdf_diff[-length(bp)] * w)
}

wilcoxon_compare <- function(bp, p, q, n = n_resample) {
  s1 <- sample(bp, n, replace = TRUE, prob = p)
  s2 <- sample(bp, n, replace = TRUE, prob = q)
  wt <- suppressWarnings(wilcox.test(s1, s2))
  c(W = unname(wt$statistic), p_value = wt$p.value)
}

compare_distributions <- function(g1, g2) {
  al <- align_grids(g1, g2)
  hd <- hellinger_distance(al$p1, al$p2)
  wd <- wasserstein_distance(al$calBP, al$p1, al$p2)
  wc <- wilcoxon_compare(al$calBP, al$p1, al$p2)
  c(hellinger   = hd,
    wasserstein = wd,
    wilcoxon_W  = as.numeric(wc["W"]),
    wilcoxon_p  = as.numeric(wc["p_value"]))
}

# build a normalised SPD (calBP/PrDens data.frame) from a set of row indices
build_spd <- function(idx, timeRange = time_range) {
  if (length(idx) == 0) return(NULL)
  sub <- C14_raw_calibration[idx]
  s <- spd(sub, timeRange = timeRange, spdnormalised = TRUE, verbose = FALSE)
  s$grid
}

# combine several material-specific SPDs into one reconstructed distribution
combine_spds <- function(spd_list, weights) {
  bp_all <- sort(unique(unlist(lapply(spd_list, `[[`, "calBP"))), decreasing = TRUE)
  mat <- sapply(seq_along(spd_list), function(i) {
    g <- spd_list[[i]]
    v <- g$PrDens[match(bp_all, g$calBP)]
    v[is.na(v)] <- 0
    v * weights[i]
  })
  combined <- rowSums(mat)
  combined <- combined / sum(combined)
  data.frame(calBP = bp_all, PrDens = combined)
}

# function to reconstruct left-out site
reconstruct_site <- function(target_site, weighting = "equal") {
  
  target_idx <- get_site_idx(target_site)
  if (length(target_idx) == 0) stop("target_site not found in c14_raw$SITE")
  
  # the site's own observed distribution, combining ALL of its dates
  original <- build_spd(target_idx)
  if (is.null(original)) stop("could not build an observed SPD for this site")
  
  # materials associated with the site, aggregated across all of its dates
  site_rows <- c14_raw[target_idx, material_cols, drop = FALSE]
  material_totals <- colSums(data.matrix(site_rows), na.rm = TRUE)  
  present <- material_cols[material_totals > 0]
  if (length(present) == 0) stop("target site has no material associations")
  
  # pool of candidate dates = every date NOT belonging to the target site
  pool <- c14_raw[-target_idx, ]
  
  spd_list <- list()
  weights  <- c()
  used_idx <- c()
  
  for (m in present) {
    ids_m <- pool$ID[pool[[m]] > 0]
    idx_m <- match(ids_m, c14_raw$ID)
    idx_m <- idx_m[!is.na(idx_m)]
    if (length(idx_m) == 0) next
    used_idx <- union(used_idx, idx_m)
    spd_m <- build_spd(idx_m)
    if (is.null(spd_m)) next
    spd_list[[m]] <- spd_m
    weights[m] <- if (weighting == "equal") material_totals[[m]] else 1
  }
  
  if (length(spd_list) == 0) stop("no material-linked dates found for this site")
  
  weights <- weights / sum(weights)
  reconstructed <- combine_spds(spd_list, weights)
  
  metrics <- compare_distributions(reconstructed, original)
  
  list(
    target_site   = target_site,
    target_idx    = target_idx,
    n_own_dates   = length(target_idx),
    materials     = present,
    n_dates_used  = length(used_idx),
    used_idx      = used_idx,
    reconstructed = reconstructed,
    original      = original,
    metrics       = metrics
  )
}

# function for random null model
random_null <- function(target_site, n_used, n_iter = n_random) {
  
  target_idx <- get_site_idx(target_site)
  original   <- build_spd(target_idx)
  pool_idx   <- setdiff(seq_len(nrow(c14_raw)), target_idx)
  
  res <- matrix(NA_real_, nrow = n_iter, ncol = 4,
                dimnames = list(NULL, c("hellinger", "wasserstein",
                                        "wilcoxon_W", "wilcoxon_p")))
  
  for (i in seq_len(n_iter)) {
    samp_idx <- sample(pool_idx, size = min(n_used, length(pool_idx)), replace = FALSE)
    rand_spd <- build_spd(samp_idx)
    if (is.null(rand_spd)) next
    res[i, ] <- compare_distributions(rand_spd, original)
  }
  
  res
}

# function to plot observed site distribution vs reconstructed
plot_reconstruction <- function(recon, results_row = NULL) {
  
  orig <- recon$original
  rec  <- recon$reconstructed
  
  xlim_range <- range(c(orig$calBP, rec$calBP))
  
  op <- par(mar = c(4.5, 4.5, 3, 1))
  on.exit(par(op))
  
  plot(orig$calBP, orig$PrDens, type = "n",
       xlim = rev(xlim_range),                 # calBP: older (larger) on the left
       ylim = c(0, max(orig$PrDens, rec$PrDens) * 1.1),
       xlab = "cal BP", ylab = "Probability density",
       main = paste0(recon$target_site,
                     " - observed vs. reconstructed date (n=",
                     recon$n_own_dates, " own dates)"))
  
  # observed (true) site distribution, filled
  polygon(c(orig$calBP, rev(orig$calBP)),
          c(orig$PrDens, rep(0, length(orig$PrDens))),
          col = adjustcolor("grey40", alpha.f = 0.5), border = NA)
  
  # reconstructed distribution, outline only
  lines(rec$calBP, rec$PrDens, col = "firebrick", lwd = 2)
  
  legend("topright", bty = "n",
         legend = c("Observed site SPD", "Reconstructed (material SPD)"),
         fill   = c(adjustcolor("grey40", alpha.f = 0.5), NA),
         border = c(NA, NA),
         lty    = c(NA, 1),
         lwd    = c(NA, 2),
         col    = c(NA, "firebrick"))
  
  if (!is.null(results_row)) {
    mtext(sprintf("Hellinger = %.3f | Wasserstein = %.1f | Wilcoxon p = %.3f",
                  results_row$hellinger, results_row$wasserstein, results_row$wilcoxon_p),
          side = 3, line = 0.3, cex = 0.8)
  }
}

# function to validate sites
validate_site <- function(target_site, make_plot = TRUE, plot_dir = "LOO_site_plots") {
  
  recon   <- reconstruct_site(target_site, weighting = weighting)
  nullmat <- random_null(target_site, n_used = recon$n_dates_used, n_iter = n_random)
  
  res_row <- data.frame(
    target_site              = target_site,
    n_own_dates               = recon$n_own_dates,
    n_dates_used              = recon$n_dates_used,
    hellinger                 = recon$metrics["hellinger"],
    wasserstein                = recon$metrics["wasserstein"],
    wilcoxon_W                = recon$metrics["wilcoxon_W"],
    wilcoxon_p                = recon$metrics["wilcoxon_p"],
    random_mean_hellinger      = mean(nullmat[, "hellinger"], na.rm = TRUE),
    random_mean_wasserstein    = mean(nullmat[, "wasserstein"], na.rm = TRUE),
    random_mean_wilcoxon_W     = mean(nullmat[, "wilcoxon_W"], na.rm = TRUE),
    random_mean_wilcoxon_p     = mean(nullmat[, "wilcoxon_p"], na.rm = TRUE),
    empirical_p_hellinger      = mean(nullmat[, "hellinger"]   <= recon$metrics["hellinger"],   na.rm = TRUE),
    empirical_p_wasserstein    = mean(nullmat[, "wasserstein"] <= recon$metrics["wasserstein"], na.rm = TRUE),
    row.names = NULL
  )
  
  if (make_plot) {
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    png(file.path(plot_dir, paste0("LOO_site_", safe_name(target_site), "_plot.png")),
        width = 1400, height = 900, res = 150)
    plot_reconstruction(recon, res_row)
    dev.off()
  }
  
  list(summary = res_row, reconstruction = recon, null_distribution = nullmat)
}


##### Loop over all dates ======================================================

loo_sites <- unique(c14_raw$SITE)   # <-- edit: e.g. unique(c14_raw$SITE) for every site

loo_runs <- setNames(vector("list", length(loo_sites)), loo_sites)
loo_failures <- list()   # <-- track failed site names and the reason

for (s in loo_sites) {
  message("Validating site: ", s)
  loo_runs[[s]] <- tryCatch(
    validate_site(s, make_plot = TRUE, plot_dir = "LOO_site_plots"),
    error = function(e) {
      msg <- conditionMessage(e)
      warning("Site ", s, " failed: ", msg)
      loo_failures[[s]] <<- msg
      NULL
    }
  )
}

# quick look at which sites were skipped, and why
if (length(loo_failures) > 0) {
  failures_df <- data.frame(
    target_site = names(loo_failures),
    reason      = unlist(loo_failures),
    row.names = NULL
  )
  cat("\n=== Sites skipped (no plot / no result) ===\n")
  print(failures_df)
  write.csv(failures_df, "LOO_sites_validation_failures.csv", row.names = FALSE)
} else {
  cat("\nAll sites in loo_sites were validated successfully.\n")
}

# combine all summary rows into one table
all_results <- do.call(rbind, lapply(loo_runs, function(x) if (!is.null(x)) x$summary))
# print(all_results)

write.csv(all_results, "LOO_sites_validation_all_summary.csv", row.names = FALSE)
# saveRDS(loo_runs, "LOO_sites_validation_all_full.rds")

# Mean of each metric across all left-out sites, for observed reconstructions
# and for the random null model, so you can compare the two at a glance.
metric_cols <- c("hellinger", "wasserstein", "wilcoxon_W", "wilcoxon_p")

overall_summary <- data.frame(
  metric = metric_cols,
  mean_observed = sapply(metric_cols, function(m) mean(all_results[[m]], na.rm = TRUE)),
  mean_random   = sapply(metric_cols, function(m) mean(all_results[[paste0("random_mean_", m)]], na.rm = TRUE)),
  row.names = NULL
)

# mean empirical p-values (how often random reconstructions beat the
# material-informed one, averaged across all left-out sites)
overall_summary_p <- data.frame(
  metric = c("hellinger", "wasserstein"),
  mean_empirical_p = c(
    mean(all_results$empirical_p_hellinger, na.rm = TRUE),
    mean(all_results$empirical_p_wasserstein, na.rm = TRUE)
  )
)

# proportion of sites where the Wilcoxon test rejected H0 (p < 0.05),
# i.e. reconstructed and observed distributions were significantly different
prop_wilcoxon_significant <- mean(all_results$wilcoxon_p < 0.05, na.rm = TRUE)


# percentage of sites where the material-informed reconstruction actually
# beat its own random-null baseline, per metric.
metric_direction <- list(
  hellinger   = "lower",
  wasserstein = "lower",
  wilcoxon_p  = "higher"
)

pct_better_than_random <- sapply(names(metric_direction), function(m) {
  obs <- all_results[[m]]
  rnd <- all_results[[paste0("random_mean_", m)]]
  if (metric_direction[[m]] == "lower") {
    mean(obs < rnd, na.rm = TRUE) * 100
  } else {
    mean(obs > rnd, na.rm = TRUE) * 100
  }
})

pct_better_df <- data.frame(
  metric = names(pct_better_than_random),
  pct_reconstruction_better_than_random = round(unname(pct_better_than_random), 1)
)

write.csv(overall_summary, "LOO_sites_validation_overall_metric_means.csv", row.names = FALSE)
write.csv(pct_better_df, "LOO_sites_validation_pct_better_than_random.csv", row.names = FALSE)

#Erase all objects created during this section
rm(list = intersect(ls(), c("material_cols","time_range","n_random","n_resample",
                            "weighting","get_site_idx","safe_name","align_grids","hellinger_distance",
                            "wasserstein_distance","wilcoxon_compare","compare_distributions","build_spd",
                            "combine_spds","reconstruct_site","random_null","plot_reconstruction",
                            "validate_site","loo_sites","loo_runs","loo_failures","s","all_results",
                            "metric_cols","overall_summary","overall_summary_p","prop_wilcoxon_significant",
                            "metric_direction","pct_better_than_random","pct_better_df","failures_df")))

gc()

## 1.1 SPD total ================================================================

#Create output folder
dir.create("Results/plots", showWarnings = FALSE)

spd_burials <- spd(
  C14_raw_calibration,
  timeRange = c(6400, 5000),   # range in BP
  spdnormalised = TRUE
)

tiff(
  "Results/1_SPD/spd_burials.tiff",
  width = 3600,
  height = 2400,
  res = 600
)

plot(
  spd_burials,
  runm = 200,                  # 200-year rolling mean
  type = "simple",
  col = "darkorange",
  lwd = 1.5,
  lty = 2
)

dev.off()

## 1.2 SPD by variable ==========================================================

# Variables to analyse
vars <- c("OBS","SIL BEDU","AXE EX","BRACE","VAR","MOL","MONT","BQ","CH")

# Function to compute normalised SPD by variable
spd_by_material <- function(data, variable, time_range = c(6400, 5000)) {
  
  # Check if the column contains at least one presence
  if (sum(data[[variable]] == 1, na.rm = TRUE) == 0) {
    message(paste("Skipping", variable, "- all values are 0"))
    return(NULL)
  }
  
  subset_data <- data[data[[variable]] == 1, ]
  
  cal <- rcarbon::calibrate(
    x = subset_data$DATE,
    errors = subset_data$SD,
    ids = subset_data$ID,
    calCurves = "intcal20"
  )
  
  spd_res <- spd(
    cal,
    timeRange = time_range,
    spdnormalised = TRUE
  )
  
  return(spd_res)
}

# Generate SPD for all variables
spd_list <- lapply(vars, function(v) spd_by_material(c14_raw, v))
names(spd_list) <- vars

# Remove variables without radiocarbon dates
spd_list <- spd_list[!sapply(spd_list, is.null)]

vars_valid <- names(spd_list)

### 1.2.1 Export combined SPD plot =================================================

tiff("Results/1_SPD/spd_styles.tiff", width = 3600, height = 2400, res = 600)

cols <- rainbow(length(spd_list))

# First SPD (just set the first plot)
plot(spd_list[[1]], runm = 200, border = cols[1], lwd = 2)

# Add remaining SPDs
for (i in 2:length(spd_list)) {
  plot(spd_list[[i]], runm = 200, add = TRUE, border = cols[i], lwd = 2)
}

# Add labels and title separately
title(main = "Summed Probability Distributions by style")

legend(
  "topright",
  legend = vars_valid,
  col = cols,
  lwd = 2,
  bty = "n"
)

dev.off()

### 1.2.2 Export individual SPD plots =============================================

pdf(
  "Results/1_SPD/spd_individual_styles.pdf",
  width = 10,
  height = 8
)

n <- length(spd_list)

# Automatic grid layout
rows <- ceiling(sqrt(n))
cols <- ceiling(n / rows)

par(mfrow = c(rows, cols))

for (i in 1:n) {
  
  plot(
    spd_list[[i]],
    runm = 200,
    border = "darkblue",
    lwd = 2
  )
  
  title(main = vars_valid[i])
}

dev.off()


png(
  "Results/1_SPD/spd_individual_styles.png",
  width = 16,
  height = 12,
  units = "in",
  res = 300
)

n <- length(spd_list)

# Automatic grid layout
rows <- ceiling(sqrt(n))
cols <- ceiling(n / rows)

par(mfrow = c(rows, cols))

for (i in 1:n) {
  
  plot(
    spd_list[[i]],
    runm = 200,
    border = "darkblue",
    lwd = 2
  )
  
  title(main = vars_valid[i])
}

dev.off()

### 1.2.3 Extract SPD grids ========================================================

spd_grids <- lapply(spd_list, function(spd_obj) {
  
  data.frame(
    calBP = spd_obj$grid$calBP,
    PrDens = spd_obj$grid$PrDens
  )
  
})

names(spd_grids) <- names(spd_list)

### 1.2.4 Extract metrics ==========================================================

# Initialize list to store results
summary_list <- list()

for (i in 1:length(spd_grids)) {
  
  spd_df <- spd_grids[[i]]
  variable_name <- names(spd_grids)[i]
  
  # Number of radiocarbon dates used
  n_dates <- sum(c14_raw[[variable_name]] == 1, na.rm = TRUE)
  
  if (n_dates < 3) {  # too few dates
    median_bp <- iqr_bp <- sd_bp <- NA
  } else {
    # Weighted median
    w_cumsum <- cumsum(spd_df$PrDens) / sum(spd_df$PrDens)
    median_bp <- spd_df$calBP[which.min(abs(w_cumsum - 0.5))]
    
    # Weighted mean & SD
    mean_bp <- sum(spd_df$calBP * spd_df$PrDens) / sum(spd_df$PrDens)
    sd_bp <- sqrt(sum(spd_df$PrDens * (spd_df$calBP - mean_bp)^2) / sum(spd_df$PrDens))
    
    # Weighted IQR
    lower <- spd_df$calBP[which.min(abs(w_cumsum - 0.25))]
    upper <- spd_df$calBP[which.min(abs(w_cumsum - 0.75))]
    iqr_bp <- upper - lower
  }
  
  # Store in list
  summary_list[[i]] <- data.frame(
    Material = variable_name,
    n = n_dates,
    `Median cal BP` = ifelse(is.na(median_bp), "—", round(median_bp)),
    `IQR (years)` = ifelse(is.na(iqr_bp), "—", round(iqr_bp)),
    `SD (years)` = ifelse(is.na(sd_bp), "—", round(sd_bp)),
    Peaks = "",                 # blank
    `Temporal structure` = "",  # blank
    stringsAsFactors = FALSE
  )
}

# Combine all results
summary_table <- do.call(rbind, summary_list)

# Export CSV
write.csv(summary_table, "Results/1_SPD/SPD_summary_table.csv", row.names = FALSE)

# 2 C14 STKDE ANALYSIS #########################################################

## 2.1 Create curves for Spatio-temporal KDE  ===================================

# Spatial analysis window
ref_raster <- rast("Data/Rasters/MDE_cat_and_100.tif")
e <- ext(ref_raster)

win <- owin(
  xrange = c(e$xmin, e$xmax),
  yrange = c(e$ymin, e$ymax)
)

# Create STKDE main output directory

dir.create("Results/2_C14_STKDE", showWarnings = FALSE, recursive = TRUE)

# Function to compute STKDE by variable

stkde_total <- function(data, win){
  
  subset_data <- data
  
  if(nrow(subset_data) == 0){
    message("No dates available")
    return(NULL)
  }
  
  # Compute weights
  weights <- 1 / subset_data$n_dates
  
  # Output directory
  out_dir <- "Results/2_C14_STKDE/stkde_total"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Calibration
  cal <- rcarbon::calibrate(
    x = subset_data$DATE,
    errors = subset_data$SD,
    ids = subset_data$ID,
    calCurves = "intcal20",
    normalised = FALSE,
    verbose = FALSE
  )
  
  # Temporal binning
  bins <- binPrep(
    sites = subset_data$SITE,
    ages = subset_data$DATE,
    h = 50
  )
  
  # Coordinates
  coords <- as.matrix(subset_data[,c("X","Y")])
  
  # Spatio-Temporal KDE
  stkde_res <- stkde(
    x = cal,                         
    coords = coords,                 
    win = win,                       
    sbw = 10000,                     # bandwidth espacial del kernel (20 km de suavizado)
    tbw = 50,                        # bandwidth temporal del kernel (±50 años)
    cellres = 1000,                  # resolución del raster de salida (celdas de 2 km)
    focalyears = seq(6400, 5000, -100), # años para los que se calcula el STKDE (cada 100 años)
    bins = bins,                     # binning temporal para reducir fechas muy próximas del mismo sitio
    backsight = 0,                 # ventana temporal hacia atrás considerada en el cálculo
    weights = weights,               # peso de cada fecha (corrige sobre-representación de sitios)
    outdir = out_dir,                
    amount = 0.5,                    # jitter espacial para evitar puntos idénticos
    verbose = FALSE                  # suprime mensajes del proceso
  )
  
  return(stkde_res)
}

# Run global STKDE

stkde_res <- stkde_total(c14_raw, win)

## 2.2. Export maps =============================================================

#Create output folder
dir.create("Results/2_C14_STKDE/maps", showWarnings = FALSE, recursive = TRUE)

for(y in seq_along(stkde_res[["impaths"]])){
    
    load(stkde_res[["impaths"]][[y]])
    
    r <- rast(focalyear[["focal"]])             
    
    min_val <- min(values(r), na.rm = TRUE)
    max_val <- max(values(r), na.rm = TRUE)
    
    r_norm <- (r - min_val) / (max_val - min_val)
    
    crs(r_norm) <- crs(ref_raster)
    
    # Reproject to reference CRS
    if (crs(r_norm) != crs(ref_raster)) {
      r_norm <- project(r_norm, ref_raster, method = "bilinear")  # continuous
    }
    
    # Extend to reference extent
    if (!all(ext(r_norm) == ext(ref_raster))) {
      r_norm <- extend(r_norm, ext(ref_raster))
    }
    
    # Match resolution
    if (!all(res(r_norm) == res(ref_raster))) {
      r_norm <- terra::resample(r_norm, ref_raster, method = "bilinear")
    }
    
    r_norm[is.na(r_norm)] <- 0
    
    # Apply mask using the reference raster
    r_norm <- mask(r_norm, ref_raster)
    
    writeRaster(r_norm,  paste0("Results/2_C14_STKDE/maps/STKDE_", focalyear[["year"]], ".tiff"), overwrite=TRUE)
    
  }


# 3 STOCHASTIC STKDE SIMULATION OF UNDATED SITES ###############################
## 3.1 Create SPD for non dated sites ===========================================

# Get unique sites
ids <- unique(no_c14_raw$ID)

# Define time range for all SPDs
time_range <- c(6400, 5000)  # adjust if needed

# Initialize a list to store summed SPDs per site
spd_per_ids <- list()

sum_spd_objects  <- function(spd_objects){
  
  # Extract the numeric probability grids
  grids <- lapply(spd_objects, function(spd) spd$grid$PrDens)
  
  # Make sure they all have the same length
  len <- sapply(grids, length)
  if(length(unique(len)) > 1) stop("SPD grids have different lengths!")
  
  # Sum the numeric grids
  summed_grid <- Reduce(`+`, grids)
  
  # Return a new CalSPD object using the first one as template
  summed_spd <- spd_objects[[1]]
  summed_spd$grid$PrDens <- summed_grid
  
  return(summed_spd)
}

for(id in ids){
  
  # Subset rows for the site
  id_data <- no_c14_raw[no_c14_raw$ID == id, ]
  
  # Find variables present at this site
  vars_present <- vars[sapply(vars, function(v) any(id_data[[v]] == 1, na.rm = TRUE))]
  
  vars_present <- vars_present[vars_present %in% names(spd_list)]
  
  if(length(vars_present) > 0){
    # Sum SPDs of variables present at this site
    spd_objects <- spd_list[vars_present]
    spd_objects <- spd_objects[!sapply(spd_objects, is.null)]
    
    summed_spd <- sum_spd_objects(spd_objects)
    
  } else {
    # No variables present -> Use general SPD
    # Use the first SPD as template if available, otherwise create a basic one
    template_spd <- spd_list[[1]]
    summed_spd <- template_spd
    summed_spd$grid$PrDens <- spd_burials$grid$PrDens
    
  }
  
  spd_per_ids[[as.character(id)]] <- summed_spd
  
}

## 3.2 Create Random spatial kernels for each date ==============================

# Function to convert CalSPD to CalDates
spd_to_calDates <- function(spd_list, ids=NULL){
  
  if(is.null(ids)){
    ids <- names(spd_list)
  }
  
  grids <- lapply(spd_list, function(x){
    
    p <- x$grid$PrDens
    p <- p / sum(p)   # normalise
    
    data.frame(
      calBP = x$grid$calBP,
      PrDens = p
    )
    
  })
  
  calBP <- grids[[1]]$calBP
  
  caldates <- list(
    metadata = data.frame(
      CRA = NA,
      Error = NA,
      ID = ids
    ),
    grids = grids,
    calBP = calBP
  )
  
  class(caldates) <- "CalDates"
  
  return(caldates)
}

# Convert all SPD to CalDate
all_cal <- spd_to_calDates(spd_per_ids, ids = names(spd_per_ids))

# Output directory
base_out <- "Results/3_RNC14_STKDE"
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)

ncores <- detectCores() - 1
clus <- makeCluster(ncores)

# Load required packages on workers
clusterEvalQ(clus, {
  library(dplyr)
  library(rcarbon)
})

# Export objects needed by workers
clusterExport(
  clus,
  varlist = c(
    "no_c14_raw",
    "spd_per_ids",
    "all_cal",
    "win",
    "base_out"
  ),
  envir = environment()
)

# Worker function
worker_stkde <- function(x3) {
  
  out_dir <- file.path(base_out, x3)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  idx <- sample(
    seq_len(nrow(no_c14_raw)),
    size = floor(0.5 * nrow(no_c14_raw)),
    replace = FALSE
  )
  
  no_c14_raw_sub <- no_c14_raw[idx, ]
  
  # Correction for site overrepresentation
  site_counts <- no_c14_raw_sub %>%
    dplyr::group_by(SITE) %>%
    dplyr::summarise(n_dates = n(), .groups = "drop")
  
  no_c14_raw_sub <- no_c14_raw_sub %>%
    dplyr::left_join(site_counts, by = "SITE")
  
  weights <- 1 / no_c14_raw_sub$n_dates
  
  # Retrieve precomputed SPD
  spd_subset <- spd_per_ids[as.character(no_c14_raw_sub$ID)]
  
  # Convert SPD list to CalDates
  cal <- list(
    metadata = all_cal$metadata[idx, ],
    grids = all_cal$grids[idx],
    calBP = all_cal$calBP
  )
  
  class(cal) <- "CalDates"
  
  coords <- as.matrix(no_c14_raw_sub[, c("X", "Y")])
  
  stkde_res <- stkde(
    x = cal,
    coords = coords,
    win = win,
    sbw = 10000,
    tbw = 50,
    cellres = 1000,
    focalyears = seq(6400, 5000, -100),
    backsight = 0,
    weights = weights,
    outdir = out_dir,
    amount = 0.5,
    verbose = FALSE
  )
  
  return(stkde_res)
}

# Run in parallel
random_stkde <- parLapply(clus, 1:100, worker_stkde)

stopCluster(clus)

## 3.3 Export maps ==============================================================

#Create output folder
out_dir <- "Results/4_RASTERS_RNC14_STKDE"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ncores <- detectCores() - 1
clus <- makeCluster(ncores)

clusterEvalQ(clus, {
  library(terra)
})

clusterExport(
  clus,
  varlist = c(
    "random_stkde"
  ),
  envir = environment()
)

worker_rasters <- function(y){
  
  for (date in seq(6400, 5000, -100)) {
    
    out_dir <- paste0("Results/4_RASTERS_RNC14_STKDE/", date)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    ref_raster <- rast("Data/Rasters/MDE_cat_and_100.tif")
    
    load(random_stkde[[y]][["impaths"]][[as.character(date)]])
    
    r <- rast(focalyear[["focal"]])
    
    min_val <- min(values(r), na.rm = TRUE)
    max_val <- max(values(r), na.rm = TRUE)
    
    r_norm <- (r - min_val) / (max_val - min_val)
    
    crs(r_norm) <- crs(ref_raster)
    
    if (crs(r_norm) != crs(ref_raster)) {
      r_norm <- project(r_norm, ref_raster, method = "bilinear")
    }
    
    if (!all(ext(r_norm) == ext(ref_raster))) {
      r_norm <- extend(r_norm, ext(ref_raster))
    }
    
    if (!all(res(r_norm) == res(ref_raster))) {
      r_norm <- terra::resample(r_norm, ref_raster, method = "bilinear")
    }
    
    r_norm[is.na(r_norm)] <- 0
    
    r_norm <- mask(r_norm, ref_raster)
    
    writeRaster(
      r_norm,
      paste0("Results/4_RASTERS_RNC14_STKDE/", focalyear[["year"]], "/", y, ".tiff"),
      overwrite = TRUE
    )
  }
  
  return(NULL)
}

parLapply(clus, seq_along(random_stkde), worker_rasters)

stopCluster(clus)

## 3.4 Mean of all resulting RNC14 STKDE  =======================================

## Mean the corrected KDEs
#List all box folders
all_folders <- list.dirs("Results/4_RASTERS_RNC14_STKDE", recursive = FALSE)

out_dir <- "Results/4_RASTERS_RNC14_STKDE/Mean_RNC14_STKDE"
dir.create(out_dir, showWarnings = FALSE)

ncores <- detectCores() - 1
clus <- makeCluster(ncores)

clusterEvalQ(clus, {
  library(terra)
})

clusterExport(
  clus,
  varlist = c("all_folders", "out_dir"),
  envir = environment()
)

worker_mean <- function(f){
  
  folder <- all_folders[f]
  
  folder_rs <- list.files(folder, pattern="\\.tif[f]?$", full.names=TRUE)
  
  # stack rasters directly from filenames
  rs_stack <- rast(folder_rs)
  
  # compute mean
  mean_raster <- app(rs_stack, mean, na.rm=TRUE)
  
  # normalization
  min_val <- global(mean_raster, "min", na.rm=TRUE)[1,1]
  max_val <- global(mean_raster, "max", na.rm=TRUE)[1,1]
  
  normalized_raster <- (mean_raster - min_val) / (max_val - min_val)
  
  writeRaster(
    normalized_raster,
    file.path(out_dir,
              paste0("Mean_RNC14_STKDE_", basename(folder), ".tiff")),
    overwrite=TRUE
  )
  
  return(NULL)
}

parLapplyLB(clus, seq_along(all_folders), worker_mean)

stopCluster(clus)


# 4 CORRECTION OF C14 STKDE ####################################################

## 4.1 Multiplication of C14 STKDE with RNC14STKDE  =============================

#Output directory
dir.create("Results/5_CORRECTED_STKDE", showWarnings = FALSE)

dates <- seq(6400, 5000, -100)

ncores <- detectCores() - 1
clus <- makeCluster(ncores)

clusterEvalQ(clus, {
  library(terra)
})

clusterExport(
  clus,
  varlist = c("dates"),
  envir = environment()
)

worker_corrected <- function(date){
  
  out_dir <- file.path("Results/5_CORRECTED_STKDE", date)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Temporal KDE
  temp_kde_rast <- rast(
    file.path("Results/2_C14_STKDE/maps",
              paste0("STKDE_", date, ".tiff"))
  )
  
  # Randomized KDE rasters for this date
  folder_rs <- list.files(
    file.path("Results/4_RASTERS_RNC14_STKDE", date),
    pattern="\\.tif[f]?$",
    full.names=TRUE
  )
  
  for (rs in seq_along(folder_rs)) {
    
    rs_kde_rast <- rast(folder_rs[rs])
    
    corrected_kde <- (0.5 * temp_kde_rast) + (0.5 * rs_kde_rast)
    
    writeRaster(
      corrected_kde,
      file.path(out_dir,
                paste0("corrected_kde_", rs, ".tiff")),
      overwrite=TRUE
    )
  }
  
  return(NULL)
}

parLapplyLB(clus, dates, worker_corrected)

stopCluster(clus)

## 4.2 Mean of all resulting STKDE ==============================================

# List folders
all_folders <- list.dirs("Results/5_CORRECTED_STKDE", recursive = FALSE)

# Output directory
out_dir <- "Results/6_MEAN_CORRECTED_STKDE"
dir.create(out_dir, showWarnings = FALSE)

ncores <- detectCores() - 1
clus <- makeCluster(ncores)

clusterEvalQ(clus, {
  library(terra)
})

clusterExport(
  clus,
  varlist = c("all_folders", "out_dir"),
  envir = environment()
)

worker_final_mean <- function(f){
  
  folder <- all_folders[f]
  
  # retrieve rasters
  folder_rs <- list.files(folder, pattern="\\.tif[f]?$", full.names=TRUE)
  
  # stack directly from filenames
  rs_stack <- rast(folder_rs)
  
  # compute mean
  mean_raster <- app(rs_stack, mean, na.rm=TRUE)
  
  # normalize
  min_val <- global(mean_raster, "min", na.rm=TRUE)[1,1]
  max_val <- global(mean_raster, "max", na.rm=TRUE)[1,1]
  
  normalized_raster <- (mean_raster - min_val) / (max_val - min_val)
  
  writeRaster(
    normalized_raster,
    file.path(out_dir,
              paste0("Mean_STKDE_", basename(folder), ".tiff")),
    overwrite=TRUE
  )
  
  return(NULL)
}

parLapplyLB(clus, seq_along(all_folders), worker_final_mean)

stopCluster(clus)

## 4.3 Convert results to GIF ===================================================

## Final model
# Folder containing rasters
folder_mkde <- "Results/6_MEAN_CORRECTED_STKDE"

#Stablish dates
dates <- seq(6400, 5000, -100)

# Read all .tif rasters in folder as single-layer SpatRasters
rasters_list <- lapply(list.files(folder_mkde, pattern = "\\.tif[f]?$", full.names = TRUE), rast)

# Ensure single-layer rasters (in case any accidentally have multiple layers)
rasters_list <- lapply(rasters_list, function(r) if (nlyr(r) > 1) r[[1]] else r)

# Find min and max across all rasters for consistent color scale
min_val <- min(sapply(rasters_list, function(x) min(values(x), na.rm = TRUE)))
max_val <- max(sapply(rasters_list, function(x) max(values(x), na.rm = TRUE)))

# Create a list of images
img_list <- lapply(seq_along(rasters_list), function(i) {
  # Open graphics device in memory
  img <- image_graph(width = 600, height = 600, res = 96)
  
  # Plot raster with consistent color scale
  plot(rasters_list[[i]], 
       main = paste(dates[i],"cal BP"),
       zlim = c(min_val, max_val),       # fixed scale
       col = viridis(100))
  
  # Close device and capture image
  dev.off()
  img
})

# Create GIF
gif <- image_animate(image_join(img_list), fps = 0.5)

# Save GIF
image_write(gif, "Results/plots/Final_model.gif")

## Base model
# Folder containing rasters
folder_mkde_2 <- "Results/2_C14_STKDE/maps"

#Stablish dates
dates <- seq(6400, 5000, -100)

# Read all .tif rasters in folder as single-layer SpatRasters
rasters_list <- lapply(list.files(folder_mkde_2, pattern = "\\.tif[f]?$", full.names = TRUE), rast)

# Ensure single-layer rasters (in case any accidentally have multiple layers)
rasters_list <- lapply(rasters_list, function(r) if (nlyr(r) > 1) r[[1]] else r)

# Find min and max across all rasters for consistent color scale
min_val <- min(sapply(rasters_list, function(x) min(values(x), na.rm = TRUE)))
max_val <- max(sapply(rasters_list, function(x) max(values(x), na.rm = TRUE)))

# Create a list of images
img_list <- lapply(seq_along(rasters_list), function(i) {
  # Open graphics device in memory
  img <- image_graph(width = 600, height = 600, res = 96)
  
  # Plot raster with consistent color scale
  plot(rasters_list[[i]], 
       main = paste(dates[i],"cal BP"),
       zlim = c(min_val, max_val),       # fixed scale
       col = viridis(100))
  
  # Close device and capture image
  dev.off()
  img
})

# Create GIF
gif_2 <- image_animate(image_join(img_list), fps = 0.5)

# Save GIF
image_write(gif_2, "Results/plots/Base_model.gif")

## RNC14_STKDE model
folder_mkde_3 <- "Results/4_RASTERS_RNC14_STKDE/Mean_RNC14_STKDE"

#Stablish dates
dates <- seq(6400, 5000, -100)

# Read all .tif rasters in folder as single-layer SpatRasters
rasters_list <- lapply(list.files(folder_mkde_3, pattern = "\\.tif[f]?$", full.names = TRUE), rast)

# Ensure single-layer rasters (in case any accidentally have multiple layers)
rasters_list <- lapply(rasters_list, function(r) if (nlyr(r) > 1) r[[1]] else r)

# Find min and max across all rasters for consistent color scale
min_val <- min(sapply(rasters_list, function(x) min(values(x), na.rm = TRUE)))
max_val <- max(sapply(rasters_list, function(x) max(values(x), na.rm = TRUE)))

# Create a list of images
img_list <- lapply(seq_along(rasters_list), function(i) {
  # Open graphics device in memory
  img <- image_graph(width = 600, height = 600, res = 96)
  
  # Plot raster with consistent color scale
  plot(rasters_list[[i]], 
       main = paste(dates[i],"cal BP"),
       zlim = c(min_val, max_val),       # fixed scale
       col = viridis(100))
  
  # Close device and capture image
  dev.off()
  img
})

# Create GIF
gif_3 <- image_animate(image_join(img_list), fps = 0.5)

# Save GIF
image_write(gif_3, "Results/plots/RNC14_model.gif")

## Unite all GIFS
gif_base  <- image_read("Results/plots/base_model.gif")
gif_random <- image_read("Results/plots/RNC14_model.gif")
gif_final <- image_read("Results/plots/Final_model.gif")

# Ensure same number of frames
n <- min(length(gif_base), length(gif_random), length(gif_final))

frames <- vector("list", n)

for(i in 1:n){
  
  frame_base  <- image_join(gif_base[i])
  frame_random  <- image_join(gif_random[i])
  frame_final <- image_join(gif_final[i])
  
  # Combine horizontally
  combined <- image_append(c(frame_base,frame_random, frame_final))
  
  # Get dimensions
  info <- image_info(combined)
  width <- info$width
  
  # Create title bar
  title_bar <- image_blank(width = width, height = 60, color = "white")
  
  title_bar <- image_annotate(
    title_bar,
    text = "Base model (C14 STKDE) vs. Random model (RNC14 STKDE) vs. Final model (Corrected STKDE)",
    size = 28,
    gravity = "center"
  )
  
  # Stack title + image
  frames[[i]] <- image_append(c(title_bar, combined), stack = TRUE)
}

combined_all <- image_join(frames)

gif_combined <- image_animate(combined_all, fps = 0.5)

image_write(gif_combined, "Results/plots/Model_comparison.gif")

# 5. SEQUENTIAL MODEL ANALYSIS #################################################

## 5.1 Raster value extraction ==================================================

#Load final models
folder_models <- "Results/6_MEAN_CORRECTED_STKDE"
final_model_list <- lapply(list.files(folder_models, pattern = "\\.tif[f]?$", full.names = TRUE), rast)
final_model_list <- rast(final_model_list)
names(final_model_list) <- seq(5000, 6400, 100) #Assign names

#Obtain values
vals <- values(final_model_list, na.rm=TRUE)
n_layers <- ncol(vals)
layer_names <- names(final_model_list)

## 5.2 Pearson correlation test =================================================

vals_cor <- cor(vals, method = "pearson")

# Save results
write.csv(vals_cor, "Results/plots/Sequential_Pearson_Correlation.csv", row.names = TRUE)

## 5.3 RMSE computation =========================================================

rmse_results <- matrix(NA, n_layers, n_layers, dimnames = list(layer_names, layer_names))
pairs <- combn(n_layers, 2)
apply(pairs, 2, function(idx){
  i <- idx[1]; j <- idx[2]
  rmse_val <- rmse(vals[,i], vals[,j])
  rmse_results[i,j] <<- rmse_val
  rmse_results[j,i] <<- rmse_val
})
diag(rmse_results) <- 0  # RMSE with itself = 0

# Save results
write.csv(rmse_results, "Results/plots/Sequential_RMSE.csv", row.names = TRUE)

## 5.4 Compute Hellinger distance and Bhattacharyya coefficient =================

# Normalize values for probability distribution 
vals_prob <- apply(vals, 2, function(x) x / sum(x, na.rm = TRUE))

# Compute Hellinger distance and Bhattacharyya coefficient
hellinger_results <- matrix(NA, n_layers, n_layers, dimnames = list(layer_names, layer_names))
bhattacharyya_results <- matrix(NA, n_layers, n_layers, dimnames = list(layer_names, layer_names))

pairs <- combn(n_layers, 2)
apply(pairs, 2, function(idx){
  i <- idx[1]; j <- idx[2]
  
  p <- vals_prob[, i]
  q <- vals_prob[, j]
  
  # Bhattacharyya coefficient
  bc <- sum(sqrt(p * q), na.rm = TRUE)
  bhattacharyya_results[i,j] <<- bc
  bhattacharyya_results[j,i] <<- bc
  
  # Hellinger distance
  h <- sqrt(1 - bc)
  hellinger_results[i,j] <<- h
  hellinger_results[j,i] <<- h
})

diag(bhattacharyya_results) <- 1   # BC with itself
diag(hellinger_results) <- 0       # Hellinger distance with itself

# Save results
write.csv(bhattacharyya_results, "Results/plots/Sequential_Bhattacharyya.csv", row.names = TRUE)
write.csv(hellinger_results, "Results/plots/Sequential_Hellinger.csv", row.names = TRUE)


## 5.5 Sequential comparison ====================================================

# Sequential values for adjacent raster layers
seq_corr <- diag(vals_cor[-1, -ncol(vals_cor)])
seq_rmse <- diag(rmse_results[-1, -ncol(rmse_results)])
seq_hellinger <- diag(hellinger_results[-1, -ncol(hellinger_results)])
seq_bhat <- diag(bhattacharyya_results[-1, -ncol(bhattacharyya_results)])

# X axis values
x_vals <- as.numeric(names(final_model_list)[-length(names(final_model_list))])

# Create dataframes for individual plots
corr_df <- data.frame(x = x_vals, y = seq_corr)
rmse_df <- data.frame(x = x_vals, y = seq_rmse)
hellinger_df <- data.frame(x = x_vals, y = seq_hellinger)
bhat_df <- data.frame(x = x_vals, y = seq_bhat)

## Individual plots
# Pearson correlation
corr_plot <- ggplot(corr_df, aes(x = x, y = y)) +
  geom_line(color = "steelblue", size = 1) + geom_point(color = "steelblue", linewidth = 2) +
  scale_x_reverse() +
  labs(x = "Sequential Model Pair (Start Value)", y = "Pearson Correlation",
       title = "Sequential Pearson Correlation") +
  theme_minimal()

ggsave("Results/plots/Sequential_Pearson_Correlation.tiff", corr_plot,
       width = 8, height = 4, dpi = 300)

# RMSE
rmse_plot <- ggplot(rmse_df, aes(x = x, y = y)) +
  geom_line(color = "firebrick", size = 1) + geom_point(color = "firebrick", size = 2) +
  scale_x_reverse() +
  labs(x = "Sequential Model Pair (Start Value)", y = "RMSE",
       title = "Sequential RMSE") +
  theme_minimal()

ggsave("Results/plots/Sequential_RMSE.tiff", rmse_plot,
       width = 8, height = 4, dpi = 300)


# Hellinger distance
hellinger_plot <- ggplot(hellinger_df, aes(x = x, y = y)) +
  geom_line(color = "darkgreen", size = 1) + geom_point(color = "darkgreen", size = 2) +
  scale_x_reverse() +
  labs(x = "Sequential Model Pair (Start Value)", y = "Hellinger Distance",
       title = "Sequential Hellinger Distance") +
  theme_minimal()

ggsave("Results/plots/Sequential_Hellinger.tiff", hellinger_plot,
       width = 8, height = 4, dpi = 300)

# Bhattacharyya coefficient
bhat_plot <- ggplot(bhat_df, aes(x = x, y = y)) +
  geom_line(color = "purple", size = 1) + geom_point(color = "purple", size = 2) +
  scale_x_reverse() +
  labs(x = "Sequential Model Pair (Start Value)", y = "Bhattacharyya Coefficient",
       title = "Sequential Bhattacharyya Coefficient") +
  theme_minimal()

ggsave("Results/plots/Sequential_Bhattacharyya.tiff", bhat_plot,
       width = 8, height = 4, dpi = 300)

## Combined plot
# Scaling factors
scale_rmse <- max(seq_corr) / max(seq_rmse)
scale_hellinger <- max(seq_corr) / max(seq_hellinger)
scale_bhat <- 1  # already 0-1

combined_df <- data.frame(
  x = x_vals,
  Correlation = seq_corr,
  RMSE = seq_rmse * scale_rmse,
  Hellinger = seq_hellinger * scale_hellinger,
  Bhattacharyya = seq_bhat
)

combined_plot <- ggplot(combined_df, aes(x = x)) +
  geom_line(aes(y = Correlation, color = "Correlation"), size = 1) +
  geom_point(aes(y = Correlation, color = "Correlation"), size = 2) +
  geom_line(aes(y = RMSE, color = "RMSE"), size = 1) +
  geom_point(aes(y = RMSE, color = "RMSE"), size = 2) +
  geom_line(aes(y = Hellinger, color = "Hellinger"), size = 1) +
  geom_point(aes(y = Hellinger, color = "Hellinger"), size = 2) +
  geom_line(aes(y = Bhattacharyya, color = "Bhattacharyya"), size = 1) +
  geom_point(aes(y = Bhattacharyya, color = "Bhattacharyya"), size = 2) +
  scale_x_reverse() +
  scale_color_manual(values = c(
    "Correlation" = "steelblue",
    "RMSE" = "firebrick",
    "Hellinger" = "darkgreen",
    "Bhattacharyya" = "purple"
  )) +
  labs(x = "Sequential Model Pair (Start Value)", y = "Scaled Values",
       title = "Sequential Comparison: Correlation, RMSE, Hellinger, Bhattacharyya",
       color = "Metric") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Save combined plot
ggsave("Results/plots/Sequential_All_Metrics.tiff", combined_plot, width = 10, height = 6)


## 5.6 Difference maps ==========================================================

#Create output folder
dir.create(file.path("Results/7_DIFFERENCE_MEAN_STKDE"))

# Sequential differences
diff_list <- lapply(1:(nlyr(final_model_list) - 1), function(i) {
  final_model_list[[i + 1]] - final_model_list[[i]]
})

# Combine into a SpatRaster
diff_rasters <- rast(diff_list)

# Name layers according to the pair
dates <- as.numeric(names(final_model_list))
names(diff_rasters) <- paste0(dates[-1], "-", dates[-length(dates)])

writeRaster(
  diff_rasters,
  filename = file.path("Results/7_DIFFERENCE_MEAN_STKDE", paste0(names(diff_rasters), ".tif")),
  overwrite = TRUE
)

# 6. STKDE BY VARIABLE #########################################################

## 6.1 Creation of STKDE of all variables =======================================

#Output directory
dir.create(file.path("Results/8_VARIABLES_STKDE"))

#Set up the function
stkde_by_material <- function(data, variable, win){
  
  subset_data <- data %>%
    filter(!is.na(.data[[variable]]) & .data[[variable]] == 1)
  
  if(nrow(subset_data) == 0){
    message(paste("Skipping", variable, "- no dates"))
    return(NULL)
  }
  # Correction for site Overrepresentation
  site_counts <- subset_data %>%
    group_by(SITE) %>%
    summarise(n_dates = n(), .groups = "drop")
  
  subset_data <- subset_data %>%
    left_join(site_counts, by = "SITE")
  
  weights <- 1 / subset_data$n_dates
  
  # Output directory
  out_dir <- paste0("Results/8_VARIABLES_STKDE/stkde_", variable)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Calibration
  cal <- rcarbon::calibrate(
    x = subset_data$DATE,
    errors = subset_data$SD,
    calCurves = "intcal20",
    normalised = FALSE,
    verbose = FALSE
  )
  
  # Temporal binning
  bins <- binPrep(
    sites = subset_data$SITE,
    ages = subset_data$DATE,
    h = 50
  )
  
  # Coordinates
  coords <- as.matrix(subset_data[,c("X","Y")])
  
  # Spatio-Temporal KDE
  stkde_res <- stkde(
    x = cal,                         
    coords = coords,                 
    win = win,                       
    sbw = 20000,                     
    tbw = 50,                        
    cellres = 2000,                  
    focalyears = seq(6400, 5000, -100), 
    bins = bins,                    
    backsight = 200,                 
    weights = weights,              
    outdir = out_dir,                
    amount = 0.5,                    
    verbose = FALSE                 
  )
  
  return(stkde_res)
}

# Run STKDE for all variables
stkde_list <- lapply(vars, function(v) stkde_by_material(c14_raw, v, win))

names(stkde_list) <- vars

# Remove variables without results
stkde_list <- stkde_list[!sapply(stkde_list, is.null)]

#Plot the resulting maps
for (v in names(stkde_list)) {
  
dir.create(file.path("Results/8_VARIABLES_STKDE",paste0("stkde_",v),"maps"), showWarnings = FALSE, recursive = TRUE)
  
for(y in seq(6400, 5000, -100)){
  
  load(stkde_list[[as.character(v)]][["impaths"]][[as.character(y)]])
  
  r <- rast(focalyear[["focal"]])             
  
  min_val <- min(values(r), na.rm = TRUE)
  max_val <- max(values(r), na.rm = TRUE)
  
  r_norm <- (r - min_val) / (max_val - min_val)
  
  crs(r_norm) <- crs(ref_raster)
  
  # Reproject to reference CRS
  if (crs(r_norm) != crs(ref_raster)) {
    r_norm <- project(r_norm, ref_raster, method = "bilinear")  # continuous
  }
  
  # Extend to reference extent
  if (!all(ext(r_norm) == ext(ref_raster))) {
    r_norm <- extend(r_norm, ext(ref_raster))
  }
  
  # Match resolution
  if (!all(res(r_norm) == res(ref_raster))) {
    r_norm <- resample(r_norm, ref_raster, method = "bilinear")
  }
  
  r_norm[is.na(r_norm)] <- 0
  
  # Apply mask using the reference raster
  r_norm <- mask(r_norm, ref_raster)
  
  writeRaster(r_norm, paste0("Results/8_VARIABLES_STKDE/",paste0("stkde_",v),"/maps/STKDE_",v,"_",y, ".tiff"), overwrite=TRUE)
  
}

}

## 6.2 Comparison with Final STKDE ==============================================

dir.create(file.path("Results/9_DIFFERENCE_VAR_STKDE"), showWarnings = FALSE, recursive = TRUE)

# Loop over variables
for (v in names(stkde_list)) {
  
  in_dir  <- file.path("Results/8_VARIABLES_STKDE", paste0("stkde_", v), "maps")
  out_dir <- file.path("Results/9_DIFFERENCE_VAR_STKDE", paste0("stkde_", v))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize metrics storage
  metrics_df <- data.frame(
    Date = numeric(),
    RMSE = numeric(),
    Pearson = numeric(),
    Hellinger = numeric(),
    Bhattacharyya = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over dates
  for (y in seq(6400, 5000, -100)) {
    
    # Read STKDE raster
    stkde_path <- file.path(in_dir, paste0("STKDE_", v, "_", y, ".tiff"))
    r_stkde <- rast(stkde_path)
    
    # Get corresponding final model raster
    r_final <- final_model_list[[as.character(y)]]
    
    # Raster subtraction 
    r_diff <- r_final - r_stkde
    
    # Save difference raster
    writeRaster(
      r_diff,
      file.path(out_dir, paste0("DIFF_FINAL_STKDE_", v, "_", y, ".tif")),
      overwrite = TRUE
    )
    
    # Extract values for metrics 
    vals_final <- values(r_final, na.rm = TRUE)
    vals_stkde <- values(r_stkde, na.rm = TRUE)
    
    # RMSE
    rmse_val <- rmse(vals_final, vals_stkde)
    
    # Pearson correlation
    pearson_val <- cor(vals_final, vals_stkde, method = "pearson")
    
    # Normalize for probability distributions
    p <- vals_final / sum(vals_final, na.rm = TRUE)
    q <- vals_stkde / sum(vals_stkde, na.rm = TRUE)
    
    # Bhattacharyya coefficient
    bc_val <- sum(sqrt(p * q), na.rm = TRUE)
    
    # Hellinger distance
    h_val <- sqrt(1 - bc_val)
    
    # Store metrics
    metrics_df <- rbind(metrics_df, data.frame(
      Date = y,
      RMSE = rmse_val,
      Pearson = pearson_val,
      Hellinger = h_val,
      Bhattacharyya = bc_val
    ))
  }
  
  # Save metrics as CSV per variable
  write.csv(metrics_df,
            file.path(out_dir, paste0("Metrics_FINAL_vs_STKDE_", v, ".csv")),
            row.names = FALSE)
}

### 6.2.1 Plot metrics ===========================================================

#Create a folder to store results
out_dir <- "Results/9_DIFFERENCE_VAR_STKDE/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over variables
for (v in names(stkde_list)) {
  
  out_dir <- file.path("Results/9_DIFFERENCE_VAR_STKDE", paste0("stkde_", v))
  
  # Read metrics CSV
  metrics_csv <- file.path(out_dir, paste0("Metrics_FINAL_vs_STKDE_", v, ".csv"))
  metrics_df <- read.csv(metrics_csv)
  
  # Scaling for RMSE and Hellinger
  rmse_max <- max(metrics_df$RMSE, na.rm = TRUE)
  corr_max <- max(metrics_df$ lyr.1, na.rm = TRUE)
  metrics_df$RMSE_scaled <- metrics_df$RMSE * (corr_max / rmse_max)
  
  # Plot 
  combined_plot <- ggplot(metrics_df, aes(x = Date)) +
    geom_line(aes(y = RMSE_scaled, color = "RMSE"), size = 1) +
    geom_point(aes(y = RMSE_scaled, color = "RMSE"), size = 2) +
    
    geom_line(aes(y =  lyr.1, color = "Correlation"), size = 1) +
    geom_point(aes(y =  lyr.1, color = "Correlation"), size = 2) +
    
    geom_line(aes(y = Hellinger, color = "Hellinger"), size = 1) +
    geom_point(aes(y = Hellinger, color = "Hellinger"), size = 2) +
    
    geom_line(aes(y = Bhattacharyya, color = "Bhattacharyya"), size = 1) +
    geom_point(aes(y = Bhattacharyya, color = "Bhattacharyya"), size = 2) +
    
    scale_x_reverse() +
    scale_color_manual(values = c(
      "RMSE" = "firebrick",
      "Correlation" = "steelblue",
      "Hellinger" = "darkgreen",
      "Bhattacharyya" = "purple"
    )) +
    labs(
      title = paste("Comparison of Metrics for Variable:", v),
      x = "Date (Start Value)",
      y = "Metric Value (scaled for visualization)",
      color = "Metric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  # Save plot 
  ggsave(
    filename = file.path("Results/9_DIFFERENCE_VAR_STKDE/plots",paste0("Combined_Metrics_", v, ".tiff")),
    plot = combined_plot,
    width = 10, height = 5, dpi = 300
  )
}

## 6.3 Plot the comparison graphics =============================================

#Set up out_dir again
out_dir <- "Results/9_DIFFERENCE_VAR_STKDE/plots"

vars <- names(stkde_list)

for (y in names(final_model_list)) {
  
  
  r_final <- final_model_list[[y]]
  
  diff_list <- lapply(vars, function(v){
    rast(file.path("Results/9_DIFFERENCE_VAR_STKDE",
                   paste0("stkde_", v),
                   paste0("DIFF_FINAL_STKDE_", v, "_", y, ".tif")))
  })
  
  diff_vals <- unlist(lapply(diff_list, function(r) {
    vals <- values(r, na.rm = TRUE)
    if(length(vals) > 0) vals else 0
  }))
  diff_max <- max(abs(diff_vals), na.rm = TRUE)
  
  png(filename = file.path(out_dir, paste0("FINAL_DIFF_", y, ".png")),
      width = 1800, 
      height = 1400 + 700*ceiling(length(vars)/2))
  
  par(oma = c(2, 2, 2, 2))
  
  n_diff <- length(diff_list)
  n_cols <- 2
  n_rows <- ceiling(n_diff / n_cols)
  layout_matrix <- rbind(matrix(1, nrow = 1, ncol = n_cols),
                         matrix(2:(n_diff+1), nrow = n_rows, ncol = n_cols, byrow = TRUE))
  layout(layout_matrix)
  
  par(mar = c(4, 4, 2, 2))
  
  n_colors <- 100
  col_final <- viridis(n_colors)
  
  # FINAL MODEL
  plot(r_final, col = col_final)
  
  e <- ext(r_final)
  x_center <- (xmin(e) + xmax(e)) / 2
  y_bottom <- ymin(e) + 0.05 * (ymax(e) - ymin(e))
  
  text(x_center, y_bottom,
       labels = paste("Final Model", y,"Cal BP"),
       cex = 5,
       font = 2)
  
  # DIFF RASTERS
  for (i in seq_along(diff_list)) {
    
    r <- diff_list[[i]]
    
    plot(r,
         col = col_div,
         zlim = c(-diff_max, diff_max))
    
    e <- ext(r)
    x_center <- (xmin(e) + xmax(e)) / 2
    y_bottom <- ymin(e) + 0.05 * (ymax(e) - ymin(e))
    
    text(x_center, y_bottom,
         labels = paste("Diff", vars[i]),
         cex = 5,
         font = 2)
  }
  
  dev.off()
}

# 7. FINAL MODEL CONTRUCTION  ##################################################

## 7.1 Probabilistic chronological assignment of undated sites  ==================

#Load final models
folder_models <- "Results/6_MEAN_CORRECTED_STKDE"
final_model_list <- lapply(list.files(folder_models, pattern = "\\.tif[f]?$", full.names = TRUE), rast)
final_model_list <- rast(final_model_list)
names(final_model_list) <- seq(5000, 6400, 100) #Assign names

#Extract values of probability of density of No C14 sites from final models
chrono_values <- terra::extract(final_model_list, no_c14_sf)

## Normalize values
# Save IDs and erase them from values (IDs refer to number of row not to IDs of site)
ids <- chrono_values$ID
prob_mat <- chrono_values[, -1]

prob_norm <- prob_mat / apply(prob_mat, 1, max, na.rm = TRUE)
chrono_norm <- cbind(ID = ids, prob_norm)

## Binary classification per raster using row-wise GMM 
# Function to classify a single row
binary_gmm_row <- function(p_row) {
  p_row <- as.numeric(p_row)
  if(length(p_row[!is.na(p_row)]) < 2) return(rep(NA_integer_, length(p_row)))
  
  model <- Mclust(p_row, G = 2, verbose = FALSE)
  
  # Component with higher mean = "occupation"
  comp_occup <- which.max(model$parameters$mean)
  
  # Binary classification for each raster in this row
  as.integer(model$classification == comp_occup)
}

# Apply across all rows
binary_matrix <- t(apply(chrono_norm[, -1], 1, binary_gmm_row))


# Assign raster names as column names
colnames(binary_matrix) <- names(final_model_list)

# Convert matrix to data frame
binary_df <- as.data.frame(binary_matrix)

# Add ID column (matching row order)
binary_df$ID <- seq_len(nrow(binary_df))

# Join
no_c14_sf_final <- no_c14_sf %>%
  mutate(ID = row_number()) %>%
  left_join(binary_df, by = "ID")

## 7.2 Assignment of C14-dated sites ============================================

# Convert raw C14 dates to integer vector
ages <- as.integer(c14_raw$DATE)

# Compute min and max calibrated dates for each sample
max_dates <- numeric(length(ages))
min_dates <- numeric(length(ages))
sample_ids <- C14_raw_calibration$metadata$DateID

for(i in seq_along(ages)) {
  cal_grid <- C14_raw_calibration[i]$grids[[as.character(i)]]$calBP
  max_dates[i] <- max(cal_grid)   # Maximum of probability distribution
  min_dates[i] <- min(cal_grid)   # Minimum of probability distribution
}

# Create data frame of calibrated dates
calibrated_dates <- data.frame(
  id = sample_ids,
  max_date = max_dates,
  min_date = min_dates
)

# Create temporal windows
intervals <- seq(6450, 4950, -100)  # From 6450 to 4950 in steps of 100
wsize <- 100
error <- as.integer(wsize*0.50)

# Compute midpoints of consecutive intervals
midpoints <- (intervals[-length(intervals)] + intervals[-1]) / 2

# Initialize vectors for window assignment
assigned_ids <- c()
assigned_window_mid <- c()

# Loop over all windows
for (j in seq_along(midpoints)) {
  start_window <- intervals[j]
  end_window <- intervals[j + 1]
  
  # Filter dates within or overlapping the window
  candidates <- filter(
    calibrated_dates,
    abs(end_window - min_date) <= wsize |
      abs(start_window - max_date) <= wsize |
      (max_date > start_window & min_date < end_window)
  )
  
  # Further filter based on error threshold
  valid <- filter(
    candidates,
    (abs(max_date - end_window) >= error & max_date <= start_window) |
      (abs(start_window - min_date) >= error & min_date >= end_window) |
      (max_date > start_window & min_date < end_window)
  )
  
  if (nrow(valid) > 0) {
    assigned_ids <- c(assigned_ids, valid$id)
    assigned_window_mid <- c(assigned_window_mid, rep(midpoints[j], nrow(valid)))
  }
}

# Create presence/absence data frame
window_data <- data.frame(
  id = assigned_ids,
  window_mid = assigned_window_mid
)

# Convert to wide format
presence_matrix <- dcast(window_data, id ~ window_mid, fun.aggregate = length)

# Convert counts to 0/1 presence
presence_matrix[,-1] <- ifelse(presence_matrix[,-1] > 0, 1, 0)

# Set rownames and remove id column
rownames(presence_matrix) <- presence_matrix$id
presence_matrix <- presence_matrix[,-1]

# Convert rownames of presence_matrix into a column for merging
presence_matrix$ID <- rownames(presence_matrix)

# Merge presence_matrix into c14_sf
c14_sf <- merge(c14_sf, presence_matrix, by = "ID", all.x = TRUE)

## 7.3 Construct and plot final models ==========================================

#crate a folder to store
dir.create("Results/10_FINAL_MODELS")
dir.create("Results/10_FINAL_MODELS/maps")
dir.create("Results/10_FINAL_MODELS/shapes")

## Merge sf
# Add a source identifier to each dataset
c14_sf$source <- "c14"
no_c14_sf_final$source <- "no_c14"

# Combine datasets
final_model <- rbind(
  c14_sf[, c("ESTRUCTURE","5000", "5100", "5200", "5300", "5400" ,"5500" ,"5600", "5700", "5800", "5900" ,"6000" ,"6100" ,"6200" ,"6300" ,"6400", "geometry", "source")],
  no_c14_sf_final[, c("ESTRUCTURE","5000", "5100", "5200", "5300", "5400" ,"5500" ,"5600", "5700", "5800", "5900" ,"6000" ,"6100" ,"6200" ,"6300" ,"6400", "geometry", "source")]
)

# Columns with raster values
raster_cols <- as.character(seq(6400, 5000, -100))

# Define colors for each source
source_colors <- c("c14" = "red", "no_c14" = "#FF7F00")

# Loop over each raster column
for(y in raster_cols){
  
  # Select points with occupancy = 1
  occupied_points <- final_model[final_model[[y]] == 1, ]
  
  # Save shapefile for all occupied points
  st_write(occupied_points, paste0("Results/10_FINAL_MODELS/shapes/points_in_", y, ".shp"), delete_layer = TRUE)
  
  # Save plot
  png(filename = paste0("Results/10_FINAL_MODELS/maps/", y, ".png"), width = 800, height = 600)
  
  # Plot the raster first
  plot(ref_raster, main = paste(y, "cal BP"))

  # Overlay points by source with different colors
  for(src in unique(occupied_points$source)){
    plot(
      st_geometry(occupied_points[occupied_points$source == src, ]),
      col = source_colors[src],
      pch = 16,
      add = TRUE
    )
  }
  
  # Optional: Add legend
  legend(x = par("usr")[1] + 0.75 * diff(par("usr")[1:2]), 
         y = par("usr")[3] + 0.10 * diff(par("usr")[3:4]),
         legend = c("C14 sites", "No C14 sites"), col = c("red", "#FF7F00"), pch = 16, bty = "n")
  
  dev.off()
}

## GIFF creation
# Folder containing rasters
folder_maps <- "Results/10_FINAL_MODELS/maps"

# Read PNG files (sorted!)
files <- list.files(folder_maps, pattern = "\\.png$", full.names = TRUE)
files <- sort(files, decreasing = TRUE)

# Read images
img_list <- lapply(seq_along(files), function(i) {
  img <- image_read(files[i])
  
  img
})

# Create GIF
gif <- image_animate(image_join(img_list), fps = 0.5)

# Save GIF
image_write(gif, "Results/10_FINAL_MODELS/Final_model.gif")

## 7.4 Summary graphic ==========================================================

## Graph
# Convert to long format
df_long <- final_model %>%
  st_drop_geometry() %>%  # remove geometry
  pivot_longer(
    cols = all_of(raster_cols),
    names_to = "year",
    values_to = "presence"
  ) %>%
  filter(presence == 1) %>%  # keep only presences
  group_by(year, source) %>%
  summarise(freq = n(), .groups = "drop")

# Ensure proper ordering
df_long$year <- as.numeric(df_long$year)

# Bar plot
summ_graph <- ggplot(df_long, aes(x = year, y = freq, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = source_colors) +
  labs(
    x = "Year (cal BP)",
    y = "Number of sites",
    fill = "Source"
  ) +
  theme_minimal()

#Save the graph
ggsave("Results/10_FINAL_MODELS/summ_graph.png", summ_graph)

## Table

# Ensure proper ordering
df_long$year <- as.numeric(df_long$year)

df_wide <- df_long %>%
  tidyr::pivot_wider(
    names_from = source,
    values_from = freq,
    values_fill = 0
  )

# Export to CSV 
write.csv(df_wide, "Results/10_FINAL_MODELS/site_frequency_by_year_and_source.csv", row.names = FALSE)
