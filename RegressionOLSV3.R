# RegressionOLSV3.R
# Run from the repository root. Requires: tidyverse, broom, readr, fs, here, stringr, car, jsonlite
# - Self-bootstraps CSVs from:
#   https://github.com/Remus753/climate-policy-regressions/tree/main/Datasets
#   https://github.com/Remus753/climate-policy-regressions/tree/main/Policies
# - Only runs 2-year lagged effects (no contemporaneous policy terms)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

library(tidyverse)
library(broom)
library(readr)
library(fs)
library(stringr)
library(here)
library(car)        # for VIFs
library(jsonlite)   # for GitHub Contents API

# --- Repo + paths ---
GITHUB_REPO   <- "Remus753/climate-policy-regressions"
base_dir      <- here::here("Datasets")
policies_dir  <- here::here("Policies")
output_dir    <- here::here("Results")

COUNTRY_ROW_START <- 3
COUNTRY_ROW_END   <- 38
YEAR_START        <- 2000
YEAR_END          <- 2023

paths <- list(
  carbon_pricing = file.path(base_dir, "Climate National Indicators EU-BRICS - Carbon Pricing.csv"),
  co2_pc         = file.path(base_dir, "Climate National Indicators EU-BRICS - CO2 Emission Per Capita.csv"),
  fossil_pct     = file.path(base_dir, "Climate National Indicators EU-BRICS - Fossil Fuel %.csv"),
  gdp_pc         = file.path(base_dir, "Climate National Indicators EU-BRICS - GDP per Capita.csv"),
  population     = file.path(base_dir, "Climate National Indicators EU-BRICS - Population.csv"),
  climate_policy = file.path(base_dir, "Climate Policy EU-BRICS - Climate Policy.csv")
)

# --- Utilities: bootstrap CSVs from GitHub if missing ---
github_list_csvs <- function(repo, subdir) {
  url <- paste0("https://api.github.com/repos/", repo, "/contents/", subdir)
  resp <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  if (is.null(resp)) return(tibble())
  as_tibble(resp) |>
    filter(type == "file", grepl("\\.csv$", name, ignore.case = TRUE)) |>
    transmute(name, download_url)
}

ensure_dir <- function(p) if (!dir_exists(p)) dir_create(p)

download_if_missing <- function(url, dest) {
  if (is.na(url) || is.na(dest) || url == "" || dest == "") return(invisible(FALSE))
  if (!file_exists(dest)) {
    message("Downloading: ", basename(dest))
    try(download.file(url, destfile = dest, mode = "wb", quiet = TRUE), silent = TRUE)
  }
  file_exists(dest)
}

bootstrap_folder_from_github <- function(repo, subdir, dest_dir) {
  ensure_dir(dest_dir)
  listing <- github_list_csvs(repo, subdir)
  if (nrow(listing) == 0) {
    message("Warning: couldn't list ", subdir, " via API; will fall back to known files if provided.")
  }
  ok <- logical(0)
  if (nrow(listing) > 0) {
    ok <- vapply(seq_len(nrow(listing)), function(i) {
      dest <- file.path(dest_dir, listing$name[i])
      download_if_missing(listing$download_url[i], dest)
    }, logical(1))
  }
  invisible(any(ok) || nrow(listing) > 0)
}

# Fallback helper to ensure specific files (for Datasets folder names used below)
ensure_known_files <- function(repo, subdir, dest_dir, files) {
  ensure_dir(dest_dir)
  for (f in files) {
    dest <- file.path(dest_dir, f)
    if (!file_exists(dest)) {
      url <- paste0("https://raw.githubusercontent.com/", repo, "/main/", subdir, "/", utils::URLencode(f))
      download_if_missing(url, dest)
    }
  }
}

# --- Bootstrap data: Datasets + Policies ---
ensure_dir(output_dir)

# Datasets via API (and fallback to the 6 known filenames)
bootstrap_folder_from_github(GITHUB_REPO, "Datasets", base_dir)
ensure_known_files(
  GITHUB_REPO, "Datasets", base_dir,
  c(
    "Climate National Indicators EU-BRICS - Carbon Pricing.csv",
    "Climate National Indicators EU-BRICS - CO2 Emission Per Capita.csv",
    "Climate National Indicators EU-BRICS - Fossil Fuel %.csv",
    "Climate National Indicators EU-BRICS - GDP per Capita.csv",
    "Climate National Indicators EU-BRICS - Population.csv",
    "Climate Policy EU-BRICS - Climate Policy.csv"
  )
)

# Policies via API (download all CSVs present)
bootstrap_folder_from_github(GITHUB_REPO, "Policies", policies_dir)

# --- Helpers for loading/labeling ---
load_matrix_csv <- function(path, value_name) {
  stopifnot(file_exists(path))
  raw <- read_csv(path, col_names = FALSE, show_col_types = FALSE)
  header_years <- raw[1, -1] |> as.character() |> as.integer()
  keep_idx <- which(header_years >= YEAR_START & header_years <= YEAR_END)
  years <- header_years[keep_idx]
  countries <- raw[COUNTRY_ROW_START:COUNTRY_ROW_END, 1, drop=FALSE] |> pull(1) |> as.character()
  values_block <- raw[COUNTRY_ROW_START:COUNTRY_ROW_END, keep_idx + 1, drop=FALSE]
  values_mat <- suppressWarnings(apply(as.matrix(values_block), 2, as.numeric))
  tibble(
    country = rep(countries, times = length(years)),
    year    = rep(years, each = length(countries)),
    !!value_name := as.vector(values_mat)
  ) |> filter(!is.na(country), !is.na(year))
}

col_label_to_index <- function(label) {
  s <- toupper(trimws(label))
  if (!grepl("^[A-Z]+$", s)) return(NA_integer_)
  n <- 0L
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    n <- n * 26L + (utf8ToInt(ch) - utf8ToInt("A") + 1L)
  }
  n
}

get_policy_name <- function(code) {
  idx <- col_label_to_index(code)
  if (is.na(idx) || idx > ncol(climate_policy_meta)) return(tolower(code))
  val <- climate_policy_meta[1, idx, drop = TRUE] |> as.character() |> trimws()
  if (is.na(val) || val == "") tolower(code) else val
}

get_category <- function(code) {
  idx <- col_label_to_index(code)
  if (is.na(idx)) return(NA_character_)
  if (idx >= col_label_to_index("Z")  && idx <= col_label_to_index("CS")) return("Policy Instrument")
  if (idx >= col_label_to_index("CU") && idx <= col_label_to_index("CZ")) return("Sector")
  if (idx >= col_label_to_index("DB") && idx <= col_label_to_index("DG")) return("Policy Type (Mitigation Area)")
  if (idx >= col_label_to_index("DI") && idx <= col_label_to_index("DQ")) return("Policy Objective")
  NA_character_
}

build_two_letter_seq <- function(start = "aa", end = "dq") {
  lv <- unlist(lapply(letters, function(x) paste0(x, letters)))
  si <- match(start, lv); ei <- match(end, lv)
  if (is.na(si) || is.na(ei) || ei < si) return(character(0))
  lv[si:ei]
}

safe_var_name <- function(path) {
  nm <- path_file(path) |> str_remove("\\.csv$") |> str_replace_all("\u00A0", " ") |> str_squish()
  m <- str_match(nm, "(?i)policy\\s*dummy[^A-Za-z]*([A-Za-z]{1,2})\\b")[,2]
  tolower(m)
}

# --- Load data ---
co2_long    <- load_matrix_csv(paths$co2_pc, "co2_pc")
gdp_pc_long <- load_matrix_csv(paths$gdp_pc, "gdp_pc")
pop_long    <- load_matrix_csv(paths$population, "population")
fossil_long <- load_matrix_csv(paths$fossil_pct, "fossil_pct")
cp_long     <- load_matrix_csv(paths$carbon_pricing, "carbon_pricing")

climate_policy_meta <- read_csv(paths$climate_policy, col_names = FALSE, show_col_types = FALSE)

cp_dummy <- cp_long |> mutate(carbon_pricing_dummy = as.integer(replace_na(carbon_pricing, 0) > 0)) |>
  select(country, year, carbon_pricing_dummy)

panel_base <- co2_long |>
  left_join(gdp_pc_long, by = c("country","year")) |>
  left_join(pop_long,    by = c("country","year")) |>
  left_join(fossil_long, by = c("country","year")) |>
  left_join(cp_dummy,    by = c("country","year")) |>
  mutate(
    ln_co2_pc  = log(co2_pc + 1e-9),
    ln_gdp_pc  = log(gdp_pc + 1e-9),
    ln_pop     = log(population + 1e-9),
    ln_co2_pc_norm   = scale(ln_co2_pc)[,1],
    ln_gdp_pc_norm   = scale(ln_gdp_pc)[,1],
    ln_pop_norm      = scale(ln_pop)[,1],
    fossil_pct_norm  = scale(fossil_pct)[,1]
  ) |>
  filter(year >= YEAR_START, year <= YEAR_END)

# Ensure outputs dir
ensure_dir(output_dir)

# --- Multicollinearity diagnostics for controls only ---
fit_controls <- lm(
  ln_co2_pc_norm ~ ln_gdp_pc_norm + ln_pop_norm + fossil_pct_norm + carbon_pricing_dummy,
  data = panel_base
)

vif_vals <- car::vif(fit_controls)
controls_vif_tbl <- tibble(
  variable  = names(vif_vals),
  vif       = as.numeric(vif_vals),
  tolerance = 1 / as.numeric(vif_vals)
)
write_csv(controls_vif_tbl, file.path(output_dir, "controls_multicollinearity_vif.csv"))

# Correlation matrices (Excel-style)
controls_mat <- panel_base |>
  select(ln_gdp_pc_norm, ln_pop_norm, fossil_pct_norm, carbon_pricing_dummy) |>
  as.data.frame()
controls_cor <- suppressWarnings(cor(controls_mat, use = "pairwise.complete.obs"))

with_outcome_mat <- panel_base |>
  select(ln_co2_pc_norm, ln_gdp_pc_norm, ln_pop_norm, fossil_pct_norm, carbon_pricing_dummy) |>
  as.data.frame()
controls_with_outcome_cor <- suppressWarnings(cor(with_outcome_mat, use = "pairwise.complete.obs"))

round_to <- 3

# Controls-only full
controls_cor_df <- data.frame(var = rownames(controls_cor),
                              round(controls_cor, round_to),
                              check.names = FALSE)
write.csv(controls_cor_df, file.path(output_dir, "controls_correlations_matrix.csv"),
          row.names = FALSE, na = "")

# Controls-only lower triangle
controls_lower <- round(controls_cor, round_to)
controls_lower[upper.tri(controls_lower, diag = FALSE)] <- NA
controls_lower_df <- data.frame(var = rownames(controls_lower),
                                controls_lower,
                                check.names = FALSE)
controls_lower_df[is.na(controls_lower_df)] <- ""
write.csv(controls_lower_df,
          file.path(output_dir, "controls_correlations_matrix_lowertriangle.csv"),
          row.names = FALSE, na = "")

# Controls + outcome full
with_outcome_df <- data.frame(var = rownames(controls_with_outcome_cor),
                              round(controls_with_outcome_cor, round_to),
                              check.names = FALSE)
write.csv(with_outcome_df,
          file.path(output_dir, "controls_with_outcome_correlations_matrix.csv"),
          row.names = FALSE, na = "")

# Controls + outcome lower triangle
with_outcome_lower <- round(controls_with_outcome_cor, round_to)
with_outcome_lower[upper.tri(with_outcome_lower, diag = FALSE)] <- NA
with_outcome_lower_df <- data.frame(var = rownames(with_outcome_lower),
                                    with_outcome_lower,
                                    check.names = FALSE)
with_outcome_lower_df[is.na(with_outcome_lower_df)] <- ""
write.csv(with_outcome_lower_df,
          file.path(output_dir, "controls_with_outcome_correlations_matrix_lowertriangle.csv"),
          row.names = FALSE, na = "")

message("Saved controls VIF and Excel-style correlation matrices.")

# --- Policy regressions (lag-2 effects ONLY) ---
stopifnot(dir_exists(policies_dir))
policy_files <- dir_ls(policies_dir, glob = "*.csv")
if (length(policy_files) == 0) stop("No policy CSVs found in the Policies folder: ", policies_dir)

allowed_codes <- c("z", build_two_letter_seq("aa","dq")) |> setdiff(c("po","ct","da","dh","ai"))

policy_codes_all  <- map_chr(policy_files, safe_var_name)
keep_idx          <- which(!is.na(policy_codes_all) & policy_codes_all %in% allowed_codes)
policy_files_kept <- policy_files[keep_idx]
policy_codes_kept <- policy_codes_all[keep_idx]

all_results <- list()
all_glance  <- list()

for (i in seq_along(policy_files_kept)) {
  pf <- policy_files_kept[i]
  var_code <- policy_codes_kept[i]
  var_name <- get_policy_name(var_code)
  cat_name <- get_category(var_code)
  message("Processing policy: ", var_code, " (", var_name, ") — Category: ", cat_name)
  
  pol_long <- load_matrix_csv(pf, var_code) |>
    mutate(!!var_code := as.integer(replace_na(.data[[var_code]], 0) != 0))
  
  df <- panel_base |>
    left_join(pol_long, by = c("country","year")) |>
    mutate(
      carbon_pricing_dummy = as.integer(replace_na(carbon_pricing_dummy, 0)),
      policy_var           = as.integer(replace_na(.data[[var_code]], 0))
    ) |>
    group_by(country) |>
    arrange(year, .by_group = TRUE) |>
    mutate(policy_var_lag2 = lag(policy_var, 2)) |>
    ungroup()
  
  # Lag-2 only
  fit <- lm(ln_co2_pc_norm ~ policy_var_lag2 +
              ln_gdp_pc_norm + ln_pop_norm + fossil_pct_norm + carbon_pricing_dummy,
            data = df)
  
  tid <- tidy(fit) |>
    mutate(policy_code = var_code, policy_name = var_name, category = cat_name)
  glc <- glance(fit) |>
    mutate(policy_code = var_code, policy_name = var_name, category = cat_name, N = nrow(model.frame(fit)))
  
  all_results[[var_code]] <- tid
  all_glance[[var_code]]  <- glc
}

results_tbl <- bind_rows(all_results)

summary_tbl <- bind_rows(all_glance) |>
  select(policy_code, policy_name, category, r.squared, adj.r.squared, sigma, AIC, BIC, N)

policy_effects <- results_tbl |>
  filter(term %in% c("policy_var_lag2")) |>
  mutate(
    effect     = "lag2",
    estimate   = estimate,
    std_error  = std.error,
    statistic  = statistic,
    p_value    = p.value
  ) |>
  select(policy_code, policy_name, category, effect, term, estimate, std_error, statistic, p_value) |>
  left_join(summary_tbl, by = c("policy_code","policy_name","category"))

policy_effects_ranked <- policy_effects |>
  mutate(
    rank_all = min_rank(estimate),
    rank_sig = if_else(p_value < 0.05, min_rank(estimate), NA_integer_)
  ) |>
  group_by(category) |>
  mutate(
    rank_cat_all = min_rank(estimate),
    rank_cat_sig = if_else(p_value < 0.05, min_rank(estimate), NA_integer_)
  ) |>
  ungroup() |>
  arrange(estimate, policy_code)

policy_levels <- c("z", build_two_letter_seq("aa","dq"))
comparison_table <- policy_effects |>
  mutate(policy_code = factor(policy_code, levels = policy_levels)) |>
  arrange(policy_code)

# --- Save outputs (lag-2 only) ---
write_csv(results_tbl,      file.path(output_dir, "policy_regression_full_coefs.csv"))
write_csv(summary_tbl,      file.path(output_dir, "policy_regression_model_summaries.csv"))
write_csv(comparison_table, file.path(output_dir, "policy_regression_lag2_only_table.csv"))
write_csv(policy_effects_ranked, file.path(output_dir, "policy_ranks_lag2.csv"))

message("Lag-2 only outputs saved in ", output_dir)
