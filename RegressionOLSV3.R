# RegressionOLSV3.R
# Run from the repository root. Requires: tidyverse, broom, readr, fs, here, stringr

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

library(tidyverse)
library(broom)
library(readr)
library(fs)
library(stringr)
library(here)

# --- File paths (relative to repo root) ---
base_dir     <- here::here("Datasets")
policies_dir <- here::here("Policies")
output_dir   <- here::here("Results")

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

# --- Helpers ---
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
  
  fit <- lm(ln_co2_pc_norm ~ policy_var + policy_var_lag2 + ln_gdp_pc_norm + ln_pop_norm + fossil_pct_norm + carbon_pricing_dummy, data = df)
  
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

policy_coefs <- results_tbl |>
  filter(term %in% c("policy_var","policy_var_lag2")) |>
  transmute(
    policy_code,
    policy_name,
    category,
    term,
    estimate,
    std_error = std.error,
    statistic,
    p_value = p.value
  )

policy_levels <- c("z", build_two_letter_seq("aa","dq"))

comparison_table <- policy_coefs |>
  left_join(summary_tbl, by = c("policy_code","policy_name","category")) |>
  mutate(policy_code = factor(policy_code, levels = policy_levels)) |>
  arrange(policy_code, term)

if (!dir_exists(output_dir)) dir_create(output_dir)

policy_corr_list <- list()
policy_corr_lag2_list <- list()

for (i in seq_along(policy_files_kept)) {
  pf <- policy_files_kept[i]
  var_code <- policy_codes_kept[i]
  var_name <- get_policy_name(var_code)
  cat_name <- get_category(var_code)
  
  pol_long <- load_matrix_csv(pf, var_code) |>
    mutate(!!var_code := as.integer(replace_na(.data[[var_code]], 0) != 0))
  
  df <- panel_base |>
    left_join(pol_long, by = c("country","year")) |>
    mutate(policy_var = as.integer(replace_na(.data[[var_code]], 0))) |>
    group_by(country) |>
    arrange(year, .by_group = TRUE) |>
    mutate(policy_var_lag2 = lag(policy_var, 2)) |>
    ungroup()
  
  n_now  <- sum(complete.cases(df$ln_co2_pc_norm, df$policy_var))
  n_lag2 <- sum(complete.cases(df$ln_co2_pc_norm, df$policy_var_lag2))
  
  c_now  <- suppressWarnings(cor(df$ln_co2_pc_norm, df$policy_var, use = "complete.obs"))
  c_lag2 <- suppressWarnings(cor(df$ln_co2_pc_norm, df$policy_var_lag2, use = "complete.obs"))
  
  policy_corr_list[[var_code]]    <- tibble(policy_code = var_code, policy_name = var_name, category = cat_name, correlation = c_now,  n_pairs = n_now)
  policy_corr_lag2_list[[var_code]] <- tibble(policy_code = var_code, policy_name = var_name, category = cat_name, correlation = c_lag2, n_pairs = n_lag2)
}

policy_corr      <- bind_rows(policy_corr_list)      |> arrange(correlation, policy_code) |> mutate(rank = row_number())
policy_corr_lag2 <- bind_rows(policy_corr_lag2_list) |> arrange(correlation, policy_code) |> mutate(rank = row_number())

write_csv(results_tbl,        file.path(output_dir, "policy_regression_full_coefs.csv"))
write_csv(summary_tbl,        file.path(output_dir, "policy_regression_model_summaries.csv"))
write_csv(comparison_table,   file.path(output_dir, "policy_regression_comparison_table.csv"))
write_csv(policy_corr,        file.path(output_dir, "policy_correlations_contemporaneous.csv"))
write_csv(policy_corr_lag2,   file.path(output_dir, "policy_correlations_lag2.csv"))

message("Saved: policy_regression_full_coefs.csv, policy_regression_model_summaries.csv, policy_regression_comparison_table.csv, policy_correlations_contemporaneous.csv, policy_correlations_lag2.csv in ", output_dir)

comparison_table |> slice_head(n = 10) |> print(n = 10)
