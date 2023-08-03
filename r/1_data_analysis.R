
# ---- libs ----


suppressPackageStartupMessages({
  library("readr")
  library("dplyr")
  library("tidyr")
  library("forcats")
  library("purrr")
  library("furrr")
  library("lubridate") # way to handle dates better than default R way
  library("tictoc")    # measure time elapsed in calcs
  library("ggplot2") 
  library("ggrepel") 
  library("knitr")
  library("gsDesign")
  library("foreach")
  library("arrow") # read/write parquet files
})


# NOTE: need to run first (only once, assumes devtools installed):
# devtools::install_github("tystan/pharmsignal") 
library("pharmsignal") # signal detection algs

# here are the functions written for these analyses
# they will be shown in the *Appendix A*
source("r/_funcs.R")



# ---- check_parallel_comp ----

# options(future.globals.maxSize = 500 * 1024 ^ 2) # = 500 MiB
options(future.globals.maxSize = 1e3 * 1024 ^ 2) # = 1 GB



# furrr parallel workers/cores setup
# change `workers = 4` based on cores available in processor being used
plan(multisession, workers = 4) 

### test parallel works
# test code from https://furrr.futureverse.org/
# sequential
tic()
dev_null <- map(c(2, 2, 2), ~Sys.sleep(.x))
toc() # ~6 sec
# parallel: should be (roughly, plus overheads) a third of the time of sequential
tic()
dev_null <- future_map(c(2, 2, 2), ~Sys.sleep(.x))
toc() # ~2 sec


# this only applies to the non-parallel (non-"future") operations
set.seed(1234) 
# this seed can be set in future_map() etc for reproducible parallel comp seeds 
furrr_seed1 <- furrr_options(seed = 5678)
furrr_seed2 <- furrr_options(seed = 9012)
furrr_seed3 <- furrr_options(seed = 3456)


# ---- consts ----

# arbitrarily, let's go with minimum cell count of 3 (should be discussed!)
arbitrary_cell_min <- 1



# ---- funcs ----


get_sig_tab <- function(nA, nB, nC, nD, alpha = 0.05, n_mcmc = 1e+05) {
  
  out_cols_of_interest <- c("est_name", "est_scale", "est", "alpha", "ci_lo", "ci_hi")
  sig_tab <- pharmsignal::bcpnn_mcmc_signal(nA, nB, nC, nD, alpha = alpha, n_mcmc = n_mcmc)
  sig_tab <- sig_tab[, out_cols_of_interest]
  # sig_tab <- bind_cols(tibble(mnth = mnth), sig_tab)
  return(sig_tab)
  
}

get_sig_tab_over_time <- function(dat, alpha = 0.05, n_mcmc = 1e+05) {
  
  n_tp <- nrow(dat)
  
  sig_tab_over_time <-
    foreach(i = 1:n_tp, .combine = bind_rows, .packages = "dplyr") %do% {
      with(
        dat, 
        get_sig_tab(
          # mnth[i], 
          nA[i], nB[i], nC[i], nD[i], 
          alpha = alpha, n_mcmc = n_mcmc
        )
      )
    }
  
  return(sig_tab_over_time)
  
}

# same as get_sig_tab_over_time(), however, alpha assumed included as column in data
get_sig_tab_over_time_2 <- function(dat, n_mcmc = 1e+05) {
  
  n_tp <- nrow(dat)
  
  sig_tab_over_time <-
    foreach(i = 1:n_tp, .combine = bind_rows, .packages = "dplyr") %do% {
      with(
        dat, 
        get_sig_tab(
          # mnth[i], 
          nA[i], nB[i], nC[i], nD[i], 
          alpha = adj_alpha[i], 
          n_mcmc = n_mcmc
        )
      )
    }
  
  return(sig_tab_over_time)
  
}


# ---- load_dat ----



sra_dat <- read_parquet("dat/sra_dat.parquet")
cumul_qtrly_dat <- read_parquet("dat/cumul_qtrly_dat.parquet")

thresholds <- sort(unique(sra_dat$thresh))



# ---- bcpnn_calcs ----






# sra_cum <- 
#   sra_dat %>%
#   dplyr::filter(dat_type == "cumulative") 
sra_cum <- 
  cumul_qtrly_dat

# make data for each combination of params nested for purrr like processing
sra_cum <-
  sra_cum %>%
  nest(data = c(mnth, nA, nB, nC, nD))

# testing/example
sra_cum$data[[1]]
get_sig_tab_over_time(sra_cum$data[[1]])




### takes ~ 90 sec
tic()
sra_cum <-
  sra_cum %>%
  mutate(
    sig_tab = 
      future_map(
        .x = data, 
        .f = get_sig_tab_over_time,
        .options = furrr_seed1
      )
  )
toc()

# check
sra_cum$sig_tab[[1]]


sra_cum_bcpnn <-
  sra_cum %>%
  unnest(cols = c(data, sig_tab)) %>%
  mutate(dte = as_date(paste0(mnth, "-01")))

sra_cum_bcpnn


# first signif
bcpnn_signif <-
  sra_cum_bcpnn %>%
  group_by(grps, dat_type, thresh) %>%
  dplyr::filter(ci_lo > 0) %>%
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(sra_cum_bcpnn)
sra_cum_bcpnn <-
  left_join(
    sra_cum_bcpnn,
    bcpnn_signif %>% select(grps, dat_type, thresh, dte_reach_sig),
    c("grps", "dat_type", "thresh")
  )
nrow(sra_cum_bcpnn)

sra_cum_bcpnn


sra_cum_bcpnn <- 
  sra_cum_bcpnn %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )



# ---- save1 ----


sra_cum_bcpnn %>%
  write_parquet(., sink = "out/sra_cum_bcpnn.parquet")







# ---- multcompar ----

# sra_cum <- 
#   sra_dat %>%
#   dplyr::filter(dat_type == "cumulative") 
sra_cum <- 
  cumul_qtrly_dat


sra_cum <-
  sra_cum %>%
  nest(data = c(mnth, nA, nB, nC, nD))


# if it's multiple comparisons central need to sparing use alpha
get_mult_compare_adj_alpha <- function(dat) {
  
  n_reports <- nrow(dat)
  
  information_fracs <-  1:n_reports / n_reports
  
  ### alternatives:
  # spend_obj <- sfLDPocock(alpha = 0.025, t = information_fracs, param = NULL)
  # spend_obj <- sfLDOF(alpha = 0.025, t = information_fracs, param = NULL)
  spend_obj <- sfExponential(alpha = 0.05, t = information_fracs, param = 0.5)
  
  # plot(1:n_reports, spend_obj$spend, main = "alpha spending func", xlab = "look")
  
  return(bind_cols(dat, adj_alpha = spend_obj$spend))

}
# test
get_mult_compare_adj_alpha(sra_cum$data[[1]])
get_sig_tab_over_time_2(get_mult_compare_adj_alpha(sra_cum$data[[1]]))

tic()
sra_cum <-
  sra_cum %>%
  mutate(
    data = 
      map(
        .x = data, 
        .f = get_mult_compare_adj_alpha
      )
  )
toc()

# test
sra_cum$data[[10]] # check adj_alpha added as column in data

### takes ~ 100 sec
tic()
sra_cum <-
  sra_cum %>%
  mutate(
    sig_tab = 
      future_map(
        .x = data, 
        .f = get_sig_tab_over_time_2, # the alpha in data version
        .options = furrr_seed1
      )
  )
toc()



# check
sra_cum$sig_tab[[1]]


sra_cum_bcpnn_mc_adj <-
  sra_cum %>%
  unnest(cols = c(data, sig_tab)) %>%
  mutate(dte = as_date(paste0(mnth, "-01")))

sra_cum_bcpnn_mc_adj


# first signif
bcpnn_mc_adj_signif <-
  sra_cum_bcpnn_mc_adj %>%
  group_by(grps, dat_type, thresh) %>%
  dplyr::filter(ci_lo > 0) %>%
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(sra_cum_bcpnn_mc_adj)
sra_cum_bcpnn_mc_adj <-
  left_join(
    sra_cum_bcpnn_mc_adj,
    bcpnn_signif %>% select(grps, dat_type, thresh, dte_reach_sig),
    c("grps", "dat_type", "thresh")
  )
nrow(sra_cum_bcpnn_mc_adj)

sra_cum_bcpnn_mc_adj


sra_cum_bcpnn_mc_adj <- 
  sra_cum_bcpnn_mc_adj %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )





# ---- save2 ----


sra_cum_bcpnn_mc_adj %>%
  write_parquet(., sink = "out/sra_cum_bcpnn_mc_adj.parquet")







# ---- maxsprt ----


# sra_cum <- 
#   sra_dat %>%
#   dplyr::filter(dat_type == "cumulative") 
sra_cum <- 
  cumul_qtrly_dat

cv_tab <-
  sra_cum %>%
  dplyr::filter(thresh < 0.070) %>%
  group_by(grps, thresh) %>%
  summarise(
    min_dte = min(mnth),
    max_dte = max(mnth),
    rows = n(),
    sum_nA = max(nA),
    sum_nC = max(nC),
    tot_n = sum_nA + sum_nC,
    .groups = "drop"
  ) %>%
  mutate(
    qtrs = interval(paste0(min_dte, "-01"), paste0(max_dte, "-01")) / months(1) / 4,
    n_per_qtr = tot_n / qtrs,
    z = sum_nC / sum_nA
  ) 

cv_tab %>%
  kable(., digits = 1)





# testing/example
row_i <- 1
cv_tab[row_i, ]
get_maxsprt_cv(cv_tab$tot_n[row_i], floor(cv_tab$n_per_qtr[row_i]), cv_tab$z[row_i])

row_i <- 50
cv_tab[row_i, ]
get_maxsprt_cv(cv_tab$tot_n[row_i], floor(cv_tab$n_per_qtr[row_i]), cv_tab$z[row_i])




### takes ~ 70 sec
tic()
cv_tab <-
  cv_tab %>%
  # dplyr::filter(row_number() < 7) %>% ### testing
  mutate(
    cv =
      future_pmap_dbl(
        .l = list(tot_n, floor(n_per_qtr), z),
        .f = ~get_maxsprt_cv(..1, ..2, ..3),
        .options = furrr_seed3
      )
  )
toc()

cv_tab



maxsprt_dat <-
  sra_cum %>%
  mutate(
    maxllr = max_sprt_stat_(c_n = nA, n = nA + nC, z = (nC + nD) / (nA + nB)),
    rre = rr_est_(c_n = nA, n = nA + nC, z = (nC + nD) / (nA + nB))
  )

# maxsprt_dat


maxsprt_dat <-
  maxsprt_dat %>%
  inner_join(
    .,
    cv_tab %>% select(grps, thresh, cv),
    c("grps", "thresh")
  ) %>%
  mutate(reached_cv = as.integer(maxllr > cv))


maxsprt_dat %>%
  select(-dat_type) %>%
  print(., n = 25)



# ---- close_future ----


## close multisession workers by switching plan
plan(sequential)

