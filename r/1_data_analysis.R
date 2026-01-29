
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
  library("ggthemes") 
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

### NB: packages required that are used in above sourced file
# Sequential
# EmpiricalCalibration


### Note setting `plan(sequential)` for Quarto doc generation,
# errors occur otherwise
plan(sequential)

# this only applies to the non-parallel (non-"future") operations
set.seed(1234) 
# this seed can be set in future_map() etc for reproducible parallel comp seeds 
furrr_seed1 <- furrr_options(seed = 5678)
furrr_seed2 <- furrr_options(seed = 9012)
furrr_seed3 <- furrr_options(seed = 3456)
furrr_seed4 <- furrr_options(seed = 7890)




# ---- check_parallel_comp ----

# options(future.globals.maxSize = 500 * 1024 ^ 2) # = 500 MiB
options(future.globals.maxSize = 2e3 * 1024 ^ 2) # = 2 GB


# furrr parallel workers/cores setup
# change `workers = 4` based on cores available in processor being used
(thread_to_use <- parallel::detectCores() - 2) # keep a core = 2 threads free
plan(multisession, workers = thread_to_use) 

### test parallel works
# test code from https://furrr.futureverse.org/
# sequential
tic()
dev_null <- map(c(2, 2, 2), ~Sys.sleep(.x))
toc() # ~6 sec
# parallel: should be (roughly, plus overheads) a third of the time of sequential
tic()
dev_null <- future_map(c(2, 2, 2), ~Sys.sleep(.x))
toc() # ~2 sec + overhead

# for fun
tic()
dev_null <- future_map(rep(2, thread_to_use), ~Sys.sleep(.x))
toc()




# ---- consts ----

# arbitrarily, let's go with minimum cell count of 1 (will change based on context/application!)
arbitrary_cell_min <- 1


# ---- funcs ----

### Note in our notation:
# nA == (target exposure, case) count
# nB == (target exposure, control) count
# nC == (comparator exposure, case) count
# nD == (comparator exposure, control) count


# do 90% CI only with lower == one sided 0.05
get_sig_tab <- function(nA, nB, nC, nD, alpha = 0.10, method = "bcpnn", n_mcmc = 1e+05) { 
  
  out_cols_of_interest <- c("est_name", "est_scale", "est", "alpha", "ci_lo") # "ci_hi" (only care about lwr)
  sig_tab <- NULL # initialise in scope
  if (method == "bcpnn") {
    sig_tab <- pharmsignal::bcpnn_mcmc_signal(nA, nB, nC, nD, alpha = alpha, n_mcmc = n_mcmc)
  } else if (method == "prr") {
    sig_tab <- pharmsignal::prr_signal(nA, nB, nC, nD, alpha = alpha)
  } else {
    stop("method for calcaultions unknown")
  }
  sig_tab <- sig_tab[, out_cols_of_interest]
  # sig_tab <- bind_cols(tibble(mnth = mnth), sig_tab)
  return(sig_tab)
  
}

get_sig_tab_over_time <- function(dat, alpha = 0.10, method = "bcpnn", n_mcmc = 1e+05) {
  
  n_tp <- nrow(dat)
  
  sig_tab_over_time <-
    foreach(i = 1:n_tp, .combine = bind_rows, .packages = "dplyr") %do% {
      with(
        dat, 
        get_sig_tab(
          # mnth[i], 
          nA[i], nB[i], nC[i], nD[i], 
          alpha = alpha, method = method, n_mcmc = n_mcmc
        )
      )
    }
  
  return(sig_tab_over_time)
  
}



# if it's multiple comparisons central need to sparing use alpha
get_mult_compare_adj_alpha <- function(dat, alpha = 0.1) {
  
  n_reports <- nrow(dat)
  
  information_fracs <- (1:n_reports) / n_reports
  
  ### alternatives:
  # spend_obj <- sfLDPocock(alpha = alpha, t = information_fracs, param = NULL)
  # spend_obj <- sfLDOF(alpha = alpha, t = information_fracs, param = NULL)
  spend_obj <- sfExponential(alpha = alpha, t = information_fracs, param = 0.5)
  
  # plot(1:n_reports, spend_obj$spend, main = "alpha spending func", xlab = "look")
  
  return(bind_cols(dat, adj_alpha = spend_obj$spend))
  
}


# same as get_sig_tab_over_time(), however, alpha assumed included as column in data
get_sig_tab_over_time_2 <- function(dat, method = "bcpnn", n_mcmc = 1e+05) {
  
  n_tp <- nrow(dat)
  
  sig_tab_over_time <-
    foreach(i = 1:n_tp, .combine = bind_rows, .packages = "dplyr") %do% {
      with(
        dat, 
        get_sig_tab(
          # mnth[i], 
          nA[i], nB[i], nC[i], nD[i], 
          alpha = adj_alpha[i], 
          method = method,
          n_mcmc = n_mcmc
        )
      )
    }
  
  return(sig_tab_over_time)
  
}

# test 
data.frame(nA = 30, nB = 5512, nC = 41, nD = 17445, adj_alpha = 0.1) %>%
  get_sig_tab_over_time_2(.)
data.frame(nA = 30, nB = 5512, nC = 41, nD = 17445, adj_alpha = 0.1) %>%
  get_sig_tab_over_time_2(., method = "prr")
2 ^ c(0.432304, 0.7942907) # similar to prr on ratio scale
log2(c(1.556277, 2.308667)) # similar to bcpnn on log2 scale

# ---- load_dat ----


### monthly for testing
sra_dat <- read_parquet("dat/sra_dat.parquet") 

### want this
cumul_qtrly_dat <- read_parquet("dat/cumul_qtrly_dat.parquet") 

(thresholds <- sort(unique(sra_dat$thresh)))

cumul_qtrly_dat


# continuity_chk <-
#   cumul_qtrly_dat %>%
#   mutate(
#     yr = as.integer(substr(mnth, 1, 4)),
#     qtr = as.integer(substr(mnth, 7, 7))
#   )
# 
# with(
#   continuity_chk,
#   table(
#     yr,
#     qtr,
#     grps,
#     thresh,
#     useNA = "ifany"
#   )
# )
# 
# cumul_qtrly_dat %>%
#   dplyr::filter(substr(grps, 1, 3) == "(b)", thresh == "0.040")


# ---- bcpnn_calcs ----








sra_cum <-
  cumul_qtrly_dat

# make data for each combination of params nested for purrr like processing
sra_cum <-
  sra_cum %>%
  nest(data = c(mnth, nA, nB, nC, nD))

sra_cum2 <-
  sra_dat %>%
  dplyr::filter(dat_type == "cumulative") %>%
  nest(data = c(mnth, nA, nB, nC, nD))


# testing/example
sra_cum$data[[9]] %>% print(., n = nrow(.))
sra_cum2$data[[9]] %>% print(., n = nrow(.))
get_sig_tab_over_time(sra_cum$data[[9]])



### for i5-8400/48GB 2133mhz memory
# takes ~ 90 sec for monthly
# takes ~ 40 sec for quarterly
### divide by a fair bit for r9-5900X
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
sra_cum$sig_tab[[9]]


sra_cum_bcpnn <-
  sra_cum %>%
  unnest(cols = c(data, sig_tab)) %>%
  mutate(
    # dte = as_date(paste0(mnth, "-01"))
    dte = 
      as_date(paste0(
        substr(mnth, 1, 5),
        sprintf("%02.0f", (as.integer(substr(mnth, 7, 7)) - 1) * 3 + 1),
        "-01"
      ))
  )

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







# ---- multcompar_bcpnn ----

# sra_cum <- 
#   sra_dat %>%
#   dplyr::filter(dat_type == "cumulative") 
sra_cum <- 
  cumul_qtrly_dat


sra_cum <-
  sra_cum %>%
  nest(data = c(mnth, nA, nB, nC, nD))



# test get_mult_compare_adj_alpha()
get_mult_compare_adj_alpha(sra_cum$data[[11]])
get_sig_tab_over_time_2(get_mult_compare_adj_alpha(sra_cum$data[[11]]))
get_sig_tab_over_time(sra_cum$data[[11]])

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
sra_cum$data[[11]] # check adj_alpha added as column in data

### takes ~ 40 sec (i5-8400 6c/6t)
### takes ~ 55 sec on laptop (i5 8th gen 4c/8t)
### takes ~ 10 sec (R9-5900X 12c/24t)
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
sra_cum$sig_tab[[11]]


sra_cum_bcpnn_mc_adj <-
  sra_cum %>%
  unnest(cols = c(data, sig_tab)) %>%
  mutate(
    # dte = as_date(paste0(mnth, "-01"))
    dte = 
      as_date(paste0(
        substr(mnth, 1, 5),
        sprintf("%02.0f", (as.integer(substr(mnth, 7, 7)) - 1) * 3 + 1),
        "-01"
      ))
  )

sra_cum_bcpnn_mc_adj

with(sra_cum_bcpnn_mc_adj, table(dte, mnth, useNA = "ifany")) %>% 
  as.data.frame() %>%
  dplyr::filter(Freq > 0) %>%
  arrange(mnth, dte)


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
    bcpnn_mc_adj_signif %>% select(grps, dat_type, thresh, dte_reach_sig),
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





# ---- multcompar_prr ----

# sra_cum <- 
#   sra_dat %>%
#   dplyr::filter(dat_type == "cumulative") 
sra_cum <- 
  cumul_qtrly_dat


sra_cum <-
  sra_cum %>%
  nest(data = c(mnth, nA, nB, nC, nD))



# test
get_mult_compare_adj_alpha(sra_cum$data[[11]])
get_sig_tab_over_time_2(get_mult_compare_adj_alpha(sra_cum$data[[11]]))
get_sig_tab_over_time_2(get_mult_compare_adj_alpha(sra_cum$data[[11]]), method = "prr")
get_sig_tab_over_time(sra_cum$data[[11]], method = "prr")

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
sra_cum$data[[11]] # check adj_alpha added as column in data

get_sig_tab_over_time_2_prr <- function(dat) {
  get_sig_tab_over_time_2(dat, method = "prr")
}


### takes ~2 sec on laptop (i5 8th gen 4c/8t)
tic()
sra_cum <-
  sra_cum %>%
  mutate(
    sig_tab = 
      future_map(
        .x = data, 
        .f = get_sig_tab_over_time_2_prr, # the alpha in data version
        .options = furrr_seed1
      )
  )
toc()



# check
sra_cum$sig_tab[[11]]


sra_cum_prr_mc_adj <-
  sra_cum %>%
  unnest(cols = c(data, sig_tab)) %>%
  mutate(
    # dte = as_date(paste0(mnth, "-01"))
    dte = 
      as_date(paste0(
        substr(mnth, 1, 5),
        sprintf("%02.0f", (as.integer(substr(mnth, 7, 7)) - 1) * 3 + 1),
        "-01"
      ))
  )

sra_cum_prr_mc_adj

with(sra_cum_prr_mc_adj, table(dte, mnth, useNA = "ifany")) %>% 
  as.data.frame() %>%
  dplyr::filter(Freq > 0) %>%
  arrange(mnth, dte)


# first signif
prr_mc_adj_signif <-
  sra_cum_prr_mc_adj %>%
  group_by(grps, dat_type, thresh) %>%
  dplyr::filter(ci_lo > 1) %>% # 1 is the critical value on ratio scale
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(sra_cum_prr_mc_adj)
sra_cum_prr_mc_adj <-
  left_join(
    sra_cum_prr_mc_adj,
    prr_mc_adj_signif %>% select(grps, dat_type, thresh, dte_reach_sig),
    c("grps", "dat_type", "thresh")
  )
nrow(sra_cum_prr_mc_adj)

sra_cum_prr_mc_adj


sra_cum_prr_mc_adj %>%
  arrange(grps, thresh, dte, mnth) %>%
  group_by(grps, thresh, dte, mnth) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  dplyr::filter(n > 1)

sra_cum_prr_mc_adj %>%
  dplyr::filter(thresh == "0.070", grepl("(a)", grps, fixed = TRUE))


sra_cum_prr_mc_adj %>%
  dplyr::filter(thresh == "0.050", grepl("(c)", grps, fixed = TRUE)) %>%
  print(., n = nrow(.))

sra_cum_prr_mc_adj <- 
  sra_cum_prr_mc_adj %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )





# ---- save3 ----


sra_cum_prr_mc_adj %>%
  write_parquet(., sink = "out/sra_cum_prr_mc_adj.parquet")





# ---- maxsprt ----


# sra_cum <- 
#   sra_dat %>%
#   dplyr::filter(dat_type == "cumulative") 
sra_cum <- 
  cumul_qtrly_dat

cv_tab <-
  sra_cum %>%
  # dplyr::filter(thresh < 0.070) %>%
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
    # qtrs = interval(paste0(min_dte, "-01"), paste0(max_dte, "-01")) / months(1) / 4,
    qtrs = rows,
    n_per_qtr = tot_n / qtrs,
    z = sum_nC / sum_nA
  ) 

cv_tab %>%
  kable(., digits = 1)


# maxsprt: create alternative CV tab 



### create CV tab for alternative n_per_qtr and z ratios
alt_mults <-
  tribble(
    ~alt_str, ~modifier, ~mult,
    "quar_n", "n_per_qtr", 0.25,
    "half_n", "n_per_qtr", 0.5,
    "doub_n", "n_per_qtr", 2  ,
    "quad_n", "n_per_qtr", 4  ,
    "quar_z",         "z", 0.25,
    "half_z",         "z", 0.5,
    "doub_z",         "z", 2  ,
    "quad_z",         "z", 4  
  )


cv_tab_alts <-
  cross_join(
    alt_mults,
    cv_tab
  ) %>%
  arrange(
    grps, thresh, modifier, mult, alt_str
  )

if (nrow(alt_mults) * nrow(cv_tab) != nrow(cv_tab_alts)) {
  stop("cross_join() has gone wrong")
}



# maxsprt: create CVs 

# testing/example
row_i <- 1
cv_tab[row_i, ]
get_maxsprt_cv(cv_tab$tot_n[row_i], floor(cv_tab$n_per_qtr[row_i]), cv_tab$z[row_i])

row_i <- 50
cv_tab[row_i, ]
get_maxsprt_cv(cv_tab$tot_n[row_i], floor(cv_tab$n_per_qtr[row_i]), cv_tab$z[row_i])



### takes ~ 70 sec (i5-8400)
# note purrr::possibly() will just catch when model fails and return as.numeric(NA) 
get_maxsprt_cv_poss <- 
  possibly(get_maxsprt_cv, otherwise = NA_real_, quiet = FALSE)

tic()
cv_tab <-
  cv_tab %>%
  # dplyr::filter(row_number() < 7) %>% ### testing
  mutate(
    cv =
      future_pmap_dbl(
        .l = list(tot_n, floor(n_per_qtr), z),
        .f = ~get_maxsprt_cv_poss(..1, ..2, ..3),
        .options = furrr_seed3
      )
  )
toc()

cv_tab
cv_tab %>% dplyr::filter(is.na(cv))
# remove analyses where thresholds don't allow enough events (extreme threshold values)
# cv_tab <- cv_tab %>% dplyr::filter(!is.na(cv))



# maxsprt: create alt CVs 

cv_tab_alts <-
  cv_tab_alts %>%
  mutate(
    n_per_qtr = if_else(modifier == "n_per_qtr", mult * n_per_qtr, n_per_qtr),
    z         = if_else(modifier ==         "z", mult * z        , z        ),
  )

cv_tab_alts


### takes ~ 270 sec (R9-5900X)
tic()
cv_tab_alts <-
  cv_tab_alts %>%
  # dplyr::filter(row_number() < 7) %>% ### testing
  mutate(
    cv =
      future_pmap_dbl(
        .l = list(tot_n, floor(n_per_qtr), z),
        .f = ~get_maxsprt_cv_poss(..1, ..2, ..3),
        .options = furrr_seed4
      )
  )
toc()

cv_tab_alts
cv_tab_alts %>% dplyr::filter(is.na(cv))


# include original CVs too
cv_tab_alts <-
  bind_rows(
    cv_tab_alts,
    cv_tab %>% mutate(alt_str = "same_n", modifier = "n_per_qtr", mult = 1),
    cv_tab %>% mutate(alt_str = "same_z", modifier = "z", mult = 1)
  ) %>%
  arrange(grps, thresh, modifier, mult, alt_str)




# maxsprt: create llr test stats 


maxsprt_dat_calcs <-
  sra_cum %>%
  mutate(
    maxllr = max_sprt_stat_(c_n = nA, n = nA + nC, z = (nC + nD) / (nA + nB)),
    rre = rr_est_(c_n = nA, n = nA + nC, z = (nC + nD) / (nA + nB))
  )

# maxsprt_dat
# maxsprt_dat %>% dplyr::filter(thresh == "0.100", substr(grps, 1, 3) == "(a)")

maxsprt_dat <-
  maxsprt_dat_calcs %>%
  left_join(
    .,
    cv_tab %>% select(grps, thresh, cv),
    c("grps", "thresh")
  ) 

maxsprt_dat <-
  maxsprt_dat %>%
  mutate(
    # some cvs don't exist so those llr never reach cv
    reached_cv = if_else(is.na(cv), 0L, as.integer(maxllr > cv)),
    # create date for start of each quarter
    dte = 
      as_date(paste0(
        substr(mnth, 1, 5),
        sprintf("%02.0f", (as.integer(substr(mnth, 7, 7)) - 1) * 3 + 1),
        "-01"
      ))
  )

maxsprt_dat %>% dplyr::filter(is.na(cv))

# have a peak
maxsprt_dat %>%
  select(-dat_type) %>%
  print(., n = 25)


# first signif
maxsprt_signif <-
  maxsprt_dat %>%
  group_by(grps, dat_type, thresh) %>%
  dplyr::filter(reached_cv > 0) %>%
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(maxsprt_dat)
maxsprt_dat <-
  left_join(
    maxsprt_dat,
    maxsprt_signif %>% select(grps, dat_type, thresh, dte_reach_sig),
    c("grps", "dat_type", "thresh")
  )
nrow(maxsprt_dat)

maxsprt_dat


maxsprt_dat <- 
  maxsprt_dat %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )

# these are where the maxllr has dropped under the CV after exceeding it previously
maxsprt_dat %>%
  dplyr::filter(
    is.na(reach_sig) | 
      is.na(reached_cv) | 
      (as.logical(reached_cv) != reach_sig)
  )



maxsprt_dat <- 
  maxsprt_dat %>%
  select(-reached_cv)




# maxsprt: create llr test stats for alt CVs


nrow(maxsprt_dat_calcs)
maxsprt_dat_alts <-
  maxsprt_dat_calcs %>%
  left_join(
    .,
    cv_tab_alts %>% select(alt_str, modifier, mult, grps, thresh, cv),
    c("grps", "thresh"),
    relationship = "many-to-many"
  ) %>%
  arrange(grps, thresh, modifier, mult, alt_str, mnth) %>%
  select(grps, thresh, modifier, mult, alt_str, everything())
nrow(maxsprt_dat_alts)


if(nrow(maxsprt_dat_alts) != (nrow(alt_mults) + 2) * nrow(maxsprt_dat_calcs)) {
  stop("many-to-many join has not worked")
}

maxsprt_dat_alts


maxsprt_dat_alts <-
  maxsprt_dat_alts %>%
  mutate(
    # some cvs don't exist so those llr never reach cv
    reached_cv = if_else(is.na(cv), 0L, as.integer(maxllr > cv)),
    # create date for start of each quarter
    dte = 
      as_date(paste0(
        substr(mnth, 1, 5),
        sprintf("%02.0f", (as.integer(substr(mnth, 7, 7)) - 1) * 3 + 1),
        "-01"
      ))
  )

maxsprt_dat %>% dplyr::filter(is.na(cv))

# have a peak
maxsprt_dat_alts %>%
  select(-dat_type) %>%
  print(., n = 25)


# first signif
maxsprt_alts_signif <-
  maxsprt_dat_alts %>%
  group_by(grps, dat_type, thresh, modifier, mult, alt_str) %>%
  dplyr::filter(reached_cv > 0) %>%
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(maxsprt_dat_alts)
maxsprt_dat_alts <-
  left_join(
    maxsprt_dat_alts,
    maxsprt_alts_signif %>% 
      select(grps, dat_type, thresh, modifier, mult, alt_str, dte_reach_sig),
    c("grps", "dat_type", "thresh", "modifier", "mult", "alt_str")
  )
nrow(maxsprt_dat_alts)

maxsprt_dat_alts


maxsprt_dat_alts <- 
  maxsprt_dat_alts %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )

# these are where the maxllr has dropped under the CV after exceeding it previously
maxsprt_dat_alts %>%
  dplyr::filter(
    is.na(reach_sig) | 
      is.na(reached_cv) | 
      (as.logical(reached_cv) != reach_sig)
  )



maxsprt_dat_alts <- 
  maxsprt_dat_alts %>%
  select(-reached_cv)





# ---- save4 ----


maxsprt_dat %>%
  write_parquet(., sink = "out/sra_cum_maxsprt.parquet")




maxsprt_dat_alts %>%
  write_parquet(., sink = "out/sra_cum_maxsprt_alt_cvs.parquet")







# ---- close_future ----


## close multisession workers by switching plan
plan(sequential)

