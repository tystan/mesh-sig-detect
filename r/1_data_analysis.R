
# ---- lib ----


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


# ---- consts ----

# arbitrarily, let's go with minimum cell count of 3 (should be discussed!)
arbitrary_cell_min <- 1



# ---- funcs ----


get_sig_tab <- function(nA, nB, nC, nD, alpha = 0.05, n_mcmc = 1e+05) {
  
  out_cols_of_interest <- c("est_name", "est_scale", "est", "ci_lo", "ci_hi")
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


# ---- load_dat ----



sra_dat <- read_parquet("dat/sra_dat.parquet")


# ---- bcpnn_calcs ----






sra_cum <- 
  sra_dat %>%
  dplyr::filter(dat_type == "cumulative") 

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



# ---- end ----

## close multisession workers by switching plan
plan(sequential)



sra_cum_bcpnn %>%
  write_parquet(., sink = "out/sra_cum_bcpnn.parquet")







# ---- maxsprt ----

sra_snp <- 
  sra_dat %>%
  dplyr::filter(dat_type == "snapshot")








# ---- multcomp ----

dat11 <-
  foreach(i = 1:length(thresholds), .combine = bind_rows, .packages = "dplyr") %do% {
    
    dat1 <-
      get_signal_dat(
        g1 = "pelvic_mesh",
        g2 = "hernia_mesh",
        pain_type = "pain_topic", 
        thresh = thresholds[i],
        cell_min = 1,
        verbose = FALSE
      ) %>%
      bind_cols(., thresh = thresholds[i])
    
    
    
    
    # so this is multiple comparisons central but let's create disproportionality stats
    # whenever a new report enters the data
    n_reports <- nrow(dat1)
    
    information_fracs <-  1:n_reports / n_reports
    # spend_obj <- sfLDPocock(alpha = 0.025, t = information_fracs, param = NULL)
    # spend_obj <- sfLDOF(alpha = 0.025, t = information_fracs, param = NULL)
    spend_obj <- sfExponential(alpha = 0.05, t = information_fracs, param = 0.5)
    
    # plot(1:n_reports, spend_obj$spend, main = "alpha spending func", xlab = "look")
    
    
    # takes ~1 sec using i5-8400
    system.time({
      da_stats <-
        foreach(i = 1:n_reports, .combine = bind_rows, .packages = "dplyr") %dopar% {
          with(
            dat1, 
            pretty_da(
              mnth[i], thresh[i], nA[i], nB[i], nC[i], nD[i],
              alpha = spend_obj$spend[i]
            )
          ) %>%
            mutate(alpha = spend_obj$spend[i])
        }
    })
    
    dat1 <-
      dat1 %>%
      inner_join(
        .,
        da_stats,
        c("mnth", "thresh")
      )
    
    n_post_join <- nrow(dat1)
    if (n_post_join != n_reports) {
      stop("join not 1-1")
    }
    
    dat1
  }




# first signif
bcpnn_mult_comp_signif <-
  dat11 %>%
  group_by(thresh) %>%
  dplyr::filter(ci_lo > 0) %>%
  arrange(mnth) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(
    dte = as_date(paste0(mnth, "-01")),
    est_name = paste0(est_name, "(MCadj)")
  )





