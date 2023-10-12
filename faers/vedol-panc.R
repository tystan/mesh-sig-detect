
# ---- libs ----

library("arrow") # reading in parquet data
library("dplyr") # piping and manipulation of data
library("tibble") # prettier data.frames
library("ggplot2") # plotting
library("ggthemes") # plotting
library("knitr") # to print kable(.) tables
library("tidyr") # pivot_wider() function
library("readr") # 
library("lubridate") #
library("purrr")
library("foreach")
library("gsDesign")


# NOTE: need to run first (only once, assumes devtools installed):
# devtools::install_github("tystan/pharmsignal") 
library("pharmsignal") # signal detection algs/disproportionality statistics


# here are the functions written for these analyses
# they will be shown in the *Appendix A*
source("r/_funcs.R")

### NB: packages required that are used in above sourced file
# Sequential
# EmpiricalCalibration



# ---- load ----

# load set_membership data set that has set identifiers for each id

# read in line data using arrow::read_parquet()
# (data saved in parquet format for speed/storage size)
system.time({
  # set_membership
  set_mem <-
    read_parquet(
      file = "faers/set_mem.parquet"
    )
})

# check data is in tibble format
is_tibble(set_mem)
class(set_mem)
nrow(set_mem)


# reports before 2015 and having [panc] indication are not required
set_mem <-
  set_mem %>% 
  dplyr::filter(ind == "not [panc]", date == ">=2015")

nrow(set_mem)


system.time({
  # fda logged dates
  fda_dts <-
    read_parquet(
      file = "faers/faers_demo_fda-date.parquet"
    )
})

# check data is in tibble format
is_tibble(fda_dts)
class(fda_dts)
nrow(fda_dts)


fda_dts <-
  fda_dts %>%
  mutate(fda_dt_lub = as_date(as.character(fda_dt)))

with(fda_dts, table(is.na(fda_dt_lub), useNA = "ifany"))
with(fda_dts, max(fda_dt_lub))
with(fda_dts, min(fda_dt_lub))
with(
  fda_dts, 
  table(substr(as.character(fda_dt), 1, 4), is.na(fda_dt_lub), useNA = "ifany")
)

fda_dts <-
  fda_dts %>%
  select(-fda_dt) %>%
  rename(fda_dt = fda_dt_lub)

fda_dts

set_mem


anti_join(
  set_mem,
  fda_dts,
  "primaryid"
)

nrow(set_mem)
set_mem <-
  inner_join(
    set_mem,
    fda_dts,
    "primaryid"
  )
nrow(set_mem)

with(
  set_mem, 
  table(year(fda_dt), date, useNA = "ifany")
)

set_mem <-
  set_mem %>%
  mutate(fda_yr = year(fda_dt), fda_qtr = quarter(fda_dt))

with(
  set_mem, 
  table(fda_yr, fda_qtr, useNA = "ifany")
)



# standardise the summary data once the dataset has been filtered as required
mk_signal_data <- function(dat, comparator_str) {
  
  dat %>%
    mutate(
      comparator = comparator_str, 
      med = if_else(med == "[vedol]", med, "not [vedol]")
    ) %>%
    select(comparator, exposure = med, outcome = reac, fda_yr, fda_qtr) %>%
    group_by(comparator, exposure, outcome, fda_yr, fda_qtr) %>%
    summarise(n = n(), .groups = "keep") %>%
    ungroup() %>%
    arrange(fda_yr, fda_qtr, comparator, exposure, outcome) 
  
}


signal_data <-
  set_mem %>%
  mk_signal_data(., "(1) All")



signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(med %in% c("[vedol]", "[IRA]", "[mAb]")) %>%
      mk_signal_data(., "(2) mAbs")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(med %in% c("[vedol]", "[IRA]")) %>%
      mk_signal_data(., "(3) IRAs")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(ind_ibd == "[IBD]") %>%
      mk_signal_data(., "(4) All IBD indi")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(ind_cd == "[CD]") %>%
      mk_signal_data(., "(4a) CD indi")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(ind_cu == "[CU]") %>%
      mk_signal_data(., "(4b) CU indi")
  )

# free up memory of patient level data
rm(list = "set_mem")
gc()

signal_data %>% sample_n(20)



# create a,b,c,d cell counts as columns  
signal_data_wide <-
  signal_data %>%
  pivot_wider(
    names_from = c("exposure", "outcome"), 
    values_from = "n",
    names_sep = "_",
    values_fill = 0
  )

# have a look
signal_data_wide



# edit exposure and outcome values for shorter expressions
# they now take form: "ExOy" = `Exposure <x> Outcome <y>` where
# <x> = "p" (positive) if exposure == "vedol", "n" (negative) otherwise
# <y> = "p" (positive) if outcome == "panc", "n" (negative) otherwise
colnames(signal_data_wide) <- gsub("\\[vedol\\]", "[E]", colnames(signal_data_wide))
colnames(signal_data_wide) <- gsub("\\[panc\\]", "[O]", colnames(signal_data_wide))
colnames(signal_data_wide) <- gsub("not \\[([EO])\\]", "[\\1]n", colnames(signal_data_wide))
colnames(signal_data_wide) <- gsub("\\](_|$)", "]p\\1", colnames(signal_data_wide))
# rm punctuation greedily
colnames(signal_data_wide) <- gsub("(\\[|\\]|_)", "", colnames(signal_data_wide)) 

# have a look
signal_data_wide



signal_data_wide <-
  signal_data_wide %>%
    rename(
      a = EpOp,
      b = EpOn,
      c = EnOp,
      d = EnOn
    )

with(signal_data_wide, table(comparator))

signal_data_wide <-
  signal_data_wide %>%
  dplyr::filter(comparator == "(3) IRAs") %>%
  arrange(comparator, fdayr, fdaqtr)


signal_data_wide <-
  signal_data_wide %>%
  mutate(
    a = cumsum(a),
    b = cumsum(b),
    c = cumsum(c),
    d = cumsum(d)
  )

signal_data_wide


cv_tab <-
  signal_data_wide %>%
  # dplyr::filter(thresh < 0.070) %>%
  summarise(
    rows = n(),
    sum_A = max(a),
    sum_C = max(c),
    tot_n = sum_A + sum_C,
    .groups = "drop"
  ) %>%
  mutate(
    # qtrs = interval(paste0(min_dte, "-01"), paste0(max_dte, "-01")) / months(1) / 4,
    qtrs = rows,
    n_per_qtr = tot_n / qtrs,
    z = sum_C / sum_A
  ) 

cv_tab %>%
  kable(., digits = 1)

# note purrr::possibly() will just catch when model fails and return as.numeric(NA) 
get_maxsprt_cv_poss <- 
  possibly(get_maxsprt_cv, otherwise = NA_real_, quiet = FALSE)



cv <- with(cv_tab, get_maxsprt_cv_poss(tot_n, floor(n_per_qtr), z))
cv



maxsprt_dat <-
  signal_data_wide %>%
  mutate(
    maxllr = max_sprt_stat_(c_n = a, n = a + c, z = (c + d) / (a + b)),
    rre = rr_est_(c_n = a, n = a + c, z = (c + d) / (a + b)),
    cv = cv
  )

# maxsprt_dat



maxsprt_dat <-
  maxsprt_dat %>%
  mutate(
    # some cvs don't exist so those llr never reach cv
    reached_cv = if_else(is.na(cv), 0L, as.integer(maxllr > cv)),
    # create date for start of each quarter
    dte = 
      as_date(paste0(
        fdayr,
        "-",
        sprintf("%02.0f", (fdaqtr - 1) * 3 + 1),
        "-01"
      ))
  ) 

maxsprt_dat %>%
  kable(.)



# if it's multiple comparisons central need to sparing use alpha
get_mult_compare_adj_alpha <- function(dat, alpha = 0.1) {
  
  n_reports <- nrow(dat)
  
  information_fracs <-  1:n_reports / n_reports
  
  ### alternatives:
  # spend_obj <- sfLDPocock(alpha = 0.025, t = information_fracs, param = NULL)
  # spend_obj <- sfLDOF(alpha = 0.025, t = information_fracs, param = NULL)
  spend_obj <- sfExponential(alpha = alpha, t = information_fracs, param = 0.5)
  
  # plot(1:n_reports, spend_obj$spend, main = "alpha spending func", xlab = "look")
  
  return(bind_cols(dat, adj_alpha = spend_obj$spend))
  
}
# test
get_mult_compare_adj_alpha(maxsprt_dat)


get_sig_tab <- function(nA, nB, nC, nD, alpha = 0.1, n_mcmc = 1e+05) {
  
  out_cols_of_interest <- c("est_name", "est_scale", "est", "alpha", "ci_lo", "ci_hi")
  sig_tab <- pharmsignal::bcpnn_mcmc_signal(nA, nB, nC, nD, alpha = alpha, n_mcmc = n_mcmc)
  sig_tab <- sig_tab[, out_cols_of_interest]
  return(sig_tab)
  
}

get_sig_tab(30,  5512, 41, 17445)


get_sig_tab_over_time <- function(dat, alpha = 0.1, n_mcmc = 1e+05) {
  
  n_tp <- nrow(dat)
  dat <- get_mult_compare_adj_alpha(dat)
  
  sig_tab_over_time <-
    foreach(i = 1:n_tp, .combine = bind_rows, .packages = "dplyr") %do% {
      with(
        dat, 
        get_sig_tab(
          # mnth[i], 
          a[i], b[i], c[i], d[i], 
          alpha = adj_alpha[i], n_mcmc = n_mcmc
        )
      )
    }
  
  return(sig_tab_over_time)
  
}

bcpnn_data <-
  bind_cols(
    maxsprt_dat %>% select(comparator, dte),
    get_sig_tab_over_time(maxsprt_dat)
  )

bcpnn_data %>%
  kable(.)

plt_dat <-
  bind_rows(
    maxsprt_dat %>% 
      select(comparator, dte, cv, reached_cv, val = maxllr) %>%
      mutate(stat = "MaxSPRT (max LLR)"),
    bcpnn_data %>% 
      select(comparator, dte, val = ci_lo) %>%
      mutate(cv = 0, reached_cv = as.integer(val > cv), stat = "IC (BCPNN, Lower 95% CI )")
  ) 

plt_dat <-
  plt_dat %>%
  arrange(stat, comparator, dte) %>%
  group_by(stat, comparator) %>%
  dplyr::filter(reached_cv == 1) %>%
  dplyr::filter(row_number() == 1) %>%
  select(stat, comparator, dte_reached = dte) %>%
  left_join(
    plt_dat,
    .,
    c("stat", "comparator")
  )

plt_dat %>%
  ggplot(., aes(x = dte, y = val, col = comparator )) +
  geom_hline(aes(yintercept = cv)) +
  geom_line(alpha = 0.5) +
  geom_point() +
  geom_vline(aes(xintercept = dte_reached, col = comparator), alpha = 0.5) +
  facet_wrap(~ stat, ncol = 1, scales = "free_y") +
  scale_colour_tableau() +
  theme_bw() +
  labs(
    x = "Quarter",
    y = "Statistic",
    col = "Comparison"
  )



