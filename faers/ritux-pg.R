
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
  sig_dat <-
    read_parquet(
      file = "faers/ritux_pg_signal_data.parquet"
    )
})

with(sig_dat, table(exposure, pt, comparator))


sig_dat %>% 
  dplyr::filter(comparator == "(1) All")

sig_dat %>% 
  sample_n(10)

# check data is in tibble format
is_tibble(sig_dat)
class(sig_dat)
nrow(sig_dat)


# ---- create_all_comp_record_data ----

sig_dat %>% 
  dplyr::filter(comparator == "(1) All") %>% 
  summarise(n = sum(record_cnt))


system.time({
  # fda logged dates
  fda_dts <-
    read_parquet(
      file = "faers/faers_demo_fda-qtr.parquet"
    )
})

fda_dts <- as_tibble(fda_dts)

# check data is in tibble format
is_tibble(fda_dts)
class(fda_dts)
nrow(fda_dts)
colnames(fda_dts)

# fda_dts <-
#  fda_dts %>%
#  mutate(fda_dt_lub = as_date(as.character(fda_dt)))

with(
  fda_dts, 
  table(substr(as.character(qtr), 1, 4), useNA = "ifany")
)


sig_dat <-
  sig_dat %>% 
  select(primaryid, comparator, exposure, pt) %>%
  dplyr::filter(comparator != "(1) All")


nrow(sig_dat)
sig_dat <-
  left_join(
    sig_dat,
    fda_dts %>% select(primaryid, qtr),
    "primaryid"
  ) # %>% dplyr::filter(is.na(qtr)) 

nrow(sig_dat)

sig_dat %>% dplyr::filter(is.na(qtr)) 




# da_dts %>%
#   dplyr::filter(substr(as.character(fda_dt), 1, 4) >= "2013") %>%
#   summarise(n = n())


rm_ids <- read_csv("faers/ritux_pg_rm_ids.csv")

fda_dts <-
  fda_dts %>%
  anti_join(
    .,
    rm_ids,
    "primaryid"
  )


pg_ids <- read_csv("faers/all_pg.csv")

fda_dts <-
  fda_dts %>%
  left_join(
    .,
    pg_ids %>% select(primaryid, pt),
    "primaryid"
  ) %>%
  mutate(pt = ifelse(is.na(pt), "not [PG]", "[PG]"))

fda_dts


fda_dts <-
  fda_dts %>%
  left_join(
    .,
    sig_dat %>% 
      dplyr::filter(comparator == "(2) mAbs", exposure == "[rituximab]") %>% 
      select(primaryid, exposure),
    "primaryid"
  ) %>%
  mutate(
    exposure = ifelse(is.na(exposure), "not [rituximab]", exposure),
    comparator = "(1) All"
  )

fda_dts






# ---- bind_all_comp_and_rest_comp ----


set_mem <-
  bind_rows(
    fda_dts,
    sig_dat
  )

set_mem %>%
  group_by(comparator, exposure, pt) %>%
  summarise(n = n())



# ---- create_counts ----



signal_data <-
  set_mem %>%
  select(comparator, qtr, exposure, outcome = pt) %>%
  group_by(comparator, qtr, exposure, outcome) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(comparator, qtr, exposure, outcome) 



# free up memory of patient level data
rm(list = c("set_mem", "sig_dat", "fda_dts", "rm_ids", "pg_ids"))
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
colnames(signal_data_wide) <- gsub("\\[rituximab\\]", "[E]", colnames(signal_data_wide))
colnames(signal_data_wide) <- gsub("\\[PG\\]", "[O]", colnames(signal_data_wide))
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
  #  dplyr::filter(comparator == "(3) IRAs") %>%
  arrange(comparator, qtr)


signal_data_wide <-
  signal_data_wide %>%
  group_by(comparator) %>%
  mutate(
    a = cumsum(a),
    b = cumsum(b),
    c = cumsum(c),
    d = cumsum(d)
  ) %>%
  ungroup()

signal_data_wide %>%
  dplyr::filter(comparator == "(3) CD20s") 

signal_data_wide %>%
  dplyr::filter(substr(comparator, 1, 3) == "(2)") %>%
  print(., n = nrow(.))

signal_data_wide %>%
  dplyr::filter(substr(comparator, 1, 3) == "(1)") 

# ---- maxsprt ----



cv_tab <-
  signal_data_wide %>%
  group_by(comparator) %>%
  summarise(
    rows = n(),
    sum_A = max(a),
    sum_C = max(c),
    tot_n = sum_A + sum_C,
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    # qtrs = interval(paste0(min_dte, "-01"), paste0(max_dte, "-01")) / months(1) / 4,
    qtrs = rows,
    n_per_qtr = floor(tot_n / qtrs),
    n_per_qtr = ifelse(n_per_qtr < 1, 1, n_per_qtr),
    z = sum_C / sum_A
  ) 

cv_tab %>%
  kable(., digits = 1)

# note purrr::possibly() will just catch when model fails and return as.numeric(NA) 
get_maxsprt_cv_poss <- 
  possibly(get_maxsprt_cv, otherwise = NA_real_, quiet = FALSE)


cv_tab_test <- cv_tab %>% dplyr::filter(comparator == "(3) CD20s") 
cv <- with(cv_tab_test, get_maxsprt_cv_poss(tot_n, floor(n_per_qtr), z))
cv


cv_tab <-
  cv_tab %>%
  mutate(
    cv = pmap_dbl(.l = list(tot_n, n_per_qtr, z), .f = get_maxsprt_cv_poss)
  )

cv_tab %>% dplyr::filter(comparator == "(3) CD20s") %>% pull(cv)

maxsprt_dat <-
  signal_data_wide %>%
  mutate(
    maxllr = max_sprt_stat_(c_n = a, n = a + c, z = (c + d) / (a + b)),
    rre = rr_est_(c_n = a, n = a + c, z = (c + d) / (a + b))
  )

for (i in 1:nrow(signal_data_wide)) {
  print(i)

signal_data_wide %>%
  dplyr::filter(row_number() == i) %>%
  mutate(
    maxllr = max_sprt_stat(c_n = a, n = a + c, z = (c + d) / (a + b)),
    rre = rr_est(c_n = a, n = a + c, z = (c + d) / (a + b))
  )
}

signal_data_wide %>%
  dplyr::filter(row_number() == 126) %>%
  mutate(
    maxllr = max_sprt_stat(c_n = a, n = a + c, z = (c + d) / (a + b)),
    rre = rr_est(c_n = a, n = a + c, z = (c + d) / (a + b))
  )




maxsprt_dat <-
  maxsprt_dat %>%
  inner_join(
    .,
    cv_tab %>% select(comparator, cv),
    "comparator"
  )


maxsprt_dat



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



# ---- mult_adj_bcpnn ----



# if it's multiple comparisons central need to sparing use alpha
get_mult_compare_adj_alpha <- function(vec, alpha = 0.1) {
  
  n_reports <- length(vec)
  
  information_fracs <-  1:n_reports / n_reports
  
  ### alternatives:
  # spend_obj <- sfLDPocock(alpha = 0.025, t = information_fracs, param = NULL)
  # spend_obj <- sfLDOF(alpha = 0.025, t = information_fracs, param = NULL)
  spend_obj <- sfExponential(alpha = alpha, t = information_fracs, param = 0.5)
  
  # plot(1:n_reports, spend_obj$spend, main = "alpha spending func", xlab = "look")
  
  # return(bind_cols(dat, adj_alpha = spend_obj$spend))
  return(spend_obj$spend)
  
}
# test
maxsprt_dat %>% 
  dplyr::filter(comparator == "(3) IRAs") %>%
  bind_cols(., adj_alpha = get_mult_compare_adj_alpha(.[["comparator"]]))


get_sig_tab <- function(nA, nB, nC, nD, alpha = 0.1, n_mcmc = 1e+05) {
  
  out_cols_of_interest <- c("est_name", "est_scale", "est", "alpha", "ci_lo")
  sig_tab <- pharmsignal::bcpnn_mcmc_signal(nA, nB, nC, nD, alpha = alpha, n_mcmc = n_mcmc)
  sig_tab <- sig_tab[, out_cols_of_interest]
  return(sig_tab)
  
}

get_sig_tab(30,  5512, 41, 17445)


get_sig_tab_over_time <- function(dat, alpha = 0.1, n_mcmc = 1e+05) {
  
  n_tp <- nrow(dat)
  # dat <- get_mult_compare_adj_alpha(dat)
  
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
  maxsprt_dat %>% 
  group_by(comparator) %>%
  mutate(adj_alpha = get_mult_compare_adj_alpha(comparator)) %>%
  ungroup()


bcpnn_data %>% dplyr::filter(comparator == "(3) IRAs")
bcpnn_data %>% dplyr::filter(comparator == "(4b) CU indi")


bcpnn_data <-
  bind_cols(
    bcpnn_data %>% select(comparator, dte),
    get_sig_tab_over_time(bcpnn_data)
  )

# truncate exceedingly negative values for plotting 
bcpnn_data <-
  bcpnn_data %>%
  mutate(
    ci_lo = if_else(ci_lo < -5, -5, ci_lo)
  )


bcpnn_data %>%
  kable(.)



# ---- plot_results ----



plt_dat <-
  bind_rows(
    maxsprt_dat %>% 
      select(comparator, dte, cv, reached_cv, val = maxllr) %>%
      mutate(stat = "MaxSPRT (max LLR)"),
    bcpnn_data %>% 
      select(comparator, dte, val = ci_lo) %>%
      mutate(cv = 0, reached_cv = as.integer(val > cv), stat = "IC (BCPNN, Lower 95% CI )")
  ) 


sig_reach_dat <-
  plt_dat %>%
  arrange(stat, comparator, dte) %>%
  group_by(stat, comparator) %>%
  dplyr::filter(reached_cv == 1) %>%
  dplyr::filter(row_number() == 1) %>%
  select(stat, comparator, dte_reached = dte) %>%
  # now create separation between reached CV values when it occurs
  group_by(stat, dte_reached) %>%
  mutate(rep_dte = 1:n()) %>%
  ungroup() %>%
  mutate(dte_reached = dte_reached + days(10 * (rep_dte - 1))) %>%
  select(-rep_dte)


plt_dat <-
  left_join(
    plt_dat,
    sig_reach_dat,
    c("stat", "comparator")
  )

plt_dat %>%
  ggplot(., aes(x = dte, y = val, col = comparator )) +
  geom_hline(aes(yintercept = cv), alpha = 0.5) +
  geom_vline(aes(xintercept = dte_reached, col = comparator), alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_point() +
  facet_wrap(~ stat, ncol = 1, scales = "free_y") +
  scale_colour_tableau() +
  theme_bw() +
  labs(
    x = "Quarter",
    y = "Statistic",
    col = "Comparison"
  )


ggsave(
  filename = "faers/vedol_panc_signal_detection_over_time.pdf",
  width = 10,
  height = 8
)

ggsave(
  filename = "faers/vedol_panc_signal_detection_over_time.png",
  width = 10,
  height = 8
)


