


# ---- libs ----


library("arrow")     # parquet files
library("binom")     # wilson confidence intervals for binomial counts
library("foreach")   # flexible looping and return amalgamation
library("episensr")  # misclassification error contingency table adjustments
# library("simdata") # NORTA method to get correlation estimates 
library("tictoc")    # time how long things take 
library("pharmsignal") # remotes::install_github("tystan/pharmsignal")
library("gsDesign")
# citation("episensr")

# tidyverse et al
library("readr")
library("dplyr")
library("purrr")
library("tidyr")
library("forcats")
library("testthat")
library("stringr")
library("lubridate")
library("furrr")
library("ggplot2")
library("ggthemes")




# ---- func ----

### get_sens_spec():
# create table of ("sens", "spec", "ppv", "npv") from data.frame containing
# observation/individual level "pain_truth" and "pain_test" columns
get_sens_spec <- function(df, pos = "Yes", alpha = 0.05) {
  sens_rws <- df[["pain_truth"]] %in% as.character(pos)
  ppv_rws <- df[["pain_test"]] %in% as.character(pos)
  tp_rws <- sens_rws & ppv_rws
  spec_rws <- !sens_rws
  npv_rws <- !ppv_rws
  tn_rws <- spec_rws & npv_rws
  cases <- c(sum(sens_rws), sum(spec_rws), sum(ppv_rws), sum(npv_rws))
  correct <- c(sum(tp_rws), sum(tn_rws), sum(tp_rws), sum(tn_rws)) 
  conf_ints <- 
    binom.confint(
      correct, 
      cases, 
      conf.level = 1 - alpha, 
      methods = "wilson"
    )
  return(tibble(
    param = c("sens", "spec", "ppv", "npv"), 
    cases = cases, 
    correct = correct, 
    est = correct / cases, 
    lo = conf_ints$lower, 
    up = conf_ints$upper
  ))
}

### testing
set.seed(12345)
test_df <- 
  tibble(
    pain_truth = sample(c("Yes", "No", "Nah"), 20, replace=TRUE), 
    pain_test  = sample(c("Yes", "No", "Nah"), 20, replace=TRUE)
  )
test_df
table(test_df)
get_sens_spec(test_df)


# the below function is for misspecification simulation later on

### pop_n_at_rand(): remove n integer values from a vector and return the
###                  edited vector and the popped values in a corresponding vector
###                  so that edited vector + popped vector == orig vector
# `vals`: assumed a vector of non-negative integers. e.g., c(6, 0, 8, 1)
# `n`: the number of integer values to randomly "pop" from `vals`
# `pop_vals`: a vector of values from `vals` that have been popped
# (default: NULL to denote it has to be intialised (NB: only used internally))
.pop_n_at_rand <- function(vals, n = 1, pop_vals = NULL) {
  
  if (any(vals < 0L)) {
    stop("No negative values allowed in `vals` vector")
  }
  if (sum(vals) < n) {
    stop("There is not enough integers in `val` to remove ", n, " integer(s)")
  }  
  if (is.null(pop_vals)) {
    pop_vals <- rep(0, length = length(vals))
    vals <- as.integer(vals)
    n <- as.integer(n)
    if (n < 1) { # only check on initialisation step of pop_vals (first call, not recursive ones)
      stop("`n` must be 1 or greater")
    }
  } else if (n < 1) {
    return(list(vals = vals, pop_vals = pop_vals))
  }
  
  poppable_vals_idx <- which(vals > 0) 
  n_poppable <- length(poppable_vals_idx)
  # the below error likely only to be invoked if some numeric/integer coercion
  # creates discordance with the other error checking
  if (n_poppable < 1) {
    stop("No positive values provided to pop a value off of")
  }
  rand_idx <- sample(poppable_vals_idx, 1)
  vals[rand_idx] <- vals[rand_idx] - 1
  pop_vals[rand_idx] <- pop_vals[rand_idx] + 1
  
  return(.pop_n_at_rand(vals, n = n - 1, pop_vals = pop_vals))
  
}

### wrapper function to remove `pop_vals` optional argument from function call
pop_n_at_rand <- function(vals, n = 1) {
  return(.pop_n_at_rand(vals, n = n))
}

### test error cases return correct errors 
expect_error(
  pop_n_at_rand(c(1, -1)), 
  "No negative values allowed in `vals` vector"
)
expect_error(
  pop_n_at_rand(rep(1, 3), n = 4), 
  "There is not enough integers in `val` to remove"
)
expect_error(
  pop_n_at_rand(rep(1, 3), n = -1), 
  "`n` must be 1 or greater"
)
expect_error(
  pop_n_at_rand(rep(1, 3), n = 0), 
  "`n` must be 1 or greater"
)
### non-error tests
pop_n_at_rand(rep(1, 3), n = 3)
pop_n_at_rand(1:5, n = 3)
pop_n_at_rand(1:5) # n = 1
pop_n_at_rand(c(6, 0, 8, 1), n = 8)
pop_n_at_rand(c(6, 0, 8, 1), n = 15)

### get_seq_cnts_from_cumul():
# if we have a vector of cumulative counts, return the original count vector
# ("sequential counts") before cumulative sum applied (inverse cumsum?)
get_seq_cnts_from_cumul <- function(cumul_n, as_prob = FALSE) {
  cumul_cnts <- cumul_n
  cnts <- c(cumul_cnts[1], diff(cumul_cnts))
  if (sum(cnts) != max(cumul_cnts)) {
    stop("the sequential counts do not sum to the max count (or negative counts exist)")
  }
  if (as_prob) {
    return(cnts / sum(cnts))
  }
  return(cnts)
}

### testing
set.seed(1234)
(these_vals <- sample(1:10))
(these_cumsum_vals <- cumsum(these_vals))
stopifnot(these_vals == get_seq_cnts_from_cumul(these_cumsum_vals))
# (these_vals <- sample((-5):5))
# (these_cumsum_vals <- cumsum(these_vals))
# stopifnot(these_vals == get_seq_cnts_from_cumul(these_cumsum_vals))

### int_matrix_to_named_vec():
# turn a matrix to column-major ordered vector keeping (row,col) combination
# names as vector names
# NB: NULL dim attributes and NA values are changed to integer indexes across dimension
int_matrix_to_named_vec <- function(mat) {
  vec <- as.integer(mat)
  rwnms <- dimnames(mat)[[1]]
  if (is.null(rwnms)) {
    rwnms <- 1:(dim(mat)[1])
  } else if (sum(is.na(rwnms)) > 0) {
    rwnms[is.na(rwnms)] <- which(is.na(rwnms))
  }
  clnms <- dimnames(mat)[[2]]
  if (is.null(clnms)) {
    clnms <- 1:(dim(mat)[2])
  } else if (sum(is.na(clnms)) > 0) {
    clnms[is.na(clnms)] <- which(is.na(clnms))
  }
  names(vec) <- 
    as.character(
      outer(rwnms, clnms, \(x, y) str_c("(", x, ",", y, ")"))
    )
  return(vec)
}
### testing
(test_mat_1 <- matrix(1:4, 2, dimnames = list(c(NA, NA), 1:2)))
int_matrix_to_named_vec(test_mat_1)
(test_mat_2 <- matrix(1:4, 2, dimnames = list(1:2, NULL)))
int_matrix_to_named_vec(test_mat_2)

# ---- import_truth_data ----


(truth_cols <- 
  cols(
    Report_ID = col_double(),
    `Pain topic` = col_character(),
    `Event description` = col_character()
  ))

pt_mg <- read_csv("dat/pain_truth_mg.csv", col_types = truth_cols)
pt_rl <- read_csv("dat/pain_truth_rl.csv", col_types = truth_cols)

(n_mg <- nrow(pt_mg))
(n_rl <- nrow(pt_rl))

pt_dat <-
  bind_rows(
    pt_mg %>% mutate(who = "mg"),
    pt_rl %>% mutate(who = "rl")
  ) %>% 
  rename(pain_truth = `Pain topic`) %>% 
  select(-`Event description`)
pt_dat


dbl_rws <-
  pt_dat %>% 
  group_by(Report_ID) %>% 
  summarise(n = n(), u_n = length(unique(pain_truth)), .groups = "drop")


dbl_rws %>% 
  filter(n > 1) 

dbl_rws %>% 
  filter(n > 1, u_n == 1) 

(n_dbl <- nrow(dbl_rws %>% filter(n > 1)))

pt_dat <-
  pt_dat  %>% 
  group_by(Report_ID) %>% 
  filter(row_number() == 1) %>% 
  ungroup() 

stopifnot(n_mg + n_rl - n_dbl == nrow(pt_dat))
  





# ---- import_topic_data ----


clean_data_cols <-
  cols(
    Report_ID = col_double(),
    Date = col_date(format = ""),
    pain_word = col_logical(),
    pain_topic = col_double(),
    type = col_character()
  )

clean_data <- read_csv("dat/clean_data.csv", col_types = clean_data_cols)

nrow(clean_data)
pt_dat <-
  clean_data %>% 
  left_join(
    .,
    pt_dat,
    "Report_ID"
  )
nrow(pt_dat)

with(pt_dat, table(pain_truth, useNA = "ifany"))







# ---- calc_sens_spec ----


thresh_use <- seq(0.010, 0.080, by = 0.01) #, 0.5, 1.01)
thresh_use_str <- sprintf("%0.3f", thresh_use)


thresh_pain_df <-
  foreach(i = seq_along(thresh_use), .combine = bind_rows) %do% {
    pt_dat_i <-
      pt_dat %>% 
      dplyr::filter(!is.na(pain_truth)) %>% 
      mutate(
        thresh = thresh_use_str[i], 
        pain_test = if_else(pain_topic >= thresh_use[i], "Yes", "No")
      )
  }

pt_dat %>% 
  dplyr::filter(!is.na(pain_truth)) %>% 
  mutate(
    thresh = "pain word", 
    pain_test = if_else(pain_word , "Yes", "No")
  )

thresh_pain_df <-
  thresh_pain_df %>% 
  mutate(
    type = fct_rev(fct_inorder(type)),
    thresh = fct_rev(fct_inorder(thresh))
  ) %>% 
  arrange(type, thresh)

thresh_pain_df


thresh_pain_df <-
  thresh_pain_df %>% 
  select(Report_ID, type, pain_truth, thresh, pain_test) %>% 
  nest(dat = -c(type, thresh))
 
thresh_pain_df

thresh_pain_df <-
  thresh_pain_df %>% 
  mutate(ss = map(dat, get_sens_spec))


thresh_pain_df


thresh_pain_df <-
  thresh_pain_df %>% 
  select(-dat) %>% 
  unnest(cols = c(ss)) %>% 
  mutate(
    param  = factor(param, levels = c("sens", "spec", "ppv", "npv"))
  ) %>% 
  arrange(type, param, thresh)


thresh_pain_df %>% 
  dplyr::filter(thresh %in% sprintf("%0.3f", 0.01 * c(4, 6, 8))) %>% 
  arrange(type, param, desc(thresh)) %>% 
  select(type, param, thresh, everything()) %>% 
  knitr::kable(., digits = 2)

(mesh_types <- rev(levels(thresh_pain_df$type)))

thresh_pain_df %>% 
  mutate(thresh = fct_rev(thresh)) %>% 
  ggplot(., aes(y = type, x = est, col = type)) +
  geom_errorbar(aes(xmin = lo, xmax = up), alpha = 0.5) + 
  geom_point(alpha = 0.95) + # aes(size = thresh)) +
  scale_x_continuous(limits = 0:1) +
  scale_y_discrete(limits = mesh_types) +
  scale_color_colorblind() +
  theme_bw() +
  facet_grid(rows = vars(thresh), cols = vars(param), labeller = label_both) +
  labs(y = "Mesh\ntype", x = "Estimate")


roc_dat <- 
  thresh_pain_df %>% 
  pivot_wider(
    # id = c(type, thresh),
    names_from = param,
    values_from = c(cases, correct,  est,   lo,   up),
    names_vary = "slowest"
    # names_glue = "{.name}_{.value}"
  )


roc_dat

roc_dat %>% 
  ggplot(., aes(y = est_sens, x = 1 - est_spec, col = type, group = type)) +
  geom_step(alpha = 0.5, orientation  = "y", direction = "hv") +  
  geom_point(alpha = 0.5) + # aes(size = thresh)) +
  scale_x_continuous(limits = 0:1) +
  scale_y_continuous(limits = 0:1) +
  scale_color_colorblind() +
  theme_bw() +
  coord_fixed(ratio = 1) +
  labs(x = "1 - specificity", y = "Sensitivity")



# ---- misclassification_effects  ----

### see:
# https://cran.r-project.org/web/packages/episensr/vignettes/episensr.html#misclassification-bias


### real qtrly data 
cumul_qtrly_dat <- read_parquet("dat/cumul_qtrly_dat.parquet") 
cumul_qtrly_dat

# distinct(cumul_qtrly_dat, grps, dat_type, thresh)
distinct(cumul_qtrly_dat, grps, dat_type)

test_cumul_qtrly_dat <-
  cumul_qtrly_dat %>% 
  filter(grps == "(b) pelvic_mesh v hernia_mesh/other_mesh", thresh == "0.050") 

test_cumul_qtrly_dat %>% filter(mnth == max(mnth))
# same, just checkin'
test_cumul_qtrly_dat %>% filter(row_number() == n())

### Note in our notation:
# nA == (target exposure, case) count
# nB == (target exposure, control) count
# nC == (comparator exposure, case) count
# nD == (comparator exposure, control) count

### while in `episensr`, the 2x2 table is transposed so that:
# a == (case, target exposure) count == nA
# b == (case, comparator exposure) count == nC
# c == (control, target exposure) count == nB
# d == (control, comparator exposure) count == nD

### therefore:
(pain_ex <-
  matrix(
    c(77, 25, 45, 85), 
    # i.e., c(nA, nB, nC, nD)  in our data; or 
    # c(a, c, b , d) in `episensr` notation
    dimnames = list(c("[pain]", "not [pain]"), c("pelvic", "hernia/other")),
    nrow = 2
  ))

pain_ex
(vec_pain_ex <- int_matrix_to_named_vec(pain_ex))

## now change counts (uniform prob) based on sims
test_cumul_qtrly_dat
sim_cumul_qtrly_dat <- test_cumul_qtrly_dat
(real_seq_a <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nA))
(real_seq_b <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nB))
(real_seq_c <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nC))
(real_seq_d <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nD))

seq_tots <- c(sum(real_seq_a), sum(real_seq_b), sum(real_seq_c), sum(real_seq_d))
names(seq_tots) <- names(vec_pain_ex)
seq_tots
expect_equal(seq_tots, vec_pain_ex)


misclass_pain_ex <-
  misclass(
    pain_ex,
    type = "outcome",
    ### Note when type = "outcome", the bias params are:
    # index_1 = Sensitivity of outcome classification among those with the exposure
    # index_2 = Sensitivity of outcome classification among those without the exposure
    # index_3 = Specificity of outcome classification among those with the exposure
    # index_4 = Specificity of outcome classification among those without the exposure
    bias_parms = c(0.99, 0.82, 0.78, 0.90)
  )

misclass_pain_ex

set.seed(123)
probsens_pain_ex <- 
  probsens(
    pain_ex,
    type = "outcome",
    reps = 1e5,
    # The sensitivity of outcome classification among those with the exposure
    seca = list("uniform", c(0.99)+ c(-0.01, 0.005)),
    # The sensitivity of outcome classification among those without the exposure
    seexp = list("uniform", c(0.82)+ c(-0.01, 0.01)),
    # The specificity of outcome classification among those with the exposure
    spca = list("uniform", c(0.78) + c(-0.01, 0.01)),
    # The specificity of outcome classification among those without the exposure
    spexp = list("uniform", c(0.90) + c(-0.01, 0.01)),
    corr_se = 0.2,
    corr_sp = 0.2
  )

probsens_pain_ex
# compare to analytical
# str(misclass_pain_ex$adj_measures)
misclass_pain_ex$adj_measures["   Misclassification Bias Corrected Odds Ratio:", " "]
median(probsens_pain_ex$sim_df$syst_OR)
median(probsens_pain_ex$sim_df$tot_OR)


tibble(probsens_pain_ex$sim_df)
sim_2x2 <- 
  probsens_pain_ex$sim_df %>% 
  select(all_of(c("ab", "cb", "bb", "db"))) %>% 
  rename(all_of(c("nA" = "ab", "nB" = "cb", "nC" = "bb", "nD" = "db")))
tibble(sim_2x2)

# check that pelvic exposure (a, b)-tuples from 2x2 table are the same total cases from orig 2x2
stopifnot(all(sum(pain_ex[, 1]) == rowSums(sim_2x2[, c("nA", "nB")])))
# check that hernia/other exposure (c, d)-tuples from 2x2 table are the same total cases from orig 2x2
stopifnot(all(sum(pain_ex[, 2]) == rowSums(sim_2x2[, c("nC", "nD")])))

# these are sim changes to contingency tabs (R is column-major operations on matrices)
i <- 5
str(sim_2x2)
sim_2x2[i, ] - vec_pain_ex
delta_tab <- t(apply(as.matrix(sim_2x2), 1, \(x) x - vec_pain_ex))
str(delta_tab)
stopifnot(all(0 == rowSums(delta_tab)))
head(delta_tab)


(delta_a_i <- delta_tab[i, "nA"])
(delta_b_i <- delta_tab[i, "nB"])
(delta_c_i <- delta_tab[i, "nC"])
(delta_d_i <- delta_tab[i, "nD"])

stopifnot(delta_a_i + delta_b_i == 0)
stopifnot(delta_c_i + delta_d_i == 0)

# initialise
seq_a_i <- real_seq_a
seq_b_i <- real_seq_b
seq_c_i <- real_seq_c
seq_d_i <- real_seq_d

# add random sim delta changes at random
if (delta_a_i < 0) {
  pop_vec_lst <- pop_n_at_rand(seq_a_i, -delta_a_i)
  seq_a_i <- pop_vec_lst$vals
  seq_b_i <- seq_b_i + pop_vec_lst$pop_vals
} else if (delta_b_i < 0) {
  pop_vec_lst <- pop_n_at_rand(seq_b_i, -delta_b_i)
  seq_a_i <- seq_a_i + pop_vec_lst$pop_vals
  seq_b_i <- pop_vec_lst$vals
}

if (delta_c_i < 0) {
  pop_vec_lst <- pop_n_at_rand(seq_c_i, -delta_c_i)
  seq_c_i <- pop_vec_lst$vals
  seq_d_i <- seq_d_i + pop_vec_lst$pop_vals
} else if (delta_d_i < 0) {
  pop_vec_lst <- pop_n_at_rand(seq_d_i, -delta_d_i)
  seq_c_i <- seq_c_i + pop_vec_lst$pop_vals
  seq_d_i <- pop_vec_lst$vals
}

c(delta_a_i, delta_b_i, delta_c_i, delta_d_i)
c(seq_a_i - real_seq_a)
c(seq_b_i - real_seq_b)
c(seq_c_i - real_seq_c)
c(seq_d_i - real_seq_d)



sim_cumul_qtrly_dat$nA <- cumsum(seq_a_i)
sim_cumul_qtrly_dat$nB <- cumsum(seq_b_i)
sim_cumul_qtrly_dat$nC <- cumsum(seq_c_i)
sim_cumul_qtrly_dat$nD <- cumsum(seq_d_i)


sim_cumul_qtrly_dat$sim_i <- as.integer(i)






perturb_real_by_misclass_sim <- function(real_df, sim_deltas_df) {
  
  real_seq_a <- get_seq_cnts_from_cumul(real_df[["nA"]])
  real_seq_b <- get_seq_cnts_from_cumul(real_df[["nB"]])
  real_seq_c <- get_seq_cnts_from_cumul(real_df[["nC"]])
  real_seq_d <- get_seq_cnts_from_cumul(real_df[["nD"]])
  
  n_sims <- nrow(sim_deltas_df)
  
  # sim_df <- NULL
  # for (i in 1:n_sims) {
    
  sim_df <-
    foreach(
      i = 1:n_sims,
      .combine = bind_rows,
      .export = c("real_df", "sim_deltas_df"),
      .errorhandling = "stop"
    ) %do% {
      
      tmp_real_df <- real_df
    
      (delta_a_i <- sim_deltas_df[i, "nA"])
      (delta_b_i <- sim_deltas_df[i, "nB"])
      (delta_c_i <- sim_deltas_df[i, "nC"])
      (delta_d_i <- sim_deltas_df[i, "nD"])
      
      # misclass only swaps outcome not exposure so parity within exposures remains
      if ((delta_a_i + delta_b_i) != 0) {
        stop("changes to target exposure case/control status not at parity")
      }
      if ((delta_c_i + delta_d_i) != 0) {
        stop("changes to comparator case/control status not at parity")
      }
      
      # initialise
      seq_a_i <- real_seq_a 
      seq_b_i <- real_seq_b
      seq_c_i <- real_seq_c
      seq_d_i <- real_seq_d
      
      # add random sim delta changes at random
      if (delta_a_i < 0) {
        pop_vec_lst <- pop_n_at_rand(seq_a_i, -delta_a_i)
        seq_a_i <- pop_vec_lst$vals
        seq_b_i <- seq_b_i + pop_vec_lst$pop_vals
      } else if (delta_b_i < 0) {
        pop_vec_lst <- pop_n_at_rand(seq_b_i, -delta_b_i)
        seq_a_i <- seq_a_i + pop_vec_lst$pop_vals
        seq_b_i <- pop_vec_lst$vals
      }
      
      if (delta_c_i < 0) {
        pop_vec_lst <- pop_n_at_rand(seq_c_i, -delta_c_i)
        seq_c_i <- pop_vec_lst$vals
        seq_d_i <- seq_d_i + pop_vec_lst$pop_vals
      } else if (delta_d_i < 0) {
        pop_vec_lst <- pop_n_at_rand(seq_d_i, -delta_d_i)
        seq_c_i <- seq_c_i + pop_vec_lst$pop_vals
        seq_d_i <- pop_vec_lst$vals
      }
      
      
      # check sequences have correct changes
      # print(c(sum(seq_a_i - real_seq_a), delta_a_i))
      if (sum(seq_a_i - real_seq_a) != delta_a_i) {
        stop("error in total changes made to target exposure case sequence")
      }
      if (sum(seq_b_i - real_seq_b) != delta_b_i) {
        stop("error in total changes made to target exposure control sequence")
      }
      if (sum(seq_c_i - real_seq_c) != delta_c_i) {
        stop("error in total changes made to comparator exposure case sequence")
      }
      if (sum(seq_d_i - real_seq_d) != delta_d_i) {
        stop("error in total changes made to comparator exposure control sequence")
      }
      
      tmp_real_df[["nA"]] <- cumsum(seq_a_i)
      tmp_real_df[["nB"]] <- cumsum(seq_b_i)
      tmp_real_df[["nC"]] <- cumsum(seq_c_i)
      tmp_real_df[["nD"]] <- cumsum(seq_d_i)
      
      tmp_real_df$sim_i <- as.integer(i)

      # sim_df <- bind_rows(sim_df, tmp_real_df)
      tmp_real_df

  }
  
  # return(relocate(sim_df, sim_i, .before = grps))
  return(sim_df)

}

# takes a minute or so to get through 10k == 1e+4 simulations
tic()
sim_cumul_qtrly_dat <-
  perturb_real_by_misclass_sim(
    real_df = test_cumul_qtrly_dat, 
    sim_deltas_df = delta_tab[1:1e+3, , drop = FALSE] ### testing starting small
    # sim_deltas_df = delta_tab
  )
toc()


# table(pull(distinct(sim_cumul_qtrly_dat, sim_i)) %in% 1:1e+4, useNA = "ifany")


## ---- disproporp_funcs ----

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





## ---- parallel_comp ----

# this seed can be set in future_map() etc for reproducible parallel comp seeds 
furrr_seed1 <- furrr_options(seed = 5678)
furrr_seed2 <- furrr_options(seed = 9012)
furrr_seed3 <- furrr_options(seed = 3456)
furrr_seed4 <- furrr_options(seed = 7890)

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




## ---- bcpnn_calcs ----




sra_cum <- sim_cumul_qtrly_dat

# make data for each combination of params nested for purrr like processing
sra_cum <-
  sra_cum %>%
  nest(data = c(mnth, nA, nB, nC, nD))



# testing/example
sra_cum$data[[5]] %>% print(., n = nrow(.))
get_sig_tab_over_time(sra_cum$data[[5]])



### for 1k sims with r9-5900X/128GB 2133mhz memory 
# takes ~ 90 sec [22 threads @ 4.4Ghz and 145W power draw]
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
sra_cum$sig_tab[[5]]


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
  group_by(grps, dat_type, thresh, sim_i) %>%
  dplyr::filter(ci_lo > 0) %>%
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(sra_cum_bcpnn)
sra_cum_bcpnn <-
  left_join(
    sra_cum_bcpnn,
    bcpnn_signif %>% select(grps, dat_type, thresh, sim_i, dte_reach_sig),
    c("grps", "dat_type", "thresh", "sim_i")
  )
nrow(sra_cum_bcpnn)

sra_cum_bcpnn


sra_cum_bcpnn <- 
  sra_cum_bcpnn %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )



## ---- save1 ----


sra_cum_bcpnn %>%
  write_parquet(., sink = "out/sra_cum_bcpnn_sim.parquet")




## ---- multcompar_bcpnn ----

# sra_cum <- 
#   sra_dat %>%
#   dplyr::filter(dat_type == "cumulative") 
sra_cum <- sim_cumul_qtrly_dat


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
  group_by(grps, dat_type, thresh, sim_i) %>%
  dplyr::filter(ci_lo > 0) %>%
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(sra_cum_bcpnn_mc_adj)
sra_cum_bcpnn_mc_adj <-
  left_join(
    sra_cum_bcpnn_mc_adj,
    bcpnn_mc_adj_signif %>% select(grps, dat_type, thresh, sim_i, dte_reach_sig),
    c("grps", "dat_type", "thresh", "sim_i")
  )
nrow(sra_cum_bcpnn_mc_adj)

sra_cum_bcpnn_mc_adj


sra_cum_bcpnn_mc_adj <- 
  sra_cum_bcpnn_mc_adj %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )





## ---- save2 ----


sra_cum_bcpnn_mc_adj %>%
  write_parquet(., sink = "out/sra_cum_bcpnn_sim_mc_adj.parquet")



## ---- plot_sim ----


# sra_cum_prr_sim_mc_adj <-  read_parquet("out/sra_cum_prr_sim_mc_adj.parquet")
sra_cum_bcpnn_sim_mc_adj <-  read_parquet("out/sra_cum_bcpnn_sim_mc_adj.parquet")
# sra_cum_maxsprt_sim <-  read_parquet("out/sra_cum_maxsprt_sim.parquet")


sra <-
  bind_rows(
    sra_cum_bcpnn_sim_mc_adj %>% mutate(stat = "BCPNN (MCadj)"),
    # sra_cum_prr_mc_adj %>% mutate(stat = "PRR (MCadj)"),
    # sra_cum_maxsprt %>% mutate(stat = "maxSPRT")
  ) %>%
  mutate(stat = fct_inorder(stat)) %>%
  select(stat, everything())

# sra %>%
#   dplyr::filter(stat == "BCPNN (MCadj)")


sra <-
  sra %>%
  mutate(
    test_stat = 
      case_when(
        # stat == "maxSPRT"     ~ maxllr, 
        # stat == "PRR (MCadj)" ~ ci_lo, 
        stat == "BCPNN (MCadj)" ~ ci_lo
      ),
    test_thresh = 
      case_when(
        # stat == "maxSPRT"     ~ cv, 
        # stat == "PRR (MCadj)" ~ 1, 
        stat == "BCPNN (MCadj)" ~ 0
      ),
    rr_stat = 
      case_when(
        # stat == "maxSPRT"     ~ rre, 
        # stat == "PRR (MCadj)" ~ est, 
        stat == "BCPNN (MCadj)" ~ 2 ^ est
      )
  )




# thresholds <- sort(unique(sra[["thresh"]]))
# length(thresholds)





plt_dat <-
  bind_rows(
    sra_cum_bcpnn_mc_adj %>% 
      mutate(
        cv = 0, 
        stat = "IC (BCPNN, Lower 95% CI)"
      ) %>%
      select(stat, grps, thresh, sim_i, dte, cv, val = ci_lo, reach_sig, dte_reach_sig) # ,
    # sra_cum_prr_mc_adj %>% 
    #   mutate(
    #     cv = 1, 
    #     stat = "RR (PRR, Lower 95% CI)"
    #   ) %>%
    #   select(stat, grps, thresh, sim_i, dte, cv, val = ci_lo, reach_sig, dte_reach_sig),
    # sra_cum_maxsprt %>%
    #   mutate(stat = "MaxSPRT (max LLR)") %>%
    #   select(stat, grps, thresh, sim_i, dte, cv, val = maxllr, reach_sig, dte_reach_sig)
  ) 


sig_reach_dat <-
  plt_dat %>%
  arrange(stat, grps, thresh, sim_i, dte) %>%
  group_by(stat, grps, thresh, sim_i) %>%
  dplyr::filter(reach_sig == 1) %>%
  dplyr::filter(row_number() == 1) %>%
  select(stat, grps, thresh, sim_i, dte_reached = dte) %>%
  # now create separation between reached CV values when it occurs
  group_by(stat, thresh, sim_i, dte_reached) %>%
  mutate(rep_dte = 1:n()) %>%
  ungroup() %>%
  mutate(dte_reached = dte_reached + days(10 * (rep_dte - 1))) %>%
  select(-rep_dte)





plt_dat <-
  left_join(
    plt_dat,
    sig_reach_dat,
    c("stat", "grps", "thresh", "sim_i")
  ) # 

plt_dat %>%
 dplyr::filter(dte_reach_sig != dte_reached)


plt_dat <-
  plt_dat %>%
  mutate(
    stat = fct_inorder(stat)
  )


# plt_dat %>%
#   dplyr::filter(thresh == "0.050", grepl("(c)", grps, fixed = TRUE))




# ---- time_to_sig_plot1 ----


col_pal <- c("darkorange", "cyan4", "purple")


date_signif_dat <-
  sra %>%
  group_by(stat, grps, dat_type, thresh, sim_i) %>%
  arrange(dte) %>%
  dplyr::filter(reach_sig) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(stat, grps, dat_type, thresh, sim_i) 


signif_plt <-
  date_signif_dat %>%
  ### only keep pelvic mesh as target vs whatever comparator
  dplyr::filter(grepl("^.*pelvic.* v ", grps)) %>%
  mutate(
    grps = gsub("\\([a-z]\\) ", "", grps),
    grps = gsub("_", " ", grps),
    grps = gsub("pelvic mesh", "Pelvic mesh", grps),
    grps = gsub("hernia mesh", "Hernia mesh", grps),
    # grps = str_to_sentence(grps),
    grps = gsub(" v ", "\nv\n", grps, fixed = TRUE),
    grps = fct_inorder(grps)
    # stat = fct_inorder(stat)
  ) 



# levels(signif_plt$grps)

signif_plt %>%
  arrange(grps, thresh) %>%
  ggplot(., aes(x = dte_reach_sig, fill = stat)) +
  geom_bar(alpha = 0.8, col = NA) +
  # geom_path(aes(group = interaction(stat, sim_i))) +
  scale_fill_manual(values = col_pal) +
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
  geom_vline(xintercept = quantile(signif_plt$dte_reach_sig, c(0.01, 0.5, 0.99), type = 3)) +
  facet_grid(stat ~ grps) +
  # facet_wrap( ~ grps, ncol = 1) +
  theme_bw() +
  theme(text = element_text(family = "serif")) +
  labs(
    x = expression("Date" ~ H[0] ~ "rejected (null hypothesis of no signal)"),
    y = "Threshold for P(topic = 'pain' | doc) to be classed a 'pain' adverse event",
    # col = "Signal detection\nmethod",
    fill = "Signal detection\nmethod"
  )

# signif_plt %>%
#   arrange(grps, thresh) %>%
#   ggplot(., aes(x = dte_reach_sig, y = as.numeric(thresh), col = stat)) +
#   geom_point(alpha = 0.2) +
#   # geom_path(aes(group = interaction(stat, sim_i))) +
#   scale_colour_manual(values = col_pal) +
#   # scale_colour_tableau(palette = "Color Blind", direction = -1) +
#   facet_grid(stat ~ grps) +
#   # facet_wrap( ~ grps, ncol = 1) +
#   theme_bw() +
#   theme(text = element_text(family = "serif")) +
#   labs(
#     x = expression("Date" ~ H[0] ~ "rejected (null hypothesis of no signal)"),
#     y = "Threshold for P(topic = 'pain' | doc) to be classed a 'pain' adverse event",
#     col = "Signal detection\nmethod"
#   )

# ggsave(
#   filename = "fig/time_to_signal_method_facets.png", 
#   dpi = 900, width = 9, height = 9
# )
ggsave(
  filename = "fig/sim_data_time_to_signal_method_facets.pdf",
  device = cairo_pdf, # embed fonts
  width = 9, height = 9, units = "in"
)








probsens_pain_ex <- 
  probsens(
    pain_ex,
    type = "exposure",
    reps = 50000,
    # mode(beta) = (a - 1) / (a + b - 2)
    seca = list("beta", c(25, 3)),
    spca = list("uniform", c(.9, .99))
  )
plot(probsens_pain_ex, "seca")
plot(probsens_pain_ex, "spca")








# ---- composite_reference_standard-crs ----

### composite reference standard 
# for when not all tests have an associated definite truth value

### package not used here but can be
# remotes::install_github("tystan/crs")
# library("crs")
# data(brenton2019)
# brenton2019

### these mappings of values make indexing values easier below
# "Yes" ~ as.integer(2)
# "No"  ~ as.integer(1)

mk_ny_12 <- function(x) {
  case_when(
    is.na(x)   ~ as.integer(NA),
    x == "Yes" ~ as.integer(2),
    x == "No"  ~ as.integer(1),
    TRUE       ~ -999L
  )
}

# test fn
table(test_df$pain_truth, mk_ny_12(test_df$pain_truth))

# required foirmat for (index, imperfect, resolver)-tuple
(test_pt_dat <- 
  pt_dat %>%  
  mutate(
    pain_test_index = if_else(pain_topic >= 0.05, "Yes", "No"),
    pain_test_imperfect = if_else(pain_topic >= 0.0001, "Yes", "No"), 
    pain_test_resolver = pain_truth
  ))

# our meshs of interest are all resolved (no CRS required)
test_pt_dat %>% 
  dplyr::filter(type %in% str_c(c("other", "hernia", "pelvic"), "_mesh")) %>% 
  with(
    ., 
    table(
      pain_test_index, 
      pain_test_imperfect, 
      pain_test_resolver,
      useNA = "ifany"
    )
  )

# our non-meshs have no resolver values (also no CRS required)
test_pt_dat %>% 
  dplyr::filter(!(type %in% str_c(c("other", "hernia", "pelvic"), "_mesh"))) %>% 
  with(
    ., 
    table(
      pain_test_index, 
      pain_test_imperfect, 
      pain_test_resolver,
      useNA = "ifany"
    )
  )

# values to cyclye through below
skeleton_grid <- 
  thresh_pain_df %>% 
  dplyr::filter(thresh %in% sprintf("%0.3f", 0.01 * c(4, 6, 8))) %>% 
  distinct(type, thresh)

# do CRS anyway to see if we get the same (sens, spec) values back. We do!!
foreach(i = seq_len(nrow(skeleton_grid)), .combine = bind_rows) %do% {
  
  type_i <- skeleton_grid[["type"]][i]
  thresh_i <- as.numeric(as.character(skeleton_grid[["thresh"]][i]))
  crs_dat <-
    pt_dat %>% 
    dplyr::filter(type == type_i) %>% 
    mutate(
      pain_test_index = if_else(pain_topic >= thresh_i, "Yes", "No"),
      pain_test_imperfect = if_else(pain_topic >= 0.0001, "Yes", "No"), 
      pain_test_resolver = pain_truth
    ) %>% 
    select(pain_test_index, pain_test_imperfect, pain_test_resolver)
  
  
  crs_dat <-
    crs_dat %>% 
    mutate(across(everything(), mk_ny_12))
  
  crs_dat %>% 
    dplyr::filter(!is.na(pain_test_resolver)) %>% 
    with(., table(pain_test_index, pain_test_imperfect, pain_test_resolver))
  
  marg <- with(crs_dat, table(pain_test_index, pain_test_imperfect))
  rslv <- 
    crs_dat %>% 
    dplyr::filter(!is.na(pain_test_resolver)) %>% 
    with(., table(pain_test_index, pain_test_imperfect, pain_test_resolver))
  
  n <- sum(marg)
  m_ij <- apply(rslv, 1:2, sum)
  m <- sum(m_ij)
  p_ij_dot <- marg / n
  r_ij <- rslv[, , 2] / (rslv[, , 1] +  rslv[, , 2])
  r_ij[is.nan(r_ij)] <- 0
  
  p_ij1 <- p_ij_dot * (1 - r_ij)
  p_ij2 <- p_ij_dot * r_ij
  
  true_pos  <- p_ij2[2, 2] + p_ij2[2, 1]
  false_neg <- p_ij2[1, 2] + p_ij2[1, 1]
  
  true_neg  <- p_ij1[1, 1] + p_ij1[1, 2]
  false_pos <- p_ij1[2, 2] + p_ij1[2, 1]
  
  (crs_ests <-
    tibble(
      type = type_i,
      theresh = thresh_i,
      sens = true_pos / (true_pos + false_neg),
      spec = true_neg / (true_neg + false_pos),
      ppv = true_pos / (true_pos + false_pos),
      npv = true_neg / (true_neg + false_neg)
    ))

}













