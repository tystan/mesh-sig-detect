


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
library("ggpubr")




# ---- func ----

# here are the functions written for these analyses
# they will be shown in the *Appendix A*
source("r/_funcs.R")



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



# ---- misclassification_effects ----

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


thresh_pain_df %>% 
  dplyr::filter(thresh %in% sprintf("%0.3f", 0.05), param %in% c("sens", "spec")) %>% 
  arrange(type, param, desc(thresh)) %>% 
  select(type, param, thresh, everything()) %>% 
  knitr::kable(., digits = 2)

### trying to replicate this approximately
  # |type        |param |thresh | cases| correct|  est|   lo|   up|
  # |:-----------|:-----|:------|-----:|-------:|----:|----:|----:|
  # |pelvic_mesh |sens  |0.050  |    70|      70| 1.00| 0.95| 1.00|
  # |pelvic_mesh |spec  |0.050  |    32|      25| 0.78| 0.61| 0.89|
  # |other_mesh  |sens  |0.050  |    31|      28| 0.90| 0.75| 0.97|
  # |other_mesh  |spec  |0.050  |    53|      44| 0.83| 0.71| 0.91|

misclass_pain_ex <-
  misclass(
    pain_ex,
    type = "outcome",
    ### Note when type = "outcome", the bias params are:
    # index_1 = Sensitivity of outcome classification among those with the exposure
    # index_2 = Sensitivity of outcome classification among those without the exposure
    # index_3 = Specificity of outcome classification among those with the exposure
    # index_4 = Specificity of outcome classification among those without the exposure
    bias_parms = c(0.99, 0.90, 0.78, 0.83)
  )

misclass_pain_ex



# note the mode of beta dist is (a - 1) / (a + b - 2)
# also the larger a + b the smaller the variance

par(mfrow = c(2, 2))
hist(rbeta(1e+5, 100, 1), main = "pelvic_mesh | sens", xlim = 0:1) # pelvic_mesh |sens| 1.00| 0.95| 1.00|
hist(rbeta(1e+5, 60, 16), main = "pelvic_mesh | spec", xlim = 0:1) # pelvic_mesh |spec| 0.78| 0.61| 0.89|
hist(rbeta(1e+5, 81, 9), main = "other_mesh | sens", xlim = 0:1) # other_mesh |sens| 0.90| 0.75| 0.97|
hist(rbeta(1e+5, 100, 21), main = "other_mesh | spec", xlim = 0:1) # other_mesh |spec| 0.83| 0.71| 0.91|
par(mfrow = c(1, 1))

set.seed(123)
probsens_pain_ex <- 
  probsens(
    pain_ex,
    type = "outcome",
    reps = 1e+5,
    # The sensitivity of outcome classification among those with the exposure
    seca = list("beta", c(100, 1)),
    # The specificity of outcome classification among those with the exposure
    spca = list("beta", c(60, 16)),
    # The sensitivity of outcome classification among those without the exposure
    seexp = list("beta", c(81, 9)),
    # The specificity of outcome classification among those without the exposure
    spexp = list("beta", c(100, 21)),
    corr_se = 0.25, # Correlations should be > 0 and < 1 (BUT more than ~0.3 brings errors [infeasibility])
    corr_sp = 0.25  # [from testing these values almost make no difference]
  )

### testing: using uniform distribution with small domain ~= point estimates. 
### I think that's clever anyway:-)
# set.seed(123)
# probsens_pain_ex <- 
#   probsens(
#     pain_ex,
#     type = "outcome",
#     reps = 1e5,
#     # The sensitivity of outcome classification among those with the exposure
#     seca = list("uniform", c(0.99)+ c(-0.01, 0.005)),
#     # The specificity of outcome classification among those with the exposure
#     spca = list("uniform", c(0.78) + c(-0.01, 0.01)),
#     # The sensitivity of outcome classification among those without the exposure
#     seexp = list("uniform", c(0.82)+ c(-0.01, 0.01)),
#     # The specificity of outcome classification among those without the exposure
#     spexp = list("uniform", c(0.90) + c(-0.01, 0.01)),
#     corr_se = 0.25,
#     corr_sp = 0.25
#   )

brk_inc <- 0.02
first_lab <- 0.1
breaks_0_1 <- seq(0 - brk_inc, 1 + brk_inc, brk_inc)
labels_0_1 <- breaks_0_1
labels_0_1[-seq(1 + 1, length(labels_0_1) - 1, round(first_lab / brk_inc))] <- ""

add_plt_extras <- function(plt_obj, which_param = 1, rm_dens_geom = TRUE) {
  if (rm_dens_geom) {
    plt_obj$layers[["geom_density"]] <- NULL
  }
  labs_vec <- 
    as.character(
      outer(
        c("Exposed", "Not exposed"), 
        c("Sensitivity", "Specificity"), 
        FUN = \(x, y) str_c(y, x, sep = " of ")
      )
    )
  plt_obj +
    scale_x_continuous(
      limits = c(-0.01, 1.01), 
      breaks = breaks_0_1, 
      labels = labels_0_1
    ) +
    labs(x = labs_vec[which_param], y = "Density of simulations") +
    theme_classic()
}

p_sens_exp <- plot(probsens_pain_ex, parms = "seca") %>% add_plt_extras(., 1)
p_spec_exp <-plot(probsens_pain_ex, parms = "spca") %>% add_plt_extras(., 2)
p_sens_nonexp <- plot(probsens_pain_ex, parms = "seexp") %>% add_plt_extras(., 3)
p_spec_nonexp <- plot(probsens_pain_ex, parms = "spexp") %>% add_plt_extras(., 4)

ggarrange(p_sens_exp, p_sens_nonexp, p_spec_exp, p_spec_nonexp, nrow = 2, ncol = 2)


probsens_pain_ex

# str(misclass_pain_ex$adj_measures)
### this is the original simple and analytical calculatiopns for misclassification
misclass_pain_ex$adj_measures["   Misclassification Bias Corrected Odds Ratio:", " "]
### compare to more complex beta priors (same median estimates roughly so working well)
### NB: the beta priors will induce more variability than the simple approach to REALLY "kick the tyres"
###     of the time-to-disproportionate-significance
median(probsens_pain_ex$sim_df$syst_OR) # very similar to the simple misclass analysis
median(probsens_pain_ex$sim_df$tot_OR)  # very similar to the simple misclass analysis too



# ---- example_induced_misclassification_effects_to_real_data ----

# re-use pain example from above
pain_ex
(vec_pain_ex <- int_matrix_to_named_vec(pain_ex))

## now change counts (uniform prob) based on sims
test_cumul_qtrly_dat
sim_cumul_qtrly_dat <- test_cumul_qtrly_dat
(real_seq_a <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nA))
(real_seq_b <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nB))
(real_seq_c <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nC))
(real_seq_d <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nD))

### checks
seq_tots <- c(sum(real_seq_a), sum(real_seq_b), sum(real_seq_c), sum(real_seq_d))
names(seq_tots) <- names(vec_pain_ex)
seq_tots
expect_equal(seq_tots, vec_pain_ex)

### these are the sims to induce on the real example
tibble(probsens_pain_ex$sim_df)
sim_2x2 <- 
  probsens_pain_ex$sim_df %>% 
  select(all_of(c("ab", "cb", "bb", "db"))) %>% 
  rename(all_of(c("nA" = "ab", "nB" = "cb", "nC" = "bb", "nD" = "db")))
tibble(sim_2x2) # each row a simulation

# check that pelvic exposure (a, b)-tuples from 2x2 table are the same total cases from orig 2x2
stopifnot(all(sum(pain_ex[, 1]) == rowSums(sim_2x2[, c("nA", "nB")])))
# check that hernia/other exposure (c, d)-tuples from 2x2 table are the same total cases from orig 2x2
stopifnot(all(sum(pain_ex[, 2]) == rowSums(sim_2x2[, c("nC", "nD")])))

# these are sim changes to contingency tabs (R is column-major operations on matrices)
i <- 5 # use the fifth sim as more interesting
str(sim_2x2)
sim_2x2[i, ] - vec_pain_ex # these are the simulated changes from the real example
delta_tab <- t(apply(as.matrix(sim_2x2), 1, \(x) x - vec_pain_ex))
str(delta_tab)
stopifnot(all(0 == rowSums(delta_tab)))
head(delta_tab)


# final time estimate for real example 
example_as_tib <- as_tibble(t(vec_pain_ex)) 
colnames(example_as_tib) <- str_c("n", LETTERS[1:4])
example_as_tib
example_as_tib %>% 
  mutate(adj_alpha = 0.1) %>%
  get_sig_tab_over_time_2(.)

# how does the sim version change the final time estimate? 
### NB much higher estimate/lo_ci as expected with less comparator outcomes
sim_2x2[i, ] 
sim_2x2[i, ] %>% 
  mutate(adj_alpha = 0.1) %>%
  get_sig_tab_over_time_2(.)


(delta_a_i <- delta_tab[i, "nA"])
(delta_b_i <- delta_tab[i, "nB"])
(delta_c_i <- delta_tab[i, "nC"])
(delta_d_i <- delta_tab[i, "nD"])

### checks logic holds
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

### checks again for fun
c(delta_a_i, delta_b_i, delta_c_i, delta_d_i)
c(seq_a_i - real_seq_a)
c(seq_b_i - real_seq_b)
c(seq_c_i - real_seq_c)
c(seq_d_i - real_seq_d)

# re-insert the simulated changes back into data structure
sim_cumul_qtrly_dat$nA <- cumsum(seq_a_i)
sim_cumul_qtrly_dat$nB <- cumsum(seq_b_i)
sim_cumul_qtrly_dat$nC <- cumsum(seq_c_i)
sim_cumul_qtrly_dat$nD <- cumsum(seq_d_i)
# label the simulation too
sim_cumul_qtrly_dat$sim_i <- as.integer(i)

### test analysis
sim_cumul_qtrly_dat %>% 
  mutate(adj_alpha = 0.1) %>%
  get_sig_tab_over_time_2(.)

# final time estimate for real example the same? (YES!!!)
example_as_tib %>% 
  mutate(adj_alpha = 0.1) %>%
  get_sig_tab_over_time_2(.)



# ---- funcs_for_induced_misclassification_effects_to_real_data ----



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

### use a subset of sims because of computational limitations at the mo
# (and seems not to linearly scale currently)
sub_sim <- 1e+3

# takes a minute or so to get through 10k == 1e+4 simulations
tic()
sim_cumul_qtrly_dat <-
  perturb_real_by_misclass_sim(
    real_df = test_cumul_qtrly_dat, 
    sim_deltas_df = delta_tab[1:sub_sim, , drop = FALSE] 
    # sim_deltas_df = delta_tab
  )
toc()


# table(pull(distinct(sim_cumul_qtrly_dat, sim_i)) %in% 1:1e+4, useNA = "ifany")





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

### for 1k sims with r9-5900X/128GB 2133mhz memory 
# takes ~ 90 sec [22 threads @ 4.4Ghz and 145W power draw]
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





## ---- save_bcpnn_sim_mc ----


sra_cum_bcpnn_mc_adj %>%
  write_parquet(., sink = "out/sra_cum_bcpnn_sim_mc_adj.parquet")







# ---- multcompar_prr ----


sra_cum <- sim_cumul_qtrly_dat



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


### takes ~5 sec 
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
sra_cum$sig_tab[[1]]


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


# deal with warnings about 0 counts that affects ests and CIs
sra_cum_prr_mc_adj <-
  sra_cum_prr_mc_adj %>%
  mutate(
    est = if_else(!is.finite(est), as.numeric(1), est),
    ci_lo = if_else(!is.finite(ci_lo), -Inf, ci_lo),
  )

sra_cum_prr_mc_adj

with(sra_cum_prr_mc_adj, table(dte, mnth, useNA = "ifany")) %>% 
  as.data.frame() %>%
  dplyr::filter(Freq > 0) %>%
  arrange(mnth, dte)


# first signif
prr_mc_adj_signif <-
  sra_cum_prr_mc_adj %>%
  group_by(grps, dat_type, thresh, sim_i) %>%
  dplyr::filter(ci_lo > 1) %>% # 1 is the critical value on ratio scale
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(sra_cum_prr_mc_adj)
sra_cum_prr_mc_adj <-
  left_join(
    sra_cum_prr_mc_adj,
    prr_mc_adj_signif %>% select(grps, dat_type, thresh, sim_i, dte_reach_sig),
    c("grps", "dat_type", "thresh", "sim_i")
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
  dplyr::filter(thresh == "0.050", grepl("(b)", grps, fixed = TRUE))


# sra_cum_prr_mc_adj %>%
#   dplyr::filter(thresh == "0.050", grepl("(b)", grps, fixed = TRUE)) %>%
#   print(., n = nrow(.))

sra_cum_prr_mc_adj <- 
  sra_cum_prr_mc_adj %>%
  mutate(
    dte_reach_sig = if_else(is.na(dte_reach_sig), as_date(today()), dte_reach_sig),
    reach_sig = dte >= dte_reach_sig
  )





# ---- save_prr_mc ----


sra_cum_prr_mc_adj %>%
  write_parquet(., sink = "out/sra_cum_prr_sim_mc_adj.parquet")





# ---- maxsprt ----


sra_cum <- sim_cumul_qtrly_dat


cv_tab <-
  sra_cum %>%
  dplyr::filter(sim_i == 1) %>%
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
  knitr::kable(., digits = 1)


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
  group_by(grps, dat_type, thresh, sim_i) %>%
  dplyr::filter(reached_cv > 0) %>%
  arrange(dte) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  rename(dte_reach_sig = dte)


nrow(maxsprt_dat)
maxsprt_dat <-
  left_join(
    maxsprt_dat,
    maxsprt_signif %>% select(grps, dat_type, thresh, sim_i, dte_reach_sig),
    c("grps", "dat_type", "thresh", "sim_i")
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
  arrange(grps, thresh, sim_i, modifier, mult, alt_str, mnth) %>%
  select(grps, thresh, sim_i, modifier, mult, alt_str, everything())
nrow(maxsprt_dat_alts)


if(nrow(maxsprt_dat_alts) != (nrow(alt_mults) + 2) * nrow(maxsprt_dat_calcs)) {
  stop("many-to-many join has not worked")
}

print(maxsprt_dat_alts, n = 30)


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
  group_by(grps, dat_type, thresh, sim_i, modifier, mult, alt_str) %>%
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
      select(grps, dat_type, thresh, sim_i, modifier, mult, alt_str, dte_reach_sig),
    c("grps", "dat_type", "thresh", "sim_i", "modifier", "mult", "alt_str")
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





# ---- save_maxsprt_sim ----


maxsprt_dat %>%
  write_parquet(., sink = "out/sra_cum_maxsprt_sim.parquet")





maxsprt_dat_alts %>%
  write_parquet(., sink = "out/sra_cum_maxsprt_sim_alt_cvs.parquet")







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













