


# ---- libs ----


library("arrow")   # parquet files
library("binom")   # wilson confidence intervals for binomial counts
library("foreach") # flexible looping and return amalgamation
library("episensr") # misclassification error contingency table adjustments
library("simdata") # NORTA method to get correlation estimates 
# citation("episensr")

# tidyverse et al
library("readr")
library("dplyr")
library("purrr")
library("tidyr")
library("forcats")
library("stringr")
library("ggplot2")
library("ggthemes")




# ---- func ----

get_sens_spec <- function (df, pos = "Yes", alpha = 0.05) {
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
set.seed(12345)
test_df <- 
  tibble(
    pain_truth = sample(c("Yes", "No", "Nah"), 20, replace=TRUE), 
    pain_test  = sample(c("Yes", "No", "Nah"), 20, replace=TRUE)
  )
test_df
table(test_df)
get_sens_spec(test_df)



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







# ---- import_topic_data ----


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

test_cumul_qtrly_dat%>% filter(mnth == max(mnth))
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
pain_ex <-
  matrix(
    c(77, 45, 25, 85), # i.e., c(nA, nC, nB, nD) or c(a, b, c, d)
    dimnames = list(c("[pain]", "not [pain]"), c("pelvic", "hernia/other")),
    nrow = 2, 
    byrow = TRUE
  )

missclass_pain_ex <-
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

missclass_pain_ex

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
# str(missclass_pain_ex$adj_measures)
missclass_pain_ex$adj_measures["   Misclassification Bias Corrected Odds Ratio:", " "]
median(probsens_pain_ex$sim_df$syst_OR)
median(probsens_pain_ex$sim_df$tot_OR)


tibble(probsens_pain_ex$sim_df)
sim_2x2 <- probsens_pain_ex$sim_df[, c("ab", "bb", "cb", "db")]
tibble(sim_2x2)
pain_ex
# check that (a, b)-tuples from 2x2 table are the same total cases from orig 2x2
stopifnot(all(sum(pain_ex[, 1]) == rowSums(sim_2x2[, c(1, 3)])))
# check that (c, d)-tules from 2x2 table are the same total cases from orig 2x2
stopifnot(all(sum(pain_ex[, 2]) == rowSums(sim_2x2[, c(2, 4)])))

# these are sim changes to contingency tabs (R is column-major operations on matrices)
str(t(sim_2x2))
t(sim_2x2)[, 5] - as.vector(t(pain_ex))
t(as.matrix(sim_2x2))
delta_tab <- t(t(sim_2x2) - as.vector(t(pain_ex)))
head(delta_tab)
stopifnot(all(0 == rowSums(delta_tab)))

## now change counts (uniform prob) based on sims
test_cumul_qtrly_dat

get_seq_cnts_from_cumul <- function(cumul_n, as_prob = FALSE) {
  cumul_cnts <- cumul_n
  cnts <- c(cumul_cnts[1], diff(cumul_cnts))
  stopifnot(sum(cnts) == max(cumul_cnts))
  if (as_prob) {
    return(cnts / sum(cnts))
  }
  return(cnts)
}
(delta_a <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nA))
sum(delta_a)
(delta_b <- get_seq_cnts_from_cumul(test_cumul_qtrly_dat$nC))
sum(delta_b)


i <- 1
(delta_a_i <- delta_tab[i, "ab"])
(delta_b_i <- delta_tab[i, "bb"])
(delta_c_i <- delta_tab[i, "cb"])
(delta_d_i <- delta_tab[i, "db"])


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








# ---- crs ----

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













