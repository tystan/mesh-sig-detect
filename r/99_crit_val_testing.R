
# ---- libs ----


library("dplyr")
library("tidyr")
library("tictoc")
library("foreach")
library("EmpiricalCalibration")
library("Sequential")
library("arrow")


# ---- consts ----

min_event <- 1
alpha <- 0.05

cv_df_0 <- 
  expand_grid(
    cntl_to_case_ratio = seq(0.5, 5, by = 0.5),
    max_n = 200, 
    look_interval = c(1:5, 8, 10)
  )

### testing
# cv_df_0 <- 
#   expand_grid(
#     cntl_to_case_ratio = c(0.5, 2, 5),
#     max_n = 200, 
#     look_interval = c(5, 10)
#   )


# ---- seq_package ----


### testing Sequential::CV.Binomial()

cv_df_seq <- cv_df_0
cv_df_seq$cv <- 0
cv_df_seq$pack <- "Sequential"
cv_df_seq



tic()
for (i in 1:nrow(cv_df_seq)) {
  
  look_i <- cv_df_seq$look_interval[i] 
  max_n_i <- cv_df_seq$max_n[i]
  z_i <- cv_df_seq$cntl_to_case_ratio[i]
  
  # has to be cumulative for Sequential::CV.Binomial()
  # i.e. test performed at 3 events then when 3 more events come in requires
  # GroupSizes = c(3, 6)
  gs_seq <- rep(look_i, floor(max_n_i / look_i))
  if (sum(gs_seq) != max_n_i) { # if doesn't go to max_n, add at end for last look
    gs_seq <- c(gs_seq, max_n_i - sum(gs_seq)) 
  }
  
  cv_df_seq$cv[i] <- 
    Sequential::CV.Binomial(
      N = max_n_i,
      alpha = alpha,
      M = min_event,
      z = z_i, 
      GroupSizes = gs_seq
    )$cv
  
}
toc()

cv_df_seq


cv_df_seq %>%
  write_parquet(., sink = "out/cv_df_seq.parquet")


# ---- empcalib_package ----


### testing EmpiricalCalibration::computeCvBinomial()

cv_df_ec <- cv_df_0
cv_df_ec$cv <- 0
cv_df_ec$pack <- "EmpiricalCalibration"
cv_df_ec


tic()
for (i in 1:nrow(cv_df_ec)) {
  
  look_i <- cv_df_ec$look_interval[i] 
  max_n_i <- cv_df_ec$max_n[i]
  z_i <- cv_df_ec$cntl_to_case_ratio[i]
  
  # has to be per period (not cumulative) for EmpiricalCalibration::computeCvBinomial()
  # i.e. test performed at 3 events then when 3 more events come in requires
  # GroupSizes = c(3, 3)
  gs_seq <- seq(look_i, max_n_i, by = look_i)
  if (max(gs_seq) != max_n_i) { # if doesn't go to max_n, add at end for last look
    gs_seq <- c(gs_seq, max_n_i) 
  }
  
  cv_df_ec$cv[i] <- 
    EmpiricalCalibration::computeCvBinomial(
      groupSizes = gs_seq,
      z = z_i,
      minimumEvents = min_event,
      alpha = alpha,
      sampleSize = 1e+07
    )
  
}
toc()

cv_df_ec



cv_df_ec %>%
  write_parquet(., sink = "out/cv_df_ec.parquet")



# ---- compare ----

inner_join(
  cv_df_ec,
  cv_df_seq,
  c("cntl_to_case_ratio", "max_n", "look_interval")
) %>%
  mutate(
    cv_diff = cv.x - cv.y,
    rel_cv_diff = cv_diff / ((cv.x - cv.y) / 2)
  )


# ---- maxsprt_calcs ----

max_sprt_stat <- function(c_n, n, z) {
  
  # the simple version of this equation is on page 70/71 of
  # Kulldorf et al. (2011) A Maximized Sequential Probability Ratio Test for
  # Drug and Vaccine Safety Surveillance. Sequential Analysis, 30(1): 58-78.
  RR0 <- 1 # null hypoth RR
  z_rr <- z / RR0
  if ((z_rr) * c_n / (n - c_n) <= 1) {
    max_llr <- 0
  } else if (c_n == n) {
    max_llr <- n * log(1 + z_rr)
  } else {
    max_llr <- 
      c_n * log(c_n / n) + 
      (n - c_n) * log((n - c_n) / n) - 
      c_n * log(1 / (z_rr + 1)) - 
      (n - c_n) * log(z_rr / (z_rr + 1))
  }
  return(max_llr)
  
}


rr_est <- function(c_n, n, z) {
  return(z * c_n / (n - c_n))
}

E_case <- function(c_n, n, z) {
  return(z * c_n / (n - c_n))
}

max_sprt_stat(4, 5, (1 + 10) / (4 + 12))
rr_est(4, 5, (1 + 10) / (4 + 12))
max_sprt_stat(7, 9, (11 + 6) / (16 + 6))
rr_est(7, 9, (11 + 6) / (16 + 6))




