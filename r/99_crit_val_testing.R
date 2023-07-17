
# ---- libs ----


library("dplyr")
library("tidyr")
library("tictoc")
library("foreach")
library("EmpiricalCalibration")
library("Sequential")
library("arrow")
library("ggplot2")


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

# ---- example_cv ----

look_i <- 10
max_n_i <- 200
z_i <- 2


# has to be per period (not cumulative) for Sequential::CV.Binomial()
# (and EmpiricalCalibration::computeCvBinomial())
# i.e. test performed at 3 events then when 4 more events come in requires
# GroupSizes = c(3, 4)
### NOT: GroupSizes = c(3, 7)
gs_seq <- rep(look_i, floor(max_n_i / look_i))
if (sum(gs_seq) != max_n_i) { # if doesn't go to max_n, add at end for last look
  gs_seq <- c(gs_seq, max_n_i - sum(gs_seq)) 
}
gs_seq

# both of the below take ~ 6 sec
tic()
Sequential::CV.Binomial(
  N = max_n_i,
  alpha = alpha,
  M = min_event,
  z = z_i, 
  GroupSizes = gs_seq
)
toc()


tic()
EmpiricalCalibration::computeCvBinomial(
  groupSizes = gs_seq,
  z = z_i,
  minimumEvents = 1,
  alpha = 0.05, # does two-tailed by default? (no)
  sampleSize = 1e+06
)
toc()

# ---- example_cv2 ----

look_i <- 10
max_n_i <- 4000
z_i <- 2


# has to be per period (not cumulative) for Sequential::CV.Binomial()
# (and EmpiricalCalibration::computeCvBinomial())
# i.e. test performed at 3 events then when 4 more events come in requires
# GroupSizes = c(3, 4)
### NOT: GroupSizes = c(3, 7)
gs_seq <- rep(look_i, floor(max_n_i / look_i))
if (sum(gs_seq) != max_n_i) { # if doesn't go to max_n, add at end for last look
  gs_seq <- c(gs_seq, max_n_i - sum(gs_seq)) 
}
gs_seq

# takes ~ 2 min but Sequential::CV.Binomial()
# warns against using for sum(groupSizes) > 500
tic()
EmpiricalCalibration::computeCvBinomial(
  groupSizes = gs_seq,
  z = z_i,
  minimumEvents = 1,
  alpha = 0.05, # does two-tailed by default?
  sampleSize = 1e+06
)
toc()
### 103.83 sec elapsed


# ---- seq_package ----


### testing Sequential::CV.Binomial()

cv_df_seq <- cv_df_0
cv_df_seq$cv <- 0
cv_df_seq$pack <- "Sequential"
cv_df_seq



tic() # takes 5 min 
for (i in 1:nrow(cv_df_seq)) {
  
  look_i <- cv_df_seq$look_interval[i] 
  max_n_i <- cv_df_seq$max_n[i]
  z_i <- cv_df_seq$cntl_to_case_ratio[i]
  
  # has to be per period (not cumulative) for Sequential::CV.Binomial()
  # i.e. test performed at 3 events then when 3 more events come in requires
  # GroupSizes = c(3, 3)
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


tic() # takes 4 hours using sampleSize = 1e+07
for (i in 1:nrow(cv_df_ec)) {
  
  look_i <- cv_df_ec$look_interval[i] 
  max_n_i <- cv_df_ec$max_n[i]
  z_i <- cv_df_ec$cntl_to_case_ratio[i]
  
  # has to be per period (not cumulative) for EmpiricalCalibration::computeCvBinomial()
  # i.e. test performed at 3 events then when 3 more events come in requires
  # GroupSizes = c(3, 3)
  gs_seq <- rep(look_i, floor(max_n_i / look_i))
  if (sum(gs_seq) != max_n_i) { # if doesn't go to max_n, add at end for last look
    gs_seq <- c(gs_seq, max_n_i - sum(gs_seq)) 
  }
  
  cv_df_ec$cv[i] <- 
    EmpiricalCalibration::computeCvBinomial(
      groupSizes = gs_seq,
      z = z_i,
      minimumEvents = min_event,
      alpha = alpha,
      sampleSize = 1e+06
    )
  
}
toc()

cv_df_ec



cv_df_ec %>%
  write_parquet(., sink = "out/cv_df_ec.parquet")



# ---- compare ----

cv_df <-
  bind_rows(
    read_parquet("out/cv_df_ec.parquet"),
    read_parquet("out/cv_df_seq.parquet")
  ) 



cv_df %>%
  mutate(
    cv_txt = sprintf("%1.1f", cv)
  ) %>%
  ggplot(., aes(x = factor(look_interval), y = cntl_to_case_ratio)) +
  geom_tile(aes(fill = cv)) +
  geom_text(aes(label = cv_txt)) +
  facet_wrap(~ pack) +
  scale_fill_viridis_b() +
  labs(
    x = expression(paste(
      "Group sequential testing: expected number of events between tests: ", 
      (a[t] - a[t-1]) + (c[t] - c[t-1]), 
      ")"
    )),
    y = expression(paste(
      "z (control:case ratio = ", 
      (c[t] + d[t]) / (a[t] + b[t])
    )),
    fill = "maxSPRT\ncritical\nvalue",
    label = "maxSPRT\ncritical\nvalue",
  ) +
  theme_bw()


cv_df_w <-
  cv_df %>%
  pivot_wider(
    .,
    id_cols = all_of(c("cntl_to_case_ratio", "max_n", "look_interval")),
    names_from = "pack",
    values_from = "cv"
  ) %>%
  mutate(
    cv_diff = EmpiricalCalibration - Sequential,
    rel_cv_diff = cv_diff / ((EmpiricalCalibration + Sequential) / 2)
  )


cv_df_w %>%
  select(
    `max n` = max_n,
    `z=E[cntl:case]` = cntl_to_case_ratio,
    `events per look` = look_interval,
    empcalib = EmpiricalCalibration,
    seq = Sequential,
    diff = cv_diff,
    `std diff` = rel_cv_diff
  ) %>%
  knitr::kable(., digits = 2)


# cv_df_w %>%
#   mutate(
#     rel_cv_diff_txt = sprintf("%0.3f", rel_cv_diff)
#   ) %>%
#   ggplot(., aes(x = factor(look_interval), y = cntl_to_case_ratio)) +
#   geom_tile(aes(fill = rel_cv_diff)) +
#   geom_text(aes(label = rel_cv_diff_txt)) +
#   scale_fill_viridis_b(option = "B") +
#   theme_bw()




