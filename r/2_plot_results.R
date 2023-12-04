
# ---- libs ----


suppressPackageStartupMessages({
  library("readr")
  library("dplyr")
  library("tidyr")
  library("forcats")
  library("lubridate") # way to handle dates better than default R way
  library("stringr") 
  library("ggplot2") 
  library("ggthemes") 
  library("ggrepel") 
  library("knitr")
  library("gsDesign")
  library("arrow")
})



col_pal <- c("darkorange", "cyan4", "purple")



# ---- load_dat ----

# maxsprt_dat <- read_parquet("out/sra_cum_maxsprt.parquet")
# bcpnn_dat <- read_parquet("out/sra_cum_bcpnn_mc_adj.parquet")
# prr_dat <- read_parquet("out/sra_cum_prr_mc_adj.parquet")


sra_cum_prr_mc_adj <-  read_parquet("out/sra_cum_prr_mc_adj.parquet")
sra_cum_bcpnn_mc_adj <-  read_parquet("out/sra_cum_bcpnn_mc_adj.parquet")
sra_cum_maxsprt <-  read_parquet("out/sra_cum_maxsprt.parquet")


sra <-
  bind_rows(
    sra_cum_bcpnn_mc_adj %>% mutate(stat = "BCPNN (MCadj)"),
    sra_cum_prr_mc_adj %>% mutate(stat = "PRR (MCadj)"),
    sra_cum_maxsprt %>% mutate(stat = "maxSPRT")
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
        stat == "maxSPRT"     ~ maxllr, 
        stat == "PRR (MCadj)" ~ ci_lo, 
        stat == "BCPNN (MCadj)" ~ ci_lo
      ),
    test_thresh = 
      case_when(
        stat == "maxSPRT"     ~ cv, 
        stat == "PRR (MCadj)" ~ 1, 
        stat == "BCPNN (MCadj)" ~ 0
      ),
    rr_stat = 
      case_when(
        stat == "maxSPRT"     ~ rre, 
        stat == "PRR (MCadj)" ~ est, 
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
      select(stat, grps, thresh, dte, cv, val = ci_lo, reach_sig, dte_reach_sig),
    sra_cum_prr_mc_adj %>% 
      mutate(
        cv = 1, 
        stat = "RR (PRR, Lower 95% CI)"
      ) %>%
      select(stat, grps, thresh, dte, cv, val = ci_lo, reach_sig, dte_reach_sig),
    sra_cum_maxsprt %>%
      mutate(stat = "MaxSPRT (max LLR)") %>%
      select(stat, grps, thresh, dte, cv, val = maxllr, reach_sig, dte_reach_sig)
  ) 


sig_reach_dat <-
  plt_dat %>%
  arrange(stat, grps, thresh, dte) %>%
  group_by(stat, grps, thresh) %>%
  dplyr::filter(reach_sig == 1) %>%
  dplyr::filter(row_number() == 1) %>%
  select(stat, grps, thresh, dte_reached = dte) %>%
  # now create separation between reached CV values when it occurs
  group_by(stat, thresh, dte_reached) %>%
  mutate(rep_dte = 1:n()) %>%
  ungroup() %>%
  mutate(dte_reached = dte_reached + days(10 * (rep_dte - 1))) %>%
  select(-rep_dte)


plt_dat <-
  left_join(
    plt_dat,
    sig_reach_dat,
    c("stat", "grps", "thresh")
  ) # %>%
# dplyr::filter(dte_reach_sig != dte_reached)


plt_dat <-
  plt_dat %>%
  mutate(
    stat = fct_inorder(stat)
  )


# plt_dat %>%
#   dplyr::filter(thresh == "0.050", grepl("(c)", grps, fixed = TRUE))





# ---- example_data ----

cols_want <- c(
  "Quarter" = "mnth",
  "$t$" = "t",
  "$a_t$" = "nA",
  "$b_t$" = "nB",
  "$c_t$" = "nC",
  "$d_t$" = "nD"
)

sra_cum_prr_mc_adj %>%
  dplyr::filter(substr(grps, 1, 3) == "(a)", thresh == "0.050") %>%
  mutate(t = 1:n()) %>%
  select(all_of(cols_want)) %>%
  kable(.)

# ---- time_to_sig_plot1 ----



date_signif_dat <-
  sra %>%
  group_by(stat, grps, dat_type, thresh) %>%
  arrange(dte) %>%
  dplyr::filter(reach_sig) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(stat, grps, dat_type, thresh) 


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
  ggplot(., aes(x = dte_reach_sig, y = as.numeric(thresh), col = stat)) +
  geom_point() +
  geom_path(aes(group = stat)) +
  scale_colour_manual(values = col_pal) +
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
  facet_grid(stat ~ grps) +
  # facet_wrap( ~ grps, ncol = 1) +
  theme_bw() +
  theme(text = element_text(family = "serif")) +
  labs(
    x = expression("Date" ~ H[0] ~ "rejected (null hypothesis of no signal)"),
    y = "Threshold for P(topic = 'pain' | doc) to be classed a 'pain' adverse event",
    col = "Signal detection\nmethod"
  )

ggsave(
  filename = "fig/time_to_signal_method_facets.png", 
  dpi = 900, width = 9, height = 9
)



# ---- time_to_sig_plot2 ----


date_signif_dat <-
  sra %>%
  group_by(stat, grps, dat_type, thresh) %>%
  arrange(dte) %>%
  dplyr::filter(reach_sig) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(stat, grps, dat_type, thresh) 


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


signif_plt %>%
  arrange(grps, thresh) %>%
  mutate(dte_reach_sig = dte_reach_sig + (as.numeric(stat) - 1) * days(10)) %>%
  ggplot(., aes(x = dte_reach_sig, y = as.numeric(thresh), col = stat)) +
  geom_point(alpha = 0.4) +
  geom_path(aes(group = stat), alpha = 0.4) +
  scale_colour_manual(values = col_pal) +
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
  # facet_grid(stat ~ grps) +
  facet_wrap( ~ grps, ncol = 1) +
  theme_bw() +
  theme(text = element_text(family = "serif")) +
  labs(
    x = expression("Date" ~ H[0] ~ "rejected (null hypothesis of no signal)"),
    y = "Threshold for P(topic = 'pain' | doc) to be classed a 'pain' adverse event",
    col = "Signal detection\nmethod"
  )


ggsave(
  filename = "fig/time_to_signal_method_overlay.png", 
  dpi = 900, width = 5, height = 9
)



# ---- stat_over_time_plot1 ----



thresh_use <- 0.040
thresh_use_str <- sprintf("%0.3f", thresh_use)

sra_stat_plt <-
  sra %>%
  # keep only subset of thresholds (too many colours otherwise)
  # dplyr::filter(thresh %in% sprintf("%0.3f", seq(0.02, 0.08, by = 0.02))) %>%
  dplyr::filter(thresh %in% thresh_use_str) %>%
  ### only keep pelvic mesh as target vs whatever comparator
  dplyr::filter(grepl("^.*pelvic.* v ", grps)) %>%
  # dplyr::filter(grepl("^\\(a\\)", grps)) %>%
  mutate(
    grps = gsub(" v ", "\nv\n", grps),
    grps = gsub("\\([a-z]\\) ", "", grps),
    grps = gsub("_", " ", grps),
    grps = gsub("pelvic mesh", "Pelvic mesh", grps),
    grps = gsub("hernia mesh", "Hernia mesh", grps),
    grps = fct_inorder(grps),
    stat = fct_inorder(stat)
  ) 

thresholds <- sort(unique(sra_stat_plt[["thresh"]]))
# length(thresholds)
# thresh_scale <- rev(hcl.colors(length(thresholds), "SunsetDark"))
# thresh_scale <- rev(hcl.colors(length(thresholds) + 1, "Inferno"))[-1] 



sra_stat_plt <-
  sra_stat_plt %>%
  mutate(
    reach_sig_alpha = ifelse(reach_sig, 1, 0.8),
    `P(topic = 'pain') threshold` = thresh,
    `Test` = 
      case_when(
        # stat == "maxSPRT"       ~ paste("Maximised LLR > CV = ", sprintf("%1.3f", cv)), 
        stat == "maxSPRT"       ~ "Maximised LLR > CV", 
        stat == "BCPNN (MCadj)" ~ "IC lower 95% > 0",
        stat == "PRR (MCadj)"   ~ "PRR lower 95% > 1"
      ),
    `Test` = fct_inorder(`Test`),
    # truncate extreme values
    test_stat = if_else((test_stat > 30) & (stat == "maxSPRT"), 30, test_stat)
  ) 


sra_stat_plt %>%
  ggplot(
    ., 
    aes(
      x = dte, 
      y = test_stat, 
      col = stat, 
      group = interaction(stat, thresh)
      # alpha = reach_sig_alpha
    )
  ) %+%
  geom_hline(aes(yintercept = test_thresh), col = "black") %+% # null value
  geom_vline(aes(xintercept = dte_reach_sig, col = stat), alpha = 0.5) %+% # sig first reached
  geom_line() %+%
  geom_point() %+%
  # geom_ribbon(alpha = 0.05, lty = 2)  %+%
  # facet_wrap(~ grps, scales = "free_y", ncol = 1) %+%
  facet_grid(
    `Test` ~ grps, 
    scales = "free_y", 
    labeller = labeller(Test = function(x) paste0("Test: ", x))
  ) %+%
  labs(
    # subtitle = "Pelvic mesh v hernia mesh",
    subtitle = paste0("P(topic = 'pain') threshold = ", thresh_use_str),
    y = "Test statistic",
    x = "Date (quarterly data accumulation)",
    col = "Signal detection\nmethod"
  ) %+%
  scale_x_continuous(
    breaks = as_date(paste0(seq(2013, 2017, by = 1), "-01-01")), 
    labels = function(x) year(x)
  ) %+%
  # datetime_scale(trans = "date", breaks = as_date(paste0(seq(2013, 2017, by = 1), "-01-01"))) %+%
  scale_colour_manual(values = col_pal)  %+%
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
  theme_bw() %+%
  theme(text = element_text(family = "serif")) 


ggsave(
  filename = "fig/multi-grps_multi-test_sig_detect_over_time_thresh-0.04.png", 
  dpi = 900, width = 10, height = 8
)


# ---- stat_over_time_plot2 ----


thresh_lablr <- function(string) paste0("P(topic = 'pain' | doc)\nthreshold = ", string)



thresh_use <- seq(0.02, 0.08, by = 0.02)
thresh_use_str <- sprintf("%0.3f", thresh_use)

sra_stat_plt <-
  sra %>%
  # keep only subset of thresholds (too many colours otherwise)
  # dplyr::filter(thresh %in% sprintf("%0.3f", seq(0.02, 0.08, by = 0.02))) %>%
  dplyr::filter(thresh %in% thresh_use_str) %>%
  ### only keep pelvic mesh as target vs whatever comparator
  # dplyr::filter(grepl("^.*pelvic.* v ", grps)) %>%
  dplyr::filter(grepl("^\\(a\\)", grps)) %>%
  mutate(
    grps = gsub(" v ", "\nv\n", grps),
    grps = gsub("\\([a-z]\\) ", "", grps),
    grps = gsub("_", " ", grps),
    grps = gsub("pelvic mesh", "Pelvic mesh", grps),
    grps = gsub("hernia mesh", "Hernia mesh", grps),
    grps = fct_inorder(grps),
    stat = fct_inorder(stat)
  ) 

thresholds <- sort(unique(sra_stat_plt[["thresh"]]))
# length(thresholds)
# thresh_scale <- rev(hcl.colors(length(thresholds), "SunsetDark"))
# thresh_scale <- rev(hcl.colors(length(thresholds) + 1, "Inferno"))[-1] 



sra_stat_plt <-
  sra_stat_plt %>%
  mutate(
    reach_sig_alpha = ifelse(reach_sig, 1, 0.8),
    `P(topic = 'pain') threshold` = thresh,
    `Test` = 
      case_when(
        # stat == "maxSPRT"       ~ paste("Maximised LLR > CV = ", sprintf("%1.3f", cv)), 
        stat == "maxSPRT"       ~ "Maximised LLR > CV", 
        stat == "BCPNN (MCadj)" ~ "IC lower 95% > 0",
        stat == "PRR (MCadj)"   ~ "PRR lower 95% > 1"
      ),
    `Test` = fct_inorder(`Test`),
    # truncate extreme values
    test_stat = if_else((test_stat > 30) & (stat == "maxSPRT"), 30, test_stat)
  ) 

levels(sra_stat_plt$stat) <- gsub(" (MCadj)", "", levels(sra_stat_plt$stat), fixed = TRUE)


sra_stat_plt %>%
  mutate(
    stat2 = 
      paste0(
        ifelse(
          stat == "BCPNN", 
          "2^IC =\nP(Pain AE & Pelvic)/{P(Pain AE)P(Pelvic)}", 
          "RR"
        ), 
        "\n[", stat, "]"
      ),
    stat2 = fct_inorder(stat2),
    reach_sig_alpha = ifelse(reach_sig, 1, 0.8)
  ) %>%
  ggplot(
    ., 
    aes(
      x = dte, 
      y = rr_stat, 
      col = stat2, 
      group = stat2
      # group = interaction(stat2, thresh)
      # alpha = reach_sig_alpha
      # ymin = ci_lo,
      # ymax = ci_hi
    )
  ) %+%
  geom_hline(aes(yintercept = 1), col = "black") %+% # null value
  geom_vline(aes(xintercept = dte_reach_sig), alpha = 0.5) %+% # sig first reached
  geom_line() %+%
  geom_point() %+%
  # geom_ribbon(alpha = 0.05)  %+%
  # facet_wrap(~ thresh, nrow = 1, labeller = as_labeller(thresh_lablr)) %+%
  facet_grid(
    stat ~ thresh,
    labeller = labeller(thresh = thresh_lablr)
  ) %+%
  labs(
    subtitle = "Pelvic mesh v hernia mesh",
    y = "Reporting ratio estimate",
    x = "Date (quarterly data accumulation)",
    col = "Reporting ratio estimate\n[signal detection method]"
  ) %+%
  scale_y_continuous(trans = "log2") %+%
  scale_colour_manual(values = col_pal)  %+%
  theme_bw() %+%
  theme(
    text = element_text(family = "serif"),
    legend.key.height = unit(2.5, units = "line")
  ) 


ggsave(
  filename = "fig/pelvic_v_hernia_rr_est_over_time.png", 
  dpi = 900, width = 12, height = 8
)




# ---- multi-grps_multi-test_signal-over-time1 ----




plt_dat %>%
  dplyr::filter(thresh == "0.030", !grepl("(e)", grps, fixed = TRUE)) %>%
  ggplot(., aes(x = dte, y = val, col = grps )) +
  geom_hline(aes(yintercept = cv), alpha = 0.5) +
  geom_vline(aes(xintercept = dte_reached, col = grps), alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_point() +
  # facet_wrap(thresh ~ stat, ncol = 1, scales = "free_y") +
  facet_grid(
    stat ~ thresh, 
    scales = "free_y", 
    labeller = labeller(thresh = function(x) paste0("Threshold: ", x))
  ) +
  scale_colour_tableau() +
  theme_bw() +
  labs(
    x = "Quarter",
    y = "Statistic",
    col = "Comparison"
  )

ggsave(
  filename = "fig/multi-grps_multi-test_signal-over-time_thresh-0.03.png", 
  dpi = 900, width = 12, height = 8
)



# ---- multi-grps_multi-test_signal-over-time2 ----


plt_dat %>%
  dplyr::filter(thresh == "0.070", !grepl("(e)", grps, fixed = TRUE)) %>%
  ggplot(., aes(x = dte, y = val, col = grps )) +
  geom_hline(aes(yintercept = cv), alpha = 0.5) +
  geom_vline(aes(xintercept = dte_reached, col = grps), alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_point() +
  # facet_wrap(thresh ~ stat, ncol = 1, scales = "free_y") +
  facet_grid(
    stat ~ thresh, 
    scales = "free_y", 
    labeller = labeller(thresh = function(x) paste0("Threshold: ", x))
  ) +
  scale_colour_tableau() +
  theme_bw() +
  labs(
    x = "Quarter",
    y = "Statistic",
    col = "Comparison"
  )

ggsave(
  filename = "fig/multi-grps_multi-test_signal-over-time_thresh-0.07.png", 
  dpi = 900, width = 12, height = 8
)


# ---- multi-grps_multi-test_signal-over-time3 ----


plt_dat %>%
  dplyr::filter(thresh %in% c("0.040", "0.050", "0.060"), !grepl("(e)", grps, fixed = TRUE)) %>%
  ggplot(., aes(x = dte, y = val, col = grps )) +
  geom_hline(aes(yintercept = cv), alpha = 0.5) +
  geom_vline(aes(xintercept = dte_reached, col = grps), alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_point() +
  # facet_wrap(thresh ~ stat, ncol = 1, scales = "free_y") +
  facet_grid(
    stat ~ thresh, 
    scales = "free_y", 
    labeller = labeller(thresh = function(x) paste0("Threshold: ", x))
  ) +
  scale_colour_tableau() +
  theme_bw() +
  labs(
    x = "Quarter",
    y = "Statistic",
    col = "Comparison"
  )




ggsave(
  filename = "fig/multi-grps_multi-test_signal-over-time_thresh-range-0.04-0.06.png", 
  dpi = 900, width = 12, height = 8
)






