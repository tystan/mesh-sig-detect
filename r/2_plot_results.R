
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

sra_cum_bcpnn <-  read_parquet("out/sra_cum_bcpnn.parquet")
sra_cum_bcpnn_mc_adj <-  read_parquet("out/sra_cum_bcpnn_mc_adj.parquet")
sra_cum_maxsprt <-  read_parquet("out/sra_cum_maxsprt.parquet")


sra <-
  bind_rows(
    sra_cum_bcpnn %>% mutate(stat = "BCPNN"),
    sra_cum_bcpnn_mc_adj %>% mutate(stat = "BCPNN (MCadj)"),
    sra_cum_maxsprt %>% mutate(stat = "maxSPRT")
  ) %>%
  select(stat, everything())

sra <-
  sra %>%
  mutate(
    test_stat = if_else(stat == "maxSPRT", maxllr, ci_lo),
    test_thresh = if_else(stat == "maxSPRT", cv, 0),
    rr_stat = if_else(stat == "maxSPRT", rre, 2 ^ est)
  )




thresholds <- sort(unique(sra[["thresh"]]))
length(thresholds)


# ---- time_to_sig_plot ----



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
  ) 

levels(signif_plt$grps)

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

signif_plt %>%
  arrange(grps, thresh) %>%
  ggplot(., aes(x = dte_reach_sig, y = as.numeric(thresh), col = stat)) +
  geom_point() +
  geom_path(aes(group = stat)) +
  scale_colour_manual(values = col_pal) +
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
  # facet_grid(stat ~ grps) +
  facet_wrap( ~ grps, ncol = 1) +
  theme_bw() +
  theme(text = element_text(family = "serif")) +
  labs(
    x = "Date H0 of no signal rejected",
    y = "Threshold for P(topic = 'pain' | doc) to be classed a 'pain' adverse event",
    col = "Signal detection\nmethod"
  )


ggsave(
  filename = "fig/time_to_signal_method_overlay.png", 
  dpi = 900, width = 5, height = 9
)



# ---- stat_over_time_plot ----



sra_stat_plt <-
  sra %>%
  # keep only subset of thresholds (too many colours otherwise)
  dplyr::filter(thresh %in% sprintf("%0.3f", seq(0.02, 0.08, by = 0.02))) %>%
  ### only keep pelvic mesh as target vs whatever comparator
  # dplyr::filter(grepl("^.*pelvic.* v ", grps)) %>%
  dplyr::filter(grepl("^\\(a\\)", grps)) %>%
  mutate(
    grps = gsub(" v ", "\nv\n", grps),
    grps = gsub("\\([a-z]\\) ", "", grps),
    grps = gsub("_", " ", grps),
    grps = fct_inorder(grps)
  ) 

thresholds <- sort(unique(sra_stat_plt[["thresh"]]))
length(thresholds)
thresh_scale <- rev(hcl.colors(length(thresholds), "SunsetDark"))
# thresh_scale <- rev(hcl.colors(length(thresholds) + 1, "Inferno"))[-1] 

sra_stat_plt %>%
  mutate(
    reach_sig_alpha = ifelse(reach_sig, 1, 0.8),
    `P(topic = 'pain') threshold` = thresh,
    `Test` = if_else(stat == "maxSPRT", "Maximised LLR > CV", "IC lower 95% > 0")
  ) %>%
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
    `Test` ~ `P(topic = 'pain') threshold`, 
    scales = "free_y", 
    labeller = label_both
  ) %+%
  labs(
    subtitle = "Pelvic mesh v hernia mesh",
    y = "Test statistic",
    x = "Date (quarterly data accumulation)",
    col = "Signal detection\nmethod"
  ) %+%
  # scale_y_continuous(limits = c(NA, 6)) %+%
  scale_colour_manual(values = col_pal)  %+%
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
  theme_bw() %+%
  theme(text = element_text(family = "serif")) 


ggsave(
  filename = "fig/pelvic_v_hernia_sig_detect_over_time.png", 
  dpi = 900, width = 10, height = 8
)



sra_stat_plt %>%
  dplyr::filter(stat != "BCPNN") %>%
  mutate(
    stat2 = 
      paste0(
        ifelse(stat == "maxSPRT", "RR", "2^IC =\nP(Pain AE)P(Pelvic)/P(Pain AE & Pelvic)"), 
        "\n(", stat, ")"
      ),
    reach_sig_alpha = ifelse(reach_sig, 1, 0.8),
    `P(topic = 'pain') threshold` = thresh,
    `Test` = if_else(stat == "maxSPRT", "Maximised LLR > CV", "IC lower 95% > 0"),
    `Statistic calculation method` = stat
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
  facet_wrap(~ `P(topic = 'pain') threshold`, nrow = 1, labeller = label_both) %+%
  # facet_grid(
  #   `Statistic calculation method` ~ `P(topic = 'pain') threshold`, 
  #   labeller = label_both
  # ) %+%
  labs(
    subtitle = "Pelvic mesh v hernia mesh",
    y = "Reporting ratio estimate",
    x = "Date (quarterly data accumulation)",
    col = "Reporting ratio estimate\n(signal detection method)"
  ) %+%
  scale_y_continuous(trans = "log2") %+%
  scale_colour_manual(values = col_pal[-2])  %+%
  theme_bw() %+%
  theme(text = element_text(family = "serif")) 


ggsave(
  filename = "fig/pelvic_v_hernia_rr_est_over_time.png", 
  dpi = 900, width = 12, height = 5
)







