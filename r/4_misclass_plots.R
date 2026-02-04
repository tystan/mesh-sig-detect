
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


col_pal <- c("darkorange", "cyan4", "purple")


## ---- plot_sim ----


sra_cum_prr_sim_mc_adj <-  read_parquet("out/sra_cum_prr_sim_mc_adj.parquet")
sra_cum_bcpnn_sim_mc_adj <-  read_parquet("out/sra_cum_bcpnn_sim_mc_adj.parquet")
sra_cum_maxsprt_sim <-  read_parquet("out/sra_cum_maxsprt_sim.parquet")


sra <-
  bind_rows(
    sra_cum_bcpnn_sim_mc_adj %>% mutate(stat = "BCPNN (MCadj)"),
    sra_cum_prr_sim_mc_adj %>% mutate(stat = "PRR (MCadj)"),
    sra_cum_maxsprt_sim %>% mutate(stat = "maxSPRT")
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
    sra_cum_bcpnn_sim_mc_adj %>% 
      mutate(
        cv = 0, 
        stat = "IC (BCPNN, Lower 95% CI)"
      ) %>%
      select(stat, grps, thresh, sim_i, dte, cv, val = ci_lo, reach_sig, dte_reach_sig),
    sra_cum_prr_sim_mc_adj %>%
      mutate(
        cv = 1,
        stat = "RR (PRR, Lower 95% CI)"
      ) %>%
      select(stat, grps, thresh, sim_i, dte, cv, val = ci_lo, reach_sig, dte_reach_sig),
    sra_cum_maxsprt_sim %>%
      mutate(stat = "MaxSPRT (max LLR)") %>%
      select(stat, grps, thresh, sim_i, dte, cv, val = maxllr, reach_sig, dte_reach_sig)
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
  # geom_vline(xintercept = quantile(signif_plt$dte_reach_sig, c(0.01, 0.5, 0.99), type = 3)) +
  facet_grid(stat ~ grps) +
  # facet_wrap( ~ grps, ncol = 1) +
  theme_bw() +
  theme(text = element_text(family = "serif")) +
  labs(
    x = expression("Date" ~ H[0] ~ "rejected (null hypothesis of no signal)"),
    y = "Number of simulations",
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






signif_heat_plt <-
  signif_plt %>%
  mutate(
    dte_reach_sig_cat =
      case_when(
        dte_reach_sig < "2014-01-01" ~ "<=2013(Q4)",
        dte_reach_sig < "2014-10-01" ~ "2014(Q1-Q3)",
        dte_reach_sig == "2014-10-01" ~ "2014(Q4)",
        dte_reach_sig < "2016-01-01" ~ "2015(Q1-Q4)",
        dte_reach_sig >= "2016-01-01" ~ ">=2016",
        TRUE ~ "ERROR!!!" # else have error
      ),
    dte_reach_sig_cat = 
      factor(
        dte_reach_sig_cat,
        levels = c("<=2013(Q4)", "2014(Q1-Q3)", "2014(Q4)", "2015(Q1-Q4)", ">=2016")
      )
  )

signif_heat_plt %>%
  arrange(grps, thresh) %>%
  ggplot(., aes(x = dte_reach_sig_cat, fill = stat)) +
  geom_bar(alpha = 0.8, col = NA) +
  # geom_path(aes(group = interaction(stat, sim_i))) +
  scale_fill_manual(values = col_pal) +
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
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


# signif_heat_plt <-
#   signif_heat_plt %>%
#   group_by(stat, grps, dat_type, thresh, dte_reach_sig_cat) %>% 
#   summarise(
#     sim_cnt = n(),
#     .groups = "drop"
#   ) %>%
#   group_by(stat, grps, dat_type, thresh) %>% 
#   mutate(sim_pct = 100 * sim_cnt / sum(sim_cnt)) %>% 
#   ungroup() %>% 
#   mutate(sim_pct_lbl = sprintf("%2.1f%%", sim_pct)) 
# 
# 
# 
# 
# 
# signif_heat_plt %>%
#   ggplot(., aes(x = dte_reach_sig_cat, y = stat, fill = sim_pct)) +
#   geom_tile() +
#   geom_text(aes(label = sim_pct_lbl)) +
#   # geom_path(aes(group = interaction(stat, sim_i))) +
#   scale_fill_viridis_c() +
#   # scale_colour_tableau(palette = "Color Blind", direction = -1) +
#   facet_grid(thresh ~ grps) +
#   # facet_wrap( ~ grps, ncol = 1) +
#   theme_bw() +
#   theme(text = element_text(family = "serif")) +
#   labs(
#     x = expression("Date" ~ H[0] ~ "rejected (null hypothesis of no signal)"),
#     # y = "Threshold for P(topic = 'pain' | doc) to be classed a 'pain' adverse event",
#     y = "Signal detection method",
#     fill = "Percent of\nmisclasificaiton\nsimulations"
#   )
# 
# 
# 
# ggsave(
#   filename = "fig/sim_data_categoryt_time_to_signal_method_facets.pdf",
#   device = cairo_pdf, # embed fonts
#   width = 9, height = 5, units = "in"
# )





# ---- stat_over_time_plot1 ----



thresh_use <- 0.050
thresh_use_str <- sprintf("%0.3f", thresh_use)

sra_stat_plt <-
  sra %>%
  # keep only subset of thresholds (too many colours otherwise)
  # dplyr::filter(thresh %in% sprintf("%0.3f", seq(0.02, 0.08, by = 0.02))) %>%
  dplyr::filter(thresh %in% thresh_use_str) %>%
  ### only keep pelvic mesh as target vs whatever comparator
  # dplyr::filter(grepl("^.*pelvic.* v ", grps)) %>%
  dplyr::filter(grepl("^\\([a-d]\\)", grps)) %>%
  mutate(
    grps = gsub(" v ", "\nv\n", grps),
    grps = gsub("\\([a-z]\\) ", "", grps),
    grps = gsub("_", " ", grps),
    grps = gsub("pelvic mesh", "Pelvic mesh", grps),
    grps = gsub("hernia mesh", "Hernia mesh", grps),
    grps = gsub("other mesh", "Other mesh", grps),
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
        # stat == "maxSPRT"       ~ "Maximised LLR > CV", 
        stat == "BCPNN (MCadj)" ~ "IC lower 95% > 0",
        # stat == "PRR (MCadj)"   ~ "PRR lower 95% > 1"
      ),
    `Test` = fct_inorder(`Test`),
    # truncate extreme values
    test_stat = if_else((test_stat > 30) & (stat == "maxSPRT"), 30, test_stat),
    test_stat = if_else((test_stat < -3) & (stat == "BCPNN (MCadj)"), -3, test_stat),
    dte_reach_sig = if_else(dte_reach_sig > "2017-12-01", as_date(NA), dte_reach_sig),
    dte_reach_insta_sig =
      if_else(
        is.na(dte_reach_sig) | (dte_reach_sig != dte), 
        as_date(NA), 
        dte_reach_sig
      )
  ) 


print(sra_stat_plt, n = 25)

sra_stat_plt %>%
  ggplot(
    ., 
    aes(
      x = dte, 
      y = test_stat, 
      col = stat, 
      group = interaction(stat, sim_i)
      # alpha = reach_sig_alpha
    )
  ) +
  geom_hline(aes(yintercept = test_thresh), col = "black") + # null value
  # geom_vline(aes(xintercept = dte_reach_sig, col = stat), alpha = 0.5) + # sig first reached
  geom_line(alpha = 0.05)+
  geom_jitter(
    data = filter(sra_stat_plt, !is.na(dte_reach_insta_sig)), 
    aes(x = dte_reach_insta_sig),
    height = 0, width = 20, alpha = 0.5
  ) +
  facet_grid(
    `Test` ~ grps, 
    scales = "free_y", 
    labeller = labeller(Test = function(x) paste0("Test: ", x))
  ) +
  labs(
    # subtitle = "Pelvic mesh v hernia mesh",
    subtitle = paste0("P(topic = 'pain') threshold = ", thresh_use_str),
    y = "Test statistic",
    x = "Date (quarterly data accumulation)",
    col = "Signal detection\nmethod"
  ) +
  scale_x_continuous(
    breaks = as_date(paste0(seq(2013, 2017, by = 1), "-01-01")), 
    labels = function(x) year(x)
  ) +
  # datetime_scale(trans = "date", breaks = as_date(paste0(seq(2013, 2017, by = 1), "-01-01"))) +
  scale_colour_manual(values = col_pal)  +
  # scale_colour_tableau(palette = "Color Blind", direction = -1) +
  theme_bw() +
  theme(text = element_text(family = "serif")) 


ggsave(
  filename = "fig/sim_multi-grps_multi-test_sig_detect_over_time_thresh-0.05.png", 
  dpi = 900, width = 10, height = 8
)

























