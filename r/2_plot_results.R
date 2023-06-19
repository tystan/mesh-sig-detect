
# ---- lib ----


suppressPackageStartupMessages({
  library("readr")
  library("dplyr")
  library("tidyr")
  library("forcats")
  library("lubridate") # way to handle dates better than default R way
  library("ggplot2") 
  library("ggrepel") 
  library("knitr")
  library("gsDesign")
})







# ---- load_dat ----

sra_cum_bcpnn <-  read_parquet("out/sra_cum_bcpnn.parquet")


bcpnn_signif <-
  sra_cum_bcpnn %>%
  group_by(grps, dat_type, thresh) %>%
  arrange(dte) %>%
  dplyr::filter(reach_sig) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup()



# ---- plot ----





bcpnn_signif_plt <-
  bcpnn_signif %>%
  # keep only multiples of 0.01 (too many colours otherwise)
  dplyr::filter(abs(100 * thresh - floor(100 * thresh)) < 1e-6) %>%
  mutate(
    grps = gsub(" v ", "\nv\n", grps),
    grps = fct_inorder(grps)
  ) 

thresholds <- sort(unique(bcpnn_signif_plt[["thresh"]]))
length(thresholds)
thresh_scale <- rev(hcl.colors(length(thresholds) + 1, "Inferno"))[-1] 
# thresh_scale <- rev(hcl.colors(length(thresholds), "SunsetDark"))

bcpnn_signif_plt %>%
  arrange(grps, thresh) %>%
  ggplot(., aes(x = dte_reach_sig, y = thresh, col = est_name)) +
  # ggplot(., aes(x = dte_reach_sig, y = est_name, col = factor(thresh))) +
  geom_point() +
  geom_path(aes(group = est_name)) +
  # scale_colour_viridis_c(option = "B", direction = -1) +
  # scale_colour_manual(values = thresh_scale) +
  facet_wrap(~ grps, ncol = 1) +
  theme_bw()


sra_cum_bcpnn_plt <-
  sra_cum_bcpnn %>%
  # keep only multiples of 0.01 (too many colours otherwise)
  dplyr::filter(abs(100 * thresh - floor(100 * thresh)) < 1e-6) %>%
  mutate(
    grps = gsub(" v ", "\nv\n", grps),
    grps = fct_inorder(grps)
  )

thresholds <- sort(unique(sra_cum_bcpnn_plt[["thresh"]]))
length(thresholds)
thresh_scale <- rev(hcl.colors(length(thresholds), "SunsetDark"))


sra_cum_bcpnn_plt %>%
  ggplot(
    ., 
    aes(
      dte, 
      est, 
      ymax = ci_hi, 
      ymin = ci_lo, 
      col = factor(thresh), 
      fill = factor(thresh),
      group = factor(thresh),
      alpha = reach_sig,
      shape = reach_sig 
    )
  ) %+%
  geom_hline(aes(yintercept = 0), col = "grey50") %+% # null value
  geom_line() %+%
  geom_point() %+%
  geom_ribbon(alpha = 0.05, lty = 2)  %+%
  facet_wrap(~ grps, scales = "free_y", ncol = 1) %+%
  labs(
    subtitle = "(Calculations made on cumulative monthly data)",
    y = "IC statistic using BCPNN MCMC method (null value = 0)",
    x = "Date",
    col = "Pain score\nthreshold",
    fill = "Pain score\nthreshold"
  ) %+%
  scale_colour_manual(values = thresh_scale, aesthetics = c("colour", "fill")) %+%
  theme_bw() 





