
# ---- lib ----


suppressPackageStartupMessages({
  library("readr")
  library("dplyr")
  library("tidyr")
  library("lubridate") # way to handle dates better than default R way
  library("ggplot2") 
  library("ggrepel") 
  library("knitr")
  library("gsDesign")
})

# set up some parallel computation for speed
library("foreach") # could use purrr equally
library("doParallel")
no_cores <- detectCores() - 1 # Calculate the number of cores (leave one free)
cl <- makeCluster(no_cores) # Create clusters
registerDoParallel(cl) # and register

# NOTE: need to run first (only once, assumes devtools installed):
# devtools::install_github("tystan/pharmsignal") 
library("pharmsignal") # signal detection algs

# here are the functions written for these analyses
# they will be shown in the *Appendix A*
source("r/_funcs.R")




# ---- consts ----

# arbitrarily, let's go with minimum cell count of 3 (should be discussed!)
arbitrary_cell_min <- 1

# these are the thresholds for pain_topic to be pain == TRUE
thresholds <- c(0.010, 0.025, 0.05, 0.075, 0.100) #, 0.150) 
thresholds <- seq(0.010, 0.100, by = 0.005) #, 0.150) 


# ---- load_dat ----


clean_data_cols <-
  cols(
    Report_ID = col_double(),
    Date = col_date(format = ""),
    pain_word = col_logical(),
    pain_topic = col_double(),
    type = col_character()
  )

clean_data <- read_csv("data/clean_data.csv", col_types = clean_data_cols)

clean_data %>%
  dplyr::filter(type == "other_mesh") %>%
  select(Report_ID)



cat("First 10 rows of raw data\n")
clean_data %>%
  arrange(Date) %>%
  dplyr::filter(row_number() < 11) %>%
  kable(.)

clean_data <-
  clean_data %>%
  dplyr::filter(
    type %in% c("pelvic_mesh", "hernia_mesh")
  )


clean_data %>%
  with(., table(type, pain_word)) %>%
  knitr::kable(.)

clean_data %>%
  with(., table(type, pain_topic >= 0.05)) %>%
  knitr::kable(.)


### Example calculations
# 4 / 42
# 1186 /   12752
# (4 / 42) /(1186 /   12752)
# ror_signal(70, 32, 1186, 12752, alpha = 0.05)



# These are the device groups and subgroups.
clean_data %>% 
  group_by(type) %>% 
  summarise(count = n()) %>%
  kable(.)

cat("\n\n## Histogram of `pain_word` (boolean) v `pain_topic` (score)")

clean_data %>%
  ggplot(., aes(pain_topic, col = pain_word, fill = pain_word)) +
  geom_histogram() +
  theme_bw()

clean_data %>% 
  dplyr::filter(
    type %in% c("pelvic_mesh", "hernia_mesh", "other_mesh", "other_device")
  ) %>%
  ggplot(., aes(pain_topic, col = pain_word, fill = pain_word)) +
  geom_histogram() +
  facet_wrap(~ type, scales = "free_y") +
  theme_bw()


# ---- wrangle ----

# Example 1:
# Use pelvic mesh as group 1 and all other mesh devices (including hernia) as group 2. 
# The value of interest is the pain topic, being above the threshold of 0.05. 
# (i.e. 5% of the document contains words from the pain topic)
# You can adjust the topic threshold if you want to balance the groups more. 
# A higher topic_threshold will look for documents that discuss "pain" more, and 
# hence find less pain documents.

# Example 2:
# group 1 is pelvic mesh devices and the comparator is hernia mesh devices
# The value of interest is pain_word (i.e. the presence of the word "pain" 
# in the event description.)

comparator_lst <-
  list(
    "hernia_mesh"
    # "other_mesh",
    # c("hernia_mesh", "other_mesh")
    # c("other_mesh", "hernia_mesh", "other_device")
  )

names(comparator_lst) <-
  paste0(
    1:length(comparator_lst),
    ". Pelvic mesh compared to ",
    c(
      "hernia mesh"
      # "all other mesh devices (EXCLUDING hernia)",
      # "all other mesh devices (including hernia)"
      # "all other devices\n(INCLUDING non-mesh devices)"
    ) #,
    #" for pain score >= threshold"
  )


i <- 1

cat(
  "\n\n\\newpage\n\n## Analysis ", 
  names(comparator_lst)[i], 
  " for pain score >= threshold\n\n\n", 
  sep = ""
)

dat1 <-
  foreach(th_i = thresholds, .combine = bind_rows, .packages = "dplyr") %do% {
    get_signal_dat(
      g1 = "pelvic_mesh",
      g2 = comparator_lst[[i]],
      pain_type = "pain_topic", 
      thresh = th_i,
      cell_min = 1,
      verbose = FALSE
    ) %>%
      bind_cols(., thresh = th_i)
  }


dat1 %>%
  dplyr::filter(thresh == 0.05) %>%
  kable()


# so this is multiple comparisons central but let's create disproportionality stats
# whenever a new report enters the data
n_reports <- nrow(dat1)
# takes ~1 sec using i5-8400
system.time({
  da_stats <-
    foreach(i = 1:n_reports, .combine = bind_rows, .packages = "dplyr") %dopar% {
      with(dat1, pretty_da(mnth[i], thresh[i], nA[i], nB[i], nC[i], nD[i]))
    }
})

dat1 <-
  dat1 %>%
  inner_join(
    .,
    da_stats,
    c("mnth", "thresh")
  )

n_post_join <- nrow(dat1)
if (n_post_join != n_reports) {
  stop("join not 1-1")
}



# first signif
bcpnn_signif <-
  dat1 %>%
  group_by(thresh) %>%
  dplyr::filter(ci_lo > 0) %>%
  arrange(mnth) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(dte = as_date(paste0(mnth, "-01")))

length(thresholds)
thresh_scale <- rev(hcl.colors(length(thresholds) + 1, "Inferno"))[-1] # SunsetDark



bcpnn_signif %>%
  ggplot(., aes(x = dte, y = est_name, col = factor(thresh))) +
  geom_point() +
  geom_label_repel(aes(label = factor(thresh))) +
  scale_colour_manual(values = thresh_scale) +
  theme_bw()





# ---- multcomp ----

dat11 <-
  foreach(i = 1:length(thresholds), .combine = bind_rows, .packages = "dplyr") %do% {
    
    dat1 <-
      get_signal_dat(
        g1 = "pelvic_mesh",
        g2 = "hernia_mesh",
        pain_type = "pain_topic", 
        thresh = thresholds[i],
        cell_min = 1,
        verbose = FALSE
      ) %>%
      bind_cols(., thresh = thresholds[i])
    
    
    
    
    # so this is multiple comparisons central but let's create disproportionality stats
    # whenever a new report enters the data
    n_reports <- nrow(dat1)
    
    information_fracs <-  1:n_reports / n_reports
    # spend_obj <- sfLDPocock(alpha = 0.025, t = information_fracs, param = NULL)
    # spend_obj <- sfLDOF(alpha = 0.025, t = information_fracs, param = NULL)
    spend_obj <- sfExponential(alpha = 0.05, t = information_fracs, param = 0.5)
    
    # plot(1:n_reports, spend_obj$spend, main = "alpha spending func", xlab = "look")
    
    
    # takes ~1 sec using i5-8400
    system.time({
      da_stats <-
        foreach(i = 1:n_reports, .combine = bind_rows, .packages = "dplyr") %dopar% {
          with(
            dat1, 
            pretty_da(
              mnth[i], thresh[i], nA[i], nB[i], nC[i], nD[i],
              alpha = spend_obj$spend[i]
            )
          ) %>%
            mutate(alpha = spend_obj$spend[i])
        }
    })
    
    dat1 <-
      dat1 %>%
      inner_join(
        .,
        da_stats,
        c("mnth", "thresh")
      )
    
    n_post_join <- nrow(dat1)
    if (n_post_join != n_reports) {
      stop("join not 1-1")
    }
    
    dat1
  }




# first signif
bcpnn_mult_comp_signif <-
  dat11 %>%
  group_by(thresh) %>%
  dplyr::filter(ci_lo > 0) %>%
  arrange(mnth) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(
    dte = as_date(paste0(mnth, "-01")),
    est_name = paste0(est_name, "(MCadj)")
  )




# thresh_scale <- rev(hcl.colors(length(thresholds) + 1, "Inferno"))[-1]
thresh_scale <- rev(hcl.colors(length(thresholds), "SunsetDark"))

bind_rows(
  bcpnn_signif,
  bcpnn_mult_comp_signif,
  tibble(
    dte = as_date("2014-12-01"),
    est_name = "MAXSPRT",
    thresh = 0.05
  )
) %>%
  ggplot(., aes(x = dte, y = est_name, col = factor(thresh))) +
  geom_point() +
  # geom_label(aes(label = factor(thresh))) +
  geom_label_repel(aes(label = factor(thresh))) +
  scale_colour_manual(values = thresh_scale) +
  theme_bw()



thresh_scale <- rev(hcl.colors(4, "Inferno"))[-1]


bind_rows(
  bcpnn_signif,
  bcpnn_mult_comp_signif,
  tribble(
    ~thresh, ~dte,      
    0.1,   "2017-03-01",
    0.075, "2014-12-01",
    0.05,  "2014-12-01",
    0.025, "2016-08-01",
    0.01,  "2017-05-01"
  ) %>%
    mutate(dte = as_date(dte), est_name = "MAXSPRT")
) %>%
  arrange(thresh) %>%
  ggplot(., aes(x = dte, y = thresh, col = est_name)) +
  geom_point() +
  geom_path(aes(group = est_name)) +
  # geom_label(aes(label = factor(thresh))) +
  # geom_label_repel(aes(label = est_name)) +
  scale_colour_manual(values = thresh_scale) +
  theme_bw()



bind_rows(
  bcpnn_signif,
  bcpnn_mult_comp_signif
)








plt <-
  dat1 %>%
  mutate(dte = as_date(paste0(mnth, "-01"))) %>% # make mnth date first of month
  ggplot(
    ., 
    aes(
      dte, 
      est, 
      ymax = ci_hi, 
      ymin = ci_lo, 
      # col = factor(thresh), 
      # fill = factor(thresh), 
      group = factor(thresh))
  ) %+%
  geom_hline(aes(yintercept = 0), col = "grey50") %+% # null value
  geom_line() %+%
  geom_point() %+%
  geom_ribbon(alpha = 0.05, lty = 2)  %+%
  labs(
    title  = "4. Pelvic mesh compared to hernia mesh (multiple comparison adjusted",
    subtitle = "(Calculations made on cumulative monthly data)",
    y = "IC statistic using BCPNN MCMC method (null value = 0)",
    x = "Date",
    col = "Pain score\nthreshold",
    fill = "Pain score\nthreshold"
  ) %+%
  # scale_colour_manual(values = thresh_scale) %+%
  # scale_fill_manual(values = thresh_scale) %+%
  theme_classic() 

print(plt)


cat("\n\n**Data for pain score >= ", 0.05, "**\n\n", sep = "")
dat1 %>%
  dplyr::filter(thresh == 0.05) %>%
  kable(., digits = 3) %>%
  print(.)




# ---- end ----

stopCluster(cl) 

sessionInfo()


