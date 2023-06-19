
# ---- libs ----


suppressPackageStartupMessages({
  library("readr")
  library("dplyr")
  library("tidyr")
  library("lubridate") # way to handle dates better than default R way
  library("knitr")
  library("foreach")
  library("arrow") # read/write parquet files
})



# here are the functions written for this project
source("r/_funcs.R")


# ---- consts ----

# arbitrarily, let's go with minimum cell count of 1 
arbitrary_cell_min <- 1

# these are the thresholds for pain_topic to be pain == TRUE
# thresholds <- c(0.010, 0.025, 0.05, 0.075, 0.100, 0.150) 
thresholds <- seq(0.010, 0.100, by = 0.005) 


col_pal <- c("cyan4", "darkorange", "purple", "dodgerblue")


target_lst <-
  list(
    "pelvic_mesh",
    "pelvic_mesh",
    "pelvic_mesh",
    "hernia_mesh",
    c("hernia_mesh", "other_mesh")
  )

compar_lst <-
  list(
    "hernia_mesh",
    c("hernia_mesh", "other_mesh"),
    c("hernia_mesh", "other_mesh", "other_device"),
    "other_mesh",
    "other_device"
  )




# ---- load_dat ----


clean_data_cols <-
  cols(
    Report_ID = col_double(),
    Date = col_date(format = ""),
    pain_word = col_logical(),
    pain_topic = col_double(),
    type = col_character()
  )

clean_data <- read_csv("dat/clean_data.csv", col_types = clean_data_cols)


# ---- check_dups ----


### all look like duplicates
inner_join(
  clean_data,
  clean_data %>% 
    group_by(Report_ID) %>% 
    summarise(n = n(), .groups = "drop") %>% 
    dplyr::filter(n > 1),
  "Report_ID"
) %>%
  arrange(Report_ID) %>%
  print(., n = nrow(.))


# make dup free
clean_data <-
  clean_data %>%  
  arrange(Report_ID, Date, desc(pain_word)) %>% # pain first in dups
  group_by(Report_ID) %>% 
  dplyr::filter(row_number() == 1) %>%
  ungroup(.) %>%  
  arrange(Date, Report_ID, desc(pain_word), desc(pain_topic))


clean_data %>%
  dplyr::filter(type == "other_mesh") %>%
  # select(Report_ID) %>%
  write_csv("out/other_mesh_ids.csv")



# ---- inspect ----


cat("First 10 rows of raw data\n")
clean_data %>%
  arrange(Date) %>%
  dplyr::filter(row_number() < 11) %>%
  kable(.)

# clean_data <-
#   clean_data %>%
#   dplyr::filter(
#     type %in% c("pelvic_mesh", "hernia_mesh")
#   )


clean_data %>%
  with(., table(type, pain_word)) %>%
  knitr::kable(.)

clean_data %>%
  with(., table(type, pain_topic >= 0.05)) %>%
  knitr::kable(.)



# These are the device groups and subgroups.
clean_data %>% 
  group_by(type) %>% 
  summarise(count = n()) %>%
  kable(.)

cat("\n\n## Histogram of `pain_word` (boolean) v `pain_topic` (score)")

clean_data %>%
  ggplot(., aes(pain_topic, fill = pain_word)) +
  geom_histogram(bins = 30) +
  scale_fill_manual(values = col_pal[1:2]) +
  theme_bw()

clean_data %>% 
  dplyr::filter(
    type %in% c("pelvic_mesh", "hernia_mesh", "other_mesh", "other_device")
  ) %>%
  ggplot(., aes(pain_topic, fill = pain_word)) +
  geom_histogram(bins = 30) +
  scale_fill_manual(values = col_pal[1:2]) +
  facet_wrap(~ type, scales = "free_y") +
  theme_bw()


# ---- create_analysis_data ----


### testing: example 1
# Use pelvic mesh as group 1 and hernia_mesh mesh devices as group 2. 
# The value of interest is the pain topic, being above the threshold of 0.05. 
# (i.e. 5% of the document contains words from the pain topic)
# You can adjust the topic threshold if you want to balance the groups more. 
# A higher topic_threshold will look for documents that discuss "pain" more, and 
# hence find less pain documents.

# get_signal_dat(
#   g1 = "pelvic_mesh",
#   g2 = "hernia_mesh",
#   pain_type = "pain_topic", 
#   thresh = 0.05,
#   cell_min = 1,
#   cumul = TRUE,
#   verbose = FALSE
# ) %>%
#   bind_cols(., thresh = 0.05)



# takes ~ 20 sec
cumul_dat <-
  foreach(i = 1:length(target_lst), .combine = bind_rows, .packages = "dplyr") %do% {
    foreach(th_j = thresholds, .combine = bind_rows, .packages = "dplyr") %do% {
      
      get_signal_dat(
        g1 = target_lst[[i]],
        g2 = compar_lst[[i]],
        pain_type = "pain_topic", 
        thresh = th_j,
        cell_min = 1,
        cumul = TRUE,
        verbose = FALSE
      ) %>%
        mutate(
          grps = 
            paste0(
              paste(target_lst[[i]], collapse = "/"), 
              " v ",
              paste(compar_lst[[i]], collapse = "/")
            ),
          dat_type = "cumulative", 
          thresh = th_j
        ) %>%
          select(grps, dat_type, thresh, everything())
      
    }
  }
cumul_dat

# takes ~ 20 sec
snpsh_dat <-
  foreach(i = 1:length(target_lst), .combine = bind_rows, .packages = "dplyr") %do% {
    foreach(th_j = thresholds, .combine = bind_rows, .packages = "dplyr") %do% {
      
      get_signal_dat(
        g1 = target_lst[[i]],
        g2 = compar_lst[[i]],
        pain_type = "pain_topic", 
        thresh = th_j,
        cell_min = 1,
        cumul = FALSE,
        verbose = FALSE
      ) %>%
        mutate(
          grps = 
            paste0(
              paste(target_lst[[i]], collapse = "/"), 
              " v ",
              paste(compar_lst[[i]], collapse = "/")
            ),
          dat_type = "snapshot", 
          thresh = th_j
        ) %>%
        select(grps, dat_type, thresh, everything())
      
    }
  }
snpsh_dat



# ---- check_analysis_data ----

nrow(cumul_dat)
if (nrow(cumul_dat) != nrow(snpsh_dat)) {
  stop("logic of creating analysis data producing different # rows in data")
}


chk_start_vals <- 
  inner_join(
    cumul_dat %>% 
      group_by(grps, dat_type, thresh) %>% 
      dplyr::filter(row_number() == 1) %>%
      ungroup(.),
    snpsh_dat %>% 
      group_by(grps, dat_type, thresh) %>% 
      dplyr::filter(row_number() == 1) %>%
      ungroup(.),
    c("grps", "thresh")
  ) %>%
    mutate(
      mnth_same = (mnth.x == mnth.y),
      counts_same = (nA.x = nA.y) & (nB.x = nB.y) & (nC.x = nC.y) & (nD.x = nD.y) 
    )

chk_start_vals %>%
  select(grps, thresh, dat_type.x, dat_type.y, mnth_same, counts_same)

with(chk_start_vals, table(mnth_same, counts_same, useNA = "ifany"))


# check the first + second row in snapshot == second row in cumulative data 
inner_join(
  cumul_dat %>% 
    group_by(grps, thresh) %>% 
    dplyr::filter(row_number() %in% 1:2) %>%
    ungroup(.),
  snpsh_dat %>% 
    group_by(grps, thresh) %>% 
    dplyr::filter(row_number() %in% 1:2) %>%
    ungroup(.),
  c("grps", "thresh", "mnth")
) 




# ---- export ----


# all spontaneous report analysis data
sra_dat <- 
  bind_rows(
    cumul_dat,
    snpsh_dat
  )

sra_dat %>%
  write_parquet(., sink = "dat/sra_dat.parquet")



