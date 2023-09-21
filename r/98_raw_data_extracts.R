library("readr")
library("stringr")


# raw data
rawd <-
  read_csv(
    "https://raw.githubusercontent.com/curtis-murray/MedicalDevicesNLP/master/data/all_reports/all_reports_df.csv"
  )


# cleaned data
clean_data_cols <-
  cols(
    Report_ID = col_double(),
    Date = col_date(format = ""),
    pain_word = col_logical(),
    pain_topic = col_double(),
    type = col_character()
  )

clean_data <- read_csv("dat/clean_data.csv", col_types = clean_data_cols)

# join sources
nrow(rawd)
nrow(clean_data)
rawd_j <-
  inner_join(
    clean_data,
    rawd,
    c("Report_ID" = "Report number")
  )
nrow(rawd_j)

rawd_j


# print options

with(rawd_j, table(`Reported event outcome`, useNA = "ifany")) %>% kable(.)
with(rawd_j, table(`Device classification`, useNA = "ifany")) %>% kable(.)
with(rawd_j, table(`Report source category`, useNA = "ifany")) %>% kable(.)


rawd_j %>%
  dplyr::filter(
    !pain_word, 
    pain_topic >= 0.01, pain_topic <= 0.03, 
    # `Device classification` == "Class IIa", 
    `Reported event outcome` == "Injury",
    `Report source category` == "Industry"
  )  %>%
  arrange(desc(`Report date`)) 


rawd_j %>%
  dplyr::filter(
    pain_topic >= 0.05, 
    type == "other_device",
    # `Device classification` == "Class IIa", 
    `Reported event outcome` == "Injury",
    `Report source category` == "Industry"
  ) %>%
  arrange(desc(`Report date`)) 


rawd_j %>%
  dplyr::filter(
    pain_topic >= 0.05, 
    !(type %in% c("pelvic_mesh", "other_device")),
    # `Device classification` == "Class IIa", 
    `Reported event outcome` == "Injury",
    `Report source category` == "Consumer"
  )%>%
  arrange(`Report date`) 


rawd_j %>%
  dplyr::filter(
    pain_topic >= 0.01, pain_topic <= 0.02 ,
    type %in% c("other_mesh"),
    # `Device classification` == "Class IIa", 
    # `Reported event outcome` == "Injury",
    `Report source category` == "Industry"
  )%>%
  arrange(desc(`Report date`)) 

rawd_j %>%
  dplyr::filter(
    pain_topic <= 0.01,
    type %in% c("other_device"),
    # `Device classification` == "Class IIa", 
    `Reported event outcome` != "Injury",
    `Report source category` != "Industry"
  )%>%
  arrange(desc(`Report date`)) 


# cat(paste(colnames(rawd_j), collapse = '",\n"'))

cprint <- c(
  "Report_ID",
  "Report date",
  "Classification" = "type",
  "Device" = "Device classification",
  "P('pain'|doc)" = "pain_topic",
  # "Trade name",
  # "Sponsor",
  # "Manufacturer",
  "ARTG no." = "ARTG number",
  # "GMDN term",
  # "Sterile",
  # "Single use",
  # "Model number",
  # "Software version",
  "Event" = "Reported event outcome",
  "Source" = "Report source category",
  "Event type",
  "Description" = "Event description"
)

ex_ids <- c(
  45265 ,  # industry, other device, pain > 0.05
  40917,  # consumer, pelvic mesh, pain > 0.05
  37537 , # industry, other device, 0.01 < pain < 0.03
  45432,  # consumer, hernia_mesh, pain > 0.05
  36797,  # industry, other_mesh,  0.01 < pain < 0.03
  45624,  # Health Professional, other device,  0.01 > pain 
  44402   # Other, other device,  0.01 > pain 
)

rawd_j %>%
  dplyr::filter(`Report_ID` %in% ex_ids) %>%
  select(all_of(cprint)) %>%
  mutate(Description = str_trunc(Description, 150)) %>%
  arrange(`Report date`) %>%
  kable(., digits = 3)
