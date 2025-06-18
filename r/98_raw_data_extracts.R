library("readr")
library("stringr")


# raw data
rawd <-
  read_csv(
    ### previously
    # "https://raw.githubusercontent.com/curtis-murray/MedicalDevicesNLP/master/data/all_reports/all_reports_df.csv"
    ### now saved locally for snapshot
    "dat/all_reports_df.csv"
  )

# rawd <- vroom::vroom("dat/all_reports_df.csv")
# vroom::problems(rawd)



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
    c("Report_ID" = "Report number"),
    relationship = "many-to-many"
  )
nrow(rawd_j)

rawd_j


# print options

with(rawd_j, table(`Reported event outcome`, useNA = "ifany")) %>% kable(.)
with(rawd_j, table(`Device classification`, useNA = "ifany")) %>% kable(.)
with(rawd_j, table(`Report source category`, useNA = "ifany")) %>% kable(.)

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
  "Event type",
  "Source" = "Report source category",
  "Description" = "Event description"
)




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
    # pain_topic >= 0.05, 
    type %in% c("pelvic_mesh"),
    # `Device classification` == "Class IIa", 
    # `Reported event outcome` == "Injury",
    `Report source category` == "Industry"
    # `pain_topic` > 0
  )%>%
  arrange(pain_topic, `Report date`) %>%
  select(all_of(cprint)) %>%
  print(., n = Inf)

rawd_j %>%
  dplyr::filter(`Report_ID` %in% c(30531, 26969, 32535)) %>%
  pull(`Report_ID`)
rawd_j %>%
  dplyr::filter(`Report_ID` %in% c(30531, 26969, 32535)) %>%
  pull(`Event description`)


rawd_j %>%
  dplyr::filter(
    # pain_topic >= 0.05, 
    type %in% c("pelvic_mesh"),
    # `Device classification` == "Class IIa", 
    # `Reported event outcome` == "Injury",
    `Report source category` == "Health Professional"
    # `pain_topic` > 0
  )%>%
  arrange(pain_topic, `Report date`) %>%
  select(all_of(cprint)) %>%
  print(., n = Inf)

rawd_j %>%
  dplyr::filter(`Report_ID` %in% c(43792, 36034 , 35033 )) %>%
  pull(`Report_ID`)
rawd_j %>%
  dplyr::filter(`Report_ID` %in% c(43792, 36034 , 35033)) %>%
  pull(`Event description`)
 

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



ex_ids <- c(
  45265,  # industry, other device, pain > 0.05
  40917,  # consumer, pelvic mesh, pain > 0.05
  30531,  # industry, pelvic mesh, pain < 0.05
  43792,  # Health Professional, pelvic mesh, pain == 0
  # 30315,  # industry, pelvic mesh, pain > 0
  37537,  # industry, other device, 0.01 < pain < 0.03
  45432,  # consumer, hernia_mesh, pain > 0.05
  36797,  # industry, other_mesh,  0.01 < pain < 0.03
  45624,  # Health Professional, other device,  0.01 > pain 
  44402   # Other, other device,  0.01 > pain 
)

# truncated
rawd_j %>%
  dplyr::filter(`Report_ID` %in% ex_ids) %>%
  select(all_of(cprint)) %>%
  mutate(Description = str_trunc(Description, 150)) %>%
  arrange(`Report date`) %>%
  kable(., digits = 3)


# summary of combos
rawd_j %>%
  dplyr::filter(`Report_ID` %in% ex_ids) %>%
  select(all_of(cprint)) %>%
  select(-Description, -`ARTG no.`) %>%
  arrange(`Report date`) %>%
  kable(., digits = 3)


rawd_j %>%
  dplyr::filter(`Report_ID` %in% ex_ids) %>%
  pull(`Event description`)


# keep plain spaces obviously!
unicode_ws_chrs <-
  c(
    "\u0009", "\u000A", "\u000B", "\u000C", "\u000D", "\u0085", "\u00A0",
    # "\u1361",
    "\u1680", "\u2000", "\u2001", "\u2002", "\u2003", "\u2004", "\u2005", "\u2006",
    "\u2007", "\u2008", "\u2009", "\u200A", "\u2028", "\u2029", "\u202F", "\u205F",
    "\u3000"
  )

unicode_ws_chrs <- paste0("[", paste(unicode_ws_chrs, collapse = ""), "]")


# description not truncated
rawd_j %>%
  dplyr::filter(`Report_ID` %in% ex_ids) %>%
  select(all_of(cprint)) %>%
  mutate(`Description` = str_replace_all(`Description`, unicode_ws_chrs, " ")) %>%
  mutate(`Description` = stringi::stri_enc_toascii(`Description`)) %>%
  mutate(Description = trimws(Description, "both", whitespace = unicode_ws_chrs)) %>%
  arrange(`Report date`) %>%
  kable(., digits = 3)
  
  



