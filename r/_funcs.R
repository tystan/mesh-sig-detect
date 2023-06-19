
# ---- funcs ----

# The below function is thanks to Curtis:
# https://github.com/curtis-murray/MedicalDevicesNLP

#' Function to get data for disproportionality analysis
#' 
#' @param group_1 any vector combination of "pelvic_mesh", "hernia_mesh",
#' "other_mesh", "other_device".
#' @param group_2 any vector combination of "pelvic_mesh", "hernia_mesh",
#' "other_mesh", "other_device".
#' @param pain_type either "pain_word" which detects if the word pain 
#' appears in a doucment, or "pain_topic", which will detect the pain 
#' topic in a document.
#' @param topic_threshold if pain_type == "pain_topic" then this is the 
#' value that the topic density must exceed for the topic to be registered.
#'  The default value of 0.05 will count documents with at least 5% of words 
#'  being from the pain topic.
#' @return Tibble of counts nA nB nC nD
#' @examples
#' get_data(
#'   group_1 = c("pelvic_mesh"), group_2 = c("hernia_mesh", "other_mesh"), 
#'   pain_type = "pain_topic", topic_threshold = 0.05
#' )
#' get_data(
#'   group_1 = c("pelvic_mesh"), group_2 = c("hernia_mesh"), 
#'   pain_type = "pain_word", topic_threshold = 0.05
#' )

get_data2 <- function(group_1, group_2, pain_type, topic_threshold = 0.05){
  
  intersection = intersect(group_1, group_2)
  if(length(intersection) > 0){
    warning(paste(intersection, " contained in both groups", sep = ""))
  }
  dat_out <-
    clean_data %>% 
    filter(type %in% c(group_1,group_2)) %>% 
    mutate(group = ifelse(type %in% group_1, "group_1", "group_2")) %>% 
    mutate(pain_topic = pain_topic >= topic_threshold) 
  
  dat_out <- dat_out[, c("Report_ID", "Date", "group", pain_type)]
  dat_out[["pain"]] <- dat_out[[pain_type]]
  dat_out[[pain_type]] <- NULL # rename and remocve
  
  # select(Report_ID, Date, group, pain = !!as.name(pain_type)) %>% 
  
  dat_out <- dat_out %>%
    mutate(group_tf = (group == "group_1")) %>% 
    mutate(nA = 1*(group_tf & pain),
           nB = 1*(group_tf & !pain),
           nC = 1*(!group_tf & pain),
           nD = 1*(!group_tf & !pain)) %>% 
    select(-group_tf)
  
}


# Ty Stanford (2022-05-24)

get_signal_dat <- 
  function(
    g1 = "pelvic_mesh", 
    g2 = c("hernia_mesh", "other_mesh"),
    pain_type = "pain_topic", 
    thresh = 0.05,
    cell_min = 1,
    cumul = TRUE, # want counts to that point in time?
    verbose = FALSE
  ) {
    
    
    d_g1_v_g2 <- 
      get_data2(
        group_1 = g1, 
        group_2 = g2, 
        pain_type = pain_type, 
        topic_threshold = thresh
      ) %>%
      mutate(Report_ID = as.character(Report_ID)) # make sure ID is discrete
    
    # check report_ids unique
    ids_df <- d_g1_v_g2 %>% select(Report_ID)
    ids <- pull(ids_df, Report_ID)
    if (length(ids) != length(unique(ids))) {
      message("The following IDs are repeated")
      ids_df %>% 
        group_by(Report_ID) %>% 
        summarise(n = n(), .groups = "drop") %>% 
        dplyr::filter(n > 1) %>%
        print(.)
      stop("duplicate IDs found")
    }
    
    # sort by date, create "YYYY-MM" variable for grouping
    d_g1_v_g2 <-
      d_g1_v_g2 %>%
      arrange(Date) %>%
      mutate(mnth = paste0(year(Date), "-", sprintf("%02.0f", month(Date))))
    
    if (verbose) {
      # peak
      d_g1_v_g2 %>% kable(.)%>% print(.)
    }
    
    # make each entry a monthly count
    d_g1_v_g2 <-
      d_g1_v_g2 %>%
      group_by(mnth) %>%
      summarise(across(matches("^n[A-D]"), sum)) %>%
      ungroup()
    
    # d_g1_v_g2 %>% print(. , n = nrow(.))
    
    # now make all counts cumulative to that point in time
    d_g1_v_g2_cumul <-
      d_g1_v_g2 %>%
      mutate(across(matches("^n[A-D]"), cumsum))
    
    if (verbose) {
      # go on, have a peak at what we're lookin at
      d_g1_v_g2 %>% kable(.) %>% print(.)
      # stop()
    }
    
    # now let's find a date with sufficient contingency counts as starting point
    date_min_dat <-
      d_g1_v_g2_cumul %>%
      pivot_longer(cols = matches("^n[A-D]")) # long format for life
    
    
    if (verbose) {
      # go on, have a peak at what we're lookin at
      date_min_dat %>% kable(.) %>% print(.)
    }
    
    date_min_dat <-
      date_min_dat %>%
      dplyr::filter(value >= cell_min) %>%
      # now need to find date when all nA to nD reach `arbitrary_cell_min`
      # so get max of mins
      group_by(name) %>%
      summarise(mnth = min(mnth)) %>%
      ungroup() 
    
    if (verbose) {
      # go on, have a peak at what we're lookin at
      date_min_dat %>% kable(.) %>% print(.)
    }
    
    date_min <-
      date_min_dat %>%
      pull(mnth) %>%
      max(.)
    
    if (verbose) {
      print(date_min)
      print(class(date_min))
    }
    
    if (cumul) { # want cumulative counts
      
      d_g1_v_g2 <- 
        d_g1_v_g2_cumul %>% 
        dplyr::filter(mnth >= date_min) # update data to date range wanted
      
    } else {  # want snapshot counts (but need cumulative until min count met)
      
      d_g1_v_g2 <-
        bind_rows(
          d_g1_v_g2_cumul %>% dplyr::filter(mnth <= date_min) %>%
            summarise(
              mnth = max(mnth),
              across(matches("^n[A-D]"), max)
            ),
          d_g1_v_g2 %>% dplyr::filter(mnth > date_min)
        )
      
    }
    
    # d_g1_v_g2 <-
    #   d_g1_v_g2 %>%
    #   arrange(mnth)
    
    # return
    return(d_g1_v_g2)
    
  }





