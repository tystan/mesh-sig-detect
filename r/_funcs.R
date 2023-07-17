
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



# ---- maxsprt_calcs ----


# null hypoth RR: RR0 = 1 
max_sprt_stat <- function(c_n, n, z, RR0 = 1) {
  
  # the simple version of this equation is on page 70/71 of
  # Kulldorf et al. (2011) A Maximized Sequential Probability Ratio Test for
  # Drug and Vaccine Safety Surveillance. Sequential Analysis, 30(1): 58-78.
  
  z_rr <- z / RR0
  if ((z_rr) * c_n / (n - c_n) <= 1) {
    max_llr <- 0
  } else if (c_n == n) {
    max_llr <- n * log(1 + z_rr)
  } else {
    max_llr <- 
      c_n * log(c_n / n) + 
      (n - c_n) * log((n - c_n) / n) - 
      c_n * log(1 / (z_rr + 1)) - 
      (n - c_n) * log(z_rr / (z_rr + 1))
  }
  return(max_llr)
  
}


max_sprt_stat_ <- Vectorize(max_sprt_stat)

rr_est <- function(c_n, n, z) {
  return(z * c_n / (n - c_n))
}

rr_est_ <- Vectorize(rr_est)


E_case <- function(c_n, n, z) {
  # return(z * c_n / (n - c_n))
  return(NULL)
}

max_sprt_stat(4, 5, (1 + 10) / (4 + 12))
rr_est(4, 5, (1 + 10) / (4 + 12))
max_sprt_stat(7, 9, (11 + 6) / (16 + 6))
rr_est(7, 9, (11 + 6) / (16 + 6))





get_maxsprt_cv <- function(max_n, per_look, z, alpha = 0.05, min_events = 1) {
  
  # has to be per period (not cumulative) for Sequential::CV.Binomial()
  # i.e. test performed at 3 events then when 3 more events come in requires
  # GroupSizes = c(3, 3)
  gs_seq <- rep(per_look, floor(max_n / per_look))
  if (sum(gs_seq) != max_n) { # if doesn't go to max_n, add at end for last look
    gs_seq <- c(gs_seq, max_n - sum(gs_seq)) 
  }
  
  if (max_n <= 500) {
    return(
      Sequential::CV.Binomial(
        N = max_n,
        alpha = alpha,
        M = min_events,
        z = z, 
        GroupSizes = gs_seq
      )$cv
    )
  } else { # for larger max n, use the sampling approach to get CV
    return(
      EmpiricalCalibration::computeCvBinomial(
        groupSizes = gs_seq,
        z = z,
        minimumEvents = min_events,
        alpha = alpha,
        sampleSize = 1e+06
      )
    )
  }
  
}




