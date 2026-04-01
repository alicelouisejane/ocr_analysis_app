library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(janitor)
# normalise and create variables OCR analysis

process_ocr_data <- function(df_input,dna_input,num_islets_per_well=70) {

  #reshape the dataframe from the traditional that was given and uplaoded to RedCap- assume this remains consistant
  df<- df_input %>%
    mutate(across(where(is.character), ~na_if(.x, ""))) %>%
    remove_empty(c("rows", "cols")) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "donor_id")

  dna<-dna_input %>%
    mutate(across(where(is.character), ~na_if(.x, ""))) %>%
    remove_empty(c("rows", "cols")) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    rownames_to_column(var = "donor_id")

  # Convert transposed data back into a dataframe and combined with dna
  ocr_clean <- df %>%
    mutate(meta_islets_well = num_islets_per_well, .after = donor_id)%>%
    merge(dna,all.x=T,by="donor_id") %>% # ================================Merge with DNA file
    select(donor_id,meta_islets_well,meta_dna_cont,everything()) %>%
    mutate(across(-donor_id, as.numeric))
  # Pivot data to long format
  ocr_clean_long <- pivot_longer(ocr_clean, cols = -c(donor_id,meta_islets_well,meta_dna_cont), names_to = "variable", values_to = "value")



  #Normalise to DNA and Normalise to DNA and baseline
  # normalize to DNA, where missing use : A typical human cell contains about 6 picograms (pg) of DNA,
  # IEQs contain approx 1500 cells meaning that 1 IEQ would contain approximately 1500*6= 9000pg of DNA
  # ~70 IEQs = 9000x70 = 0.63 µg of DNA
  #put OCR data into phases in order to calculate AUC per phase
  ocr_clean_long_norm<-ocr_clean_long %>%
    group_by(donor_id) %>%
    separate(variable,c("variable","time"),sep = "_") %>%
    mutate(across(time,as.numeric)) %>%
    mutate(value_norm_dna=value/ifelse(is.na(meta_dna_cont),(meta_islets_well*9000)/1e6,meta_dna_cont)) %>% #
    mutate(value_norm_dna_baseline=value_norm_dna/value_norm_dna[time==26.89]) %>%
    mutate(phase=ifelse(time<35.5,"2.8 mM glucose",
                        ifelse(time>=35.5 & time<86.62,"16.7 mM glucose",
                               ifelse(time>=86.62 & time<154.75,"Oligomycin",
                                      ifelse(time>=154.75 & time<205.88,"FCCP",
                                             ifelse(time>=205.88 & time<=265.42,"Rotenone/Antimycin",NA))))))

  ### data QC - identify outliters on the triplicate and flag on raw value only
  outliers_time_withintrip <- ocr_clean_long %>%
    separate(donor_id, c("donor_id", "replicate"), sep = "_", remove = FALSE) %>%
    separate(variable, c("measurement", "time"), sep = "_") %>%
    mutate(time = as.numeric(time)) %>%
    select(donor_id, replicate, time, meta_dna_cont, measurement, value) %>%
    arrange(donor_id, time, replicate) %>%
    group_by(donor_id, time, measurement) %>%
    mutate(
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      Lower_Bound = Q1 - 1.5 * IQR,
      Upper_Bound = Q3 + 1.5 * IQR,
      Is_Outlier = if_else(value < Lower_Bound | value > Upper_Bound, TRUE, FALSE)
    ) %>%
    ungroup()

  flag_donor_id_outliers_time_withintrip <- outliers_time_withintrip %>%
    filter(Is_Outlier) %>%
    select(donor_id, replicate, time, measurement, value, meta_dna_cont, Is_Outlier)

  # removal of within triplicate outliers ?- need to have a better method here see ppt
  if(nrow(flag_donor_id_outliers_time_withintrip)>0){
    print(paste(
      "Donor IDs with Tukey-identified outliers within triplicate:",
      paste(unique(flag_donor_id_outliers_time_withintrip$donor_id), collapse = ", ")
    ))
    }else{
    print("No donors ids have a Tukey identified outlier within the triplicate.")
  }

  #Merge flag of any outliers identifed in triplicate -
  ocr_clean_long_norm_qc<-ocr_clean_long_norm %>%
    pivot_longer(
      cols = starts_with("value"),
      names_to = "measurement",
      values_to = "value"
    ) %>%
    left_join(flag_donor_id_outliers_time_withintrip)

  # maybe at this point show a map of any potential error replicates?
  # need a way to do this error replicate properly


  ############# CALCULATIONS
  #performed on the long dataframe
  # calculations should only be performed on the raw or normalised to DNA values:
  # Normalised to DNA and baseline data is essentially a fold change so it would be inappropriate to put it into one of the equations to calculate the OCR bioenergetic metrics as those equations were derived on the actual data not fold changes in order to give the output.. if that makes sense
  #
  # links below supporting this:
  # “Calculation of bioenergetic parameters requires normalization to correct for differences in numbers of cells, organoids, or amount of biological material.”  https://www.sciencedirect.com/science/article/pii/S2666166721000939
  # “The Wave software used to view XF data has a built in “Baseline” feature that transforms absolute XF rate data to a relative (%) scale. Most often, the baseline is set to the rate just prior to the first injection. Baselining data is most appropriate when attempting to minimize slight well to well differences in rate due to variations in cell seeding or proliferation, and is helpful to visualize changes in rates from acute treatments/injections.” https://www.agilent.com/cs/library/technicaloverviews/public/Methods_and_Strategies_for_Normalizing_Tech_Overview_022118.pdf

  calculations<- ocr_clean_long_norm_qc %>%
    mutate(measurement=gsub("value_","",measurement)) %>%
    pivot_wider(names_from = measurement, values_from = c(value,Is_Outlier)) %>%
    rename(value=value_value) %>%
    group_by(donor_id) %>%

    # Non-Mitochondrial Oxygen Consumption

    # minimum of rotinone pase
    mutate(calc_nonmito_oc_raw=min(value[phase=="Rotenone/Antimycin"],na.rm = T)) %>%
    mutate(calc_nonmito_oc_norm_dna=min(value_norm_dna[phase=="Rotenone/Antimycin"],na.rm = T)) %>%

    # Basal Respiration
    mutate(calc_basal_resp_raw=value[time==26.89]-calc_nonmito_oc_raw) %>%
    mutate(calc_basal_resp_norm_dna=value_norm_dna[time==26.89]-calc_nonmito_oc_norm_dna) %>%


    # ATP-Linked Respiration
    mutate(calc_atp_resp_raw=value[time==26.89]-min(value[phase=="Oligomycin"],na.rm = T)) %>%
    mutate(calc_atp_resp_norm_dna=value_norm_dna[time==26.89]-min(value_norm_dna[phase=="Oligomycin"],na.rm = T)) %>%

    # Proton Leak
    mutate(calc_proton_leak_raw=min(value[phase=="Oligomycin"],na.rm = T)-calc_nonmito_oc_raw) %>%
    mutate(calc_proton_leak_norm_dna=min(value_norm_dna[phase=="Oligomycin"],na.rm = T)-calc_nonmito_oc_norm_dna) %>%

    # Maximal Glucose Respiration
    mutate(calc_max_gluc_resp_raw=max(value[phase=="16.7 mM glucose"],na.rm = T)-calc_nonmito_oc_raw) %>%
    mutate(calc_max_gluc_resp_norm_dna=max(value_norm_dna[phase=="16.7 mM glucose"],na.rm = T)-calc_nonmito_oc_norm_dna) %>%

    # Stimulated glucose respiration
    mutate(calc_stim_gluc_resp_raw=max(value[phase=="16.7 mM glucose"],na.rm = T)-value[time==26.89]) %>%
    mutate(calc_stim_gluc_resp_norm_dna=max(value_norm_dna[phase=="16.7 mM glucose"],na.rm = T)-value_norm_dna[time==26.89]) %>%

    # Maximal Respiration
    mutate(calc_max_resp_raw=max(value[phase=="FCCP"],na.rm = T)-calc_nonmito_oc_raw) %>%
    mutate(calc_max_resp_norm_dna=max(value_norm_dna[phase=="FCCP"],na.rm = T)-calc_nonmito_oc_norm_dna) %>%

    # Spare Respiratory Capacity
    mutate(calc_spare_cap_raw=max(value[phase=="FCCP"],na.rm = T)-value[time==26.89]) %>%
    mutate(calc_spare_cap_norm_dna=max(value_norm_dna[phase=="FCCP"],na.rm = T)-value_norm_dna[time==26.89]) %>%

    # Glucose-Stimulated Respiration (if glucose is added at a specific time and measured subsequently)
    mutate(calc_gluc_stim_oci_raw=max(value[phase=="16.7 mM glucose"],na.rm = T)/value[time==26.89]) %>%
    mutate(calc_gluc_stim_oci_norm_dna=max(value_norm_dna[phase=="16.7 mM glucose"],na.rm = T)/value_norm_dna[time==26.89]) %>%

    # Glucose stimulated spared capacity
    mutate(calc_gluc_stim_sparecap_raw=max(value[phase=="FCCP"],na.rm = T)-max(value[phase=="16.7 mM glucose"],na.rm = T)) %>%
    mutate(calc_gluc_stim_sparecap_norm_dna=max(value_norm_dna[phase=="16.7 mM glucose"],na.rm = T)-max(value_norm_dna[phase=="FCCP"],na.rm = T)) %>%

    # ATP-linked in the glucose-stimulated state (last rate before Oligomycin) # ie max time is 78.01
    mutate(calc_atp_resp_gluc_raw =value[time == max(time[phase == "16.7 mM glucose"], na.rm = TRUE)] - min(value[phase == "Oligomycin"], na.rm = TRUE)) %>%
    mutate(calc_atp_resp_gluc_norm_dna = value_norm_dna[time == max(time[phase == "16.7 mM glucose"], na.rm = TRUE)] - min(value_norm_dna[phase == "Oligomycin"], na.rm = TRUE)) %>%

    # bioenergetic health index
    # https://www.agilent.com/Library/usermanuals/Public/BHI_Report_Generator_User_Guide_RevA.pdf
    mutate(calc_bhi_raw=(calc_spare_cap_raw*calc_atp_resp_raw)/(calc_nonmito_oc_raw*calc_proton_leak_raw)) %>%
    mutate(calc_bhi_norm_dna=(calc_spare_cap_norm_dna*calc_atp_resp_norm_dna)/(calc_nonmito_oc_norm_dna*calc_proton_leak_norm_dna)) %>%

    # glucose stimulated bioenergetic health index
    mutate(calc_bhi_gluc_raw=(calc_gluc_stim_sparecap_raw*calc_atp_resp_gluc_raw)/(calc_nonmito_oc_raw*calc_proton_leak_raw)) %>%
    mutate(calc_bhi_gluc_norm_dna=(calc_spare_cap_norm_dna*calc_atp_resp_gluc_norm_dna)/(calc_nonmito_oc_norm_dna*calc_proton_leak_norm_dna)) %>%

    select(donor_id,starts_with("calc_")) %>%
    unique()

  nice_names <- c(
    calc_atp_resp          = "ATP-linked respiration",
    calc_atp_resp_gluc     = "ATP-linked (glucose-stim.)",
    calc_basal_resp        = "Basal respiration",
    calc_bhi               = "BHI",
    calc_bhi_gluc          = "BHI (glucose stim)",
    calc_gluc_stim_oci     = "Glucose-stimulated OCI (max/basal)",
    calc_gluc_stim_sparecap = "Spare capacity (glucose-stim.)",
    calc_stim_gluc_resp    = "Stimulated glucose respiration",
    calc_max_resp          = "Max respiration (FCCP)",
    calc_max_gluc_resp    = "Max glucose respiration (16.7Mm)",
    calc_nonmito_oc        = "Non-mitochondrial OCR",
    calc_proton_leak       = "Proton leak",
    calc_spare_cap         = "Spare capacity"
  )

  ################## AVERAGING TRIPLICATES to give final
  # Calculated parameters
  ocr_final_calculations<- calculations %>%
    separate(donor_id,c("donor_id","repeat"),"_") %>%
    pivot_longer(
      cols = starts_with("calc"),
      names_to = "measurement",
      values_to = "value"
    ) %>%
    mutate(normalisation_type = case_when(
      grepl("value$|raw$", measurement) ~ "Raw",
      grepl("value_norm_dna$|_norm_dna$", measurement) ~ "Norm DNA",
      grepl("value_norm_dna_baseline|_norm_dna_baseline", measurement) ~ "Norm DNA and Baseline",
      TRUE ~ "Unknown"  # Catch-all for any unexpected cases
    )) %>%
    mutate(measurement=gsub("_value.*|_raw.*|_norm.*","",measurement)) %>%
    group_by(donor_id,measurement,normalisation_type) %>%
    summarise(value=mean(value)) %>%
    ungroup()

  # phase on raw ocr data in order to calculate AUC
  ocr_final<- ocr_clean_long_norm_qc %>%
    separate(donor_id,c("donor_id","repeat"),"_") %>%
    group_by(donor_id,measurement,time) %>%
    summarise(mean=mean(value)) %>%
    ungroup() %>%
    mutate(normalisation_type = case_when(
      grepl("value$|raw$", measurement) ~ "Raw",
      grepl("value_norm_dna$|_norm_dna$", measurement) ~ "Norm DNA",
      grepl("value_norm_dna_baseline|_norm_dna_baseline", measurement) ~ "Norm DNA and Baseline",
      TRUE ~ "Unknown"  # Catch-all for any unexpected cases
    )) %>%
    mutate(measurement=gsub("_value.*|_raw.*|_norm.*","",measurement)) %>%
    ungroup() %>%
    mutate(phase=ifelse(time<=35.5,"2.8 mM glucose",
                        ifelse(time>35.5 & time<=86.62,"16.7 mM glucose",
                               ifelse(time>86.62 & time<=154.75,"Oligomycin",
                                      ifelse(time>154.75 & time<=205.88,"FCCP",
                                             ifelse(time>205.88 & time<=265.42,"Rotenone/Antimycin",NA)))))) %>%

    group_by(donor_id,phase,normalisation_type) %>%
    mutate(auc=pracma::trapz(time,mean)) %>%
    ungroup() %>%
    mutate(measurement=time) %>%
    select(-time)

  ################## Outliers after averaging
  #Flag if baseline time is <50 or >450- raw trace no normalisation
  outliers_baseline_ocr_final<-ocr_final %>%
    filter(normalisation_type=="Raw") %>%
    filter(measurement==1.38) %>%
    unique() %>%
    group_by(donor_id) %>%
    mutate(FLAG_firstvallessthat50_raw=ifelse(mean<50,TRUE,FALSE)) %>%
    mutate(FLAG_firstvalmorethan450_raw=ifelse(mean>450,TRUE,FALSE)) %>%
    select(donor_id,FLAG_firstvallessthat50_raw,FLAG_firstvalmorethan450_raw) %>%
    unique()

  # Tukey outlier on phase AUC
  outliers_phase_ocr_final<-ocr_final %>%
    filter(normalisation_type=="Norm DNA and Baseline") %>%
    select(donor_id,phase,auc) %>%
    unique() %>%
    group_by(phase) %>%
    mutate(Q1 = quantile(auc, 0.25, na.rm = TRUE),
           Q3 = quantile(auc, 0.75, na.rm = TRUE),
           IQR = Q3 - Q1,
           lower_bound = Q1 - 1.5 * IQR, # use IQR methods for small within sample testing as moer robust
           upper_bound = Q3 + 1.5 * IQR,
           FLAG_outlierphase_whennormdnabaseline = ifelse(auc < lower_bound | auc > upper_bound,TRUE,FALSE))%>%
    ungroup() %>%
    select(donor_id,phase,auc,FLAG_outlierphase_whennormdnabaseline)

  # Tukey outlier on calculated parameters
    outliers_ocr_final_calculations<-ocr_final_calculations %>%
      filter(normalisation_type=="Norm DNA") %>%
      group_by(measurement) %>%
      mutate(Q1 = quantile(value, 0.25, na.rm = TRUE),
             Q3 = quantile(value, 0.75, na.rm = TRUE),
             IQR = Q3 - Q1,
             lower_bound = Q1 - 1.5 * IQR, # use IQR methods for small within sample testing as moer robust
             upper_bound = Q3 + 1.5 * IQR,
             FLAG_outliercalc_whennormdna = ifelse(value < lower_bound | value > upper_bound,TRUE,FALSE)) %>%
      ungroup() %>%
      select(donor_id,measurement,value,FLAG_outliercalc_whennormdna)

# plot distribution of calculations to enable visualisation of outliers
set.seed(1)

    outliers_ocr_final_calculations2 <- outliers_ocr_final_calculations %>%
      dplyr::mutate(x_jit = 1 + runif(dplyr::n(), -0.12, 0.12))

    set.seed(1)

    outliers_ocr_final_calculations2 <- outliers_ocr_final_calculations %>%
      dplyr::mutate(
        x_jit = 1 + runif(dplyr::n(), -0.10, 0.10),
        donor_label = ifelse(FLAG_outliercalc_whennormdna, donor_id, NA_character_)
      )

    p1 <- ggplot(outliers_ocr_final_calculations2, aes(x = "", y = value)) +
      geom_boxplot(
        width = 0.25,
        outlier.shape = NA,
        fill = "grey92",
        color = "grey40"
      ) +
      geom_jitter(
        aes(color = FLAG_outliercalc_whennormdna),
        width = 0.12,
        height = 0,
        size = 2,
        alpha = 0.8
      ) +
      facet_wrap(
        ~measurement,
        scales = "free_y",
        ncol = 3,
        labeller = labeller(measurement = nice_names)
      ) +
      scale_color_manual(
        values = c(`TRUE` = "red3", `FALSE` = "forestgreen"),
        name = "Outlier?"
      ) +
      labs(
        x = NULL,
        y = "Value",
        subtitle = "Calculated metrics normalised to DNA"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing = unit(1, "lines")
      )

# plot distribution of phases
    outliers_phase_ocr_final2 <- outliers_phase_ocr_final %>%
      dplyr::mutate(
        donor_label = ifelse(FLAG_outlierphase_whennormdnabaseline, donor_id, NA_character_)
      )

    p2 <- ggplot(outliers_phase_ocr_final2, aes(x = "", y = auc)) +
      geom_boxplot(
        width = 0.25,
        outlier.shape = NA,
        fill = "grey92",
        color = "grey40"
      ) +
      geom_jitter(
        aes(color = FLAG_outlierphase_whennormdnabaseline),
        width = 0.12,
        height = 0,
        size = 2,
        alpha = 0.8
      ) +
      facet_wrap(~phase, scales = "free_y", ncol = 3) +
      scale_color_manual(
        values = c(`TRUE` = "red3", `FALSE` = "forestgreen"),
        name = "Outlier?"
      ) +
      labs(
        x = NULL,
        y = "AUC",
        subtitle = "Phase AUC using data normalised to DNA and baseline"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing = unit(1, "lines")
      )
########################## FINAL DATA SET OUTPUTS
  ocr_final_traces<-ocr_final %>%
    mutate(across(measurement,as.numeric)) %>%
    dplyr::group_by(donor_id) %>%
    rename(time=measurement) %>%
    arrange(time) %>%
    arrange(donor_id) %>%
    ungroup() %>%
    full_join(outliers_baseline_ocr_final)

  ocr_final_calculations<- ocr_final_calculations %>%
    full_join(outliers_baseline_ocr_final)

# PLOT OF TRACES
  # Segment data
  segments_data <- data.frame(
    start_time = c(0, 35.5, 86.62, 154.75, 205.88),
    end_time = c(35.5, 86.62, 154.75, 205.88, 265.42),
    condition = c("2.8mM Glucose", "16.7mM Glucose", "Oligomycin (5uM)",
                  "FCCP (3uM)", "Rotenone/Antimycin A (5uM)"),
    color = c("#7774FF",
              "#92C797",
              "#C274FF",
              "#85A3BE",
              "#F974FF")
    ,
    text_position = c(18, 60, 120, 180, 235)
  )
  plot_data <- ocr_final_traces %>%
    filter(normalisation_type == "Norm DNA and Baseline") %>%   # change if needed
    mutate(donor_id = as.factor(donor_id)) %>%
    mutate(flag=ifelse(FLAG_firstvallessthat50_raw==T | FLAG_firstvalmorethan450_raw==T,T,F))
  flag_labels <- plot_data %>%
    filter(flag == TRUE) %>%
    group_by(donor_id) %>%
    slice_max(time, n = 1) %>%
    ungroup()

trace_plot<- ggplot(plot_data, aes(x = time, y = mean, group = donor_id)) +

    # background segments
    geom_rect(
      data = segments_data,
      aes(xmin = start_time, xmax = end_time,
          ymin = -Inf, ymax = Inf),
      fill = segments_data$color,
      alpha = 0.25,
      inherit.aes = FALSE
    ) +

    # traces
    geom_line(aes(color = flag), linewidth = 0.7, alpha = 0.8) +
  geom_line(data=filter(plot_data,flag==T),color="magenta", linewidth = 0.7, alpha = 0.8) +

    scale_color_manual(
      values = c(`FALSE` = "black", `TRUE` = "magenta"),
      name = "Flag"
    ) +

    labs(
      title = "Norm DNA and baseline",
      x = "Time (minutes)",
      y = "OCR (pmol/min)"
    ) +

    coord_cartesian(clip = "off") +

    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(10, 40, 10, 10),
      legend.position = "bottom"
     ) +
    # geom_text_repel(
    #   data = flag_labels,
    #   aes(label = donor_id, color = flag),
    #   direction = "y",
    #   nudge_x = 10,
    #   hjust = 0,
    #   segment.color = "grey40",
    #   box.padding = 0.3,
    #   point.padding = 0.2,
    #   size = 3,
    #   show.legend = FALSE
    # ) +
  geom_text(
    data = segments_data,
    aes(x = text_position, y = Inf, label = condition),
    inherit.aes = FALSE,
    vjust = 1.5,
    size = 4
  )


# retun list of outputs
return(list(
  ocr_traces = ocr_final_traces,
  ocr_calculations = ocr_final_calculations,
  outlier_calcs = outliers_ocr_final_calculations,
  outlier_phase = outliers_phase_ocr_final,
  outlier_baseline = outliers_baseline_ocr_final,
  outlier_triplicate = flag_donor_id_outliers_time_withintrip,
  outlier_calcs.plot = p1,
  outlier_phase.plot = p2,
  ocr_final = ocr_final
))

}
