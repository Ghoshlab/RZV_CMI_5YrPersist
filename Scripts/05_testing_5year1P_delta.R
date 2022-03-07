# ==============================================================================
# 05_testing_5year1P_delta.R
# calculate additional deltas, but focused on changes btwn peak/1yr --> 5yrs
# ==============================================================================




# calculate additional deltas ==================================================

fiveyrslope_dat <- joineddat_wide

# visualize changes over time in slope plot
ggplot(joineddat %>%
         filter(Month %in% c(12, 60)) %>%
         arrange(-Month) %>%
         mutate(Time = (Time)),
       aes(Time, `IL-2`)) +
  geom_point() +
  geom_line(aes(group = PID), alpha = 0.5) +
  facet_grid(Stim ~ Vaccine) +
  theme_cowplot()

# calculate "Delta5M1P": five years minus 'baseline' (1 yr or peak)
generate_5minus1yr_deltas <-
  function(Outcome, Baseline) {
    fiveyrslope_dat <<-
      fiveyrslope_dat %>%
      rowwise() %>%
      mutate(!!sym(paste0("Delta5M1P.", Outcome, "")) :=
               log10(!!sym(Outcome) + 1) -
               log10(!!sym(Baseline) + 1) )
  }




# calculate additional deltas ==================================================

# calculates five-year delta relative to argument 'time_point'
run_fiveyear_change_models <-
  function(time_point) {
  
  # modeling matrix, relative to 1 year
  modeling_matrix_5minus1yr <-
    expand_grid(Time = c(60),
                Stim = c("VZV", "gE"),
                Marker = c("IFN-γ", "IL-2", "DP")) %>%
    rowwise() %>%
    mutate(Outcome = paste0(c(Marker, Time, Stim), collapse = "_"),
           Baseline =  paste0(c(Marker, time_point, Stim), collapse = "_"))
  
  # calculate the deltas
  pmap(modeling_matrix_5minus1yr %>%
         select(Outcome, Baseline),
       generate_5minus1yr_deltas)
  
  # run models
  fiveyrslope_dat <-
    fiveyrslope_dat %>%
    select(1:10, contains("Delta5M1P")) %>%
    pivot_longer(cols = names(fiveyrslope_dat) %>% .[grepl("Delta5M1P", .)],
                 names_to = c("Measurement", "Month", "Stim"), names_sep = "_",
                 names_transform = list(Month = as.double)) %>%
    mutate(Measurement = gsub("Delta5M1P.", "", Measurement) %>%
             factor %>% fct_rev)
    
  # sanity check plots (prints, but not returned by function)
  delta_by_vac_fig <-
      fiveyrslope_dat %>% 
      group_by(Measurement, Month, Stim) %>%
      group_split() %>%
      map(~ ggplot(data = .x, aes(Vaccine, value, color = Vaccine)) +
            geom_quasirandom(width = 0.2) +
            ylab(paste0("Delta (Time: 5yr - ", time_point, ")\n",
                        first(.x$Measurement), ", ",
                        first(.x$Stim), "-Stim")) +
            theme_classic() +
            geom_hline(yintercept = 0, lty = 2) +
            scale_color_manual(values = palette_vaccine) +
            theme(legend.position = "none")
      )
    
    plot5y1p <- plot_grid(delta_by_vac_fig[[1]],
                      delta_by_vac_fig[[2]],
                      delta_by_vac_fig[[3]],
                      nrow = 1)
    print(plot5y1p)
    
    # run t-test on differences
    fiveyrslope_dat_summary <-
      fiveyrslope_dat %>%
      group_by(Measurement, Month, Stim) %>%
      filter(!is.na(value)) %>%
      group_split() %>%
      map_dbl(~ lm(value ~ Vaccine, data = .x) %>% summary %>% coef %>% .[2, 4])
  
  fiveyr_pvals <-
    fiveyrslope_dat %>%
    group_by(Measurement, Month, Stim) %>%
    group_keys() %>%
    bind_cols(Pvalue = fiveyrslope_dat_summary)
  
  fiveyr_FDR <-
    p.adjust(fiveyr_pvals$Pvalue, "fdr")
  
  # arrange results and summarize
  fiveyrslope_dat %>%
    group_by(Vaccine, Measurement, Month, Stim) %>%
    summarise(MeanSD = summaryMeanSD(value, na.rm = T))  %>% 
    ungroup() %>%
    pivot_wider(id_cols = c("Measurement", "Month", "Stim"),
                names_from = Vaccine, values_from = MeanSD) %>%
    bind_cols(Pvalue = fiveyr_pvals$Pvalue,
              FDR = fiveyr_FDR)  %>%
    select(-Month)

}




# run tests ====================================================================

fiveyear_to_oneyear <- run_fiveyear_change_models(12) # 1 year
fiveyear_to_peak <- run_fiveyear_change_models(999) # peak

# final results table joining 1-yr and peak
left_join(fiveyear_to_peak,
          fiveyear_to_oneyear,
          by = c("Stim", "Measurement"),
          suffix = c(".Peak", ".1Year")) %>%
  mutate_if(.predicate = is.numeric,
            .funs = function(x) { formatC(x, digits = 2, format = "fg")}) %>%
  arrange(Stim, Measurement) %>%
  View("5-Year Decline, Suppl. Table 4")



