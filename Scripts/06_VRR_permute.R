# ==============================================================================
# 06_VRR_permute.R
# calculate vaccine response rates (VRR) and test differences by vaccine
# ==============================================================================




# calculate VRR ================================================================

# calculate outcome at time T / outcome at baseline
# store with variable name "RATIO."
VR_dat <- joineddat_wide

generate_ratios <-
  function(Outcome, Baseline) {
    VR_dat <<- 
      VR_dat %>%
      rowwise() %>%
      mutate(!!sym(paste0("RATIO.", Outcome, "")) :=  !!sym(Outcome) /
               (!!sym(Baseline) + 0.1) )
  }

pmap(modeling_matrix_add0 %>%
       select(Outcome, Baseline),
     generate_ratios)

VR_dat <-
  VR_dat %>%
  select(1:10, contains("RATIO")) %>%
  pivot_longer(cols = names(VR_dat) %>% .[grepl("RATIO", .)],
               names_to = c("Measurement", "Month", "Stim"), names_sep = "_",
               names_transform = list(Month = as.double)) %>%
  mutate(Measurement = gsub("RATIO.", "", Measurement) %>%
           factor %>% fct_rev)




# calculate VRRs ===============================================================

# define VRR at time T as the proportion of participants with ratio >= 2
VR_dat_summary <-
  VR_dat %>%
  group_by(Vaccine, Measurement, Month, Stim) %>%
  filter(!is.na(value)) %>%
  summarize(N = n(),
            Ratio = median(value),
            VRR = sum(value >= 2) / N * 100)

# consider Delta := VRR (RZV) - VRR (ZVL) at time T
delta_vrr <- 
  VR_dat_summary %>%
  pivot_wider(id_cols = c("Measurement", "Month", "Stim"),
              names_from = Vaccine, values_from = VRR) %>%
  mutate(Delta = RZV - ZVL)

# vector of observed deltas
delta_vrr_obs <- 
  delta_vrr %>%
  .$Delta




# compare VRRs at time T by vaccine ============================================

# compare observed delta to 10,000 permuted values of delta
# in which the vaccine labels are shuffled
VR_dat_permutation <-
  map_dfc(1:10000,
          function(seed) {
            delta_vrr_perm <-
              VR_dat %>%
              group_by(Measurement, Month, Stim) %>%
              filter(!is.na(value)) %>%
              group_split() %>%
              map(.f = function(df) {
                set.seed(seed)
                df$Vaccine <- sample(df$Vaccine)
                df }) %>%
              bind_rows() %>%
              group_by(Vaccine, Measurement, Month, Stim) %>%
              summarize(N = n(),
                        VRR = sum(value >= 2) / N * 100) %>%
              pivot_wider(id_cols = c("Measurement", "Month", "Stim"),
                          names_from = Vaccine, values_from = VRR) %>%
              mutate(Delta = RZV - ZVL) %>%
              .$Delta

  (abs(delta_vrr_perm) >= abs(delta_vrr_obs))
    })

# convert the TRUE/FALSE values to p-values
delta_vrr$Pval <-
  VR_dat_permutation %>% rowSums() %>% `+`(., 1) %>% `/`(., 10001)

# omit baseline, then FDR correct
delta_vrr_results <-
  delta_vrr %>%
  filter(Month != 0)

delta_vrr_results$FDR <- p.adjust(delta_vrr_results$Pval, 'fdr')
delta_vrr_results$Signif <- if_else(delta_vrr_results$FDR < 0.05, "*", "")




# final table of VRR results ===================================================

delta_vrr_results %>%
  mutate(Month = gsub("999", "Peak", Month)) %>%
  mutate_if(is.numeric,
            function(val) {
              if_else(val == 1/10001,
                      "<0.0001",
                      formatC(val, digits = 3, format = "g", flag = "#") ) }) %>%
  arrange(Stim, Measurement, Month) %>%
  select(3, 2, 1, 4:9) %>%
  View("VRR Values and Permutation Results, Suppl. Table 5")



