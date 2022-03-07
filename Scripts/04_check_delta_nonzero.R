# ==============================================================================
# 04_check_delta_nonzero.R
# among recipients of a given vaccine, do values at time T differ from baseline?
# ==============================================================================




# key functions ================================================================

# define delta = log10[time T] - log10[baseline val])
# null hypothesis: delta = 0
# vaccname = vaccine name in indicator variable (if RZV, will return ZVL)

run_models_delta0 <-
  function(Outcome, Baseline, vaccname) {
    
    joineddat_wide$Delta <-
      log10(joineddat_wide[ , Outcome, drop = T] + 1) -
      log10(joineddat_wide[ , Baseline, drop = T] + 1) 
    
    joineddat_wide <-
      joineddat_wide %>% filter(Vaccine == vaccname)
    
    base_formula <-
      paste0("Delta ~ 1") %>%
      as.formula
    
    tryCatch(
      lm(base_formula,
         joineddat_wide) %>%
        summary %>%
        coef %>%
        .[1, ] %>%
        formatC(digits = 3, format = "g", flag = "#") %>%
        c(., `Δ Mean (SD)` = summaryMeanSD(joineddat_wide$Delta, na.rm = T, digits = 3)),
      error = function(e) { rep(NA_real_, 4) %>%
          set_names(c("Estimate", "Std. Error", "t value", "Pr(>|t|)")) })
  }



# run models ===================================================================
# modeling matrix from script 02,
# but exclude peak and add vaccine

modmatrix_delta0 <-
  modeling_matrix %>% 
  filter(Time != 999) %>%
  mutate(vaccname = "RZV") %>%
  bind_rows(., mutate(., vaccname = "ZVL"))

# run models
results_returnto0 <-
  pmap(modmatrix_delta0 %>% 
         select(c("Outcome", "Baseline", "vaccname")),
       run_models_delta0) %>%
  bind_rows %>%
  bind_cols(modmatrix_delta0, .)

results_returnto0$FDR <- p.adjust(results_returnto0$`Pr(>|t|)`, "fdr")
results_returnto0$Signif <- if_else(results_returnto0$FDR < 0.05, "*", "")

results_returnto0 <-
  results_returnto0 %>%
  mutate_if(is.numeric, formatC, digits = 3, format = "g", flag = "#")



# final table ==================================================================

results_returnto0 %>%
  transmute(Stim, Marker, Time = as.numeric(Time) %>% as.character(),
            Vaccine = vaccname, `Δ Mean (SD)`, `p-value` = `Pr(>|t|)`, FDR, Signif) %>%
  pivot_wider(id_cols = 1:3, values_from = c(8, 5:7),
              names_from = Vaccine) %>%
  arrange(Stim, Marker != "IL-2", Marker != "IFN-γ", as.numeric(Time)) %>%
  View("Time T - Baseline Values, Supplementary Table 3")



