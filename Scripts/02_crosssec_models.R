# ==============================================================================
# 02_crossec_models.R
# do vaccines (RZV vs ZVL) differ at a given cross-sectional time T?
# ==============================================================================



# preliminaries ================================================================

# do baseline fluorospot values differ based on prior zoster vacc?
# (yes, but weak association)
joineddat_wide %>%
  names %>%
  .[grepl("_0_", .)] %>%
  paste0("`", ., "` ~ PriorZosterVac") %>%
  map(function(x) {
    lm(x, joineddat_wide) %>%
      summary()
    } ) %>%
  set_names(
    joineddat_wide %>%
      names %>%
      .[grepl("_0_", .)]
  )

# are fluorospot values at 1yr associated with baseline values?
# (yes, strong association)
joineddat_wide %>%
  names %>%
  .[grepl("_12_", .)] %>%
  paste0("`", ., "` ~ `", gsub("_12_", "_0_", .), "`") %>%
  map(function(x) {
    lm(x, joineddat_wide) %>%
      summary()
  } ) %>%
  set_names(
    joineddat_wide %>%
      names %>%
      .[grepl("_0_", .)]
  )




# grid to setup all cross-sectional comparisons ================================

# grid of all outcomes, adjusting for baseline
modeling_matrix <-
  expand_grid(Time =
              c(1, 3, 12, 24, 48, 60, 999),
              Stim = c("VZV", "gE"),
              Marker = c("IFN-γ", "IL-2", "DP")) %>%
  rowwise() %>%
  mutate(Outcome = paste0(c(Marker, Time, Stim), collapse = "_"),
         Baseline =  paste0(c(Marker, 0, Stim), collapse = "_"))

# above but with 0
modeling_matrix_add0 <-
  expand_grid(Time =
                c(0, 1, 3, 12, 24, 48, 60, 999),
              Stim = c("VZV", "gE"),
              Marker = c("IFN-γ", "IL-2", "DP")) %>%
  rowwise() %>%
  mutate(Outcome = paste0(c(Marker, Time, Stim), collapse = "_"),
         Baseline =  paste0(c(Marker, 0, Stim), collapse = "_"))




# key functions ================================================================
# functions to model outcome at post-vaccination T, adjusting for baseline  

run_models <-
  function(Outcome, Baseline, covariates) {
    
    base_formula <-
      paste0("log10(`", Outcome,
             "` + 1) ~ Vaccine + log10(`", Baseline, "` + 1)",
             covariates) %>%
      as.formula
    
    tryCatch(
      lm(base_formula,
         joineddat_wide) %>%
        summary %>%
        coef %>%
        .[2, ],
      error = function(e) { rep(NA_real_, 4) %>%
          set_names(c("Estimate", "Std. Error", "t value", "Pr(>|t|)")) })
  }

# function to clean-up results for publication
clean_results_df <-
  function(results) {
    results
  results$T_FDR <- p.adjust(results$`Pr(>|t|)`, "fdr")
  results$SignifStudyVac1 <- if_else(results$T_FDR < 0.05, "*", "")
  
  results <-
    results %>%
    rowwise() %>%
    transmute(Time, Stim, Marker,
              Estimate = 10^(-Estimate),
              P = `Pr(>|t|)`,
              FDR = T_FDR,
              Signif = SignifStudyVac1) %>%
    mutate_if(is.numeric, formatC, digits = 3, format = "g", flag = "#")
  
  results
}




# core set of results, no covariate adjustment =================================

# run models comparing vaccine (RZV vs ZVL)
results_core <-
  pmap(modeling_matrix %>% 
         select(c("Outcome", "Baseline")) %>%
         mutate(covariates = ""),
       run_models) %>%
  bind_rows %>%
  bind_cols(modeling_matrix, .) %>%
  clean_results_df

# descriptive stats to join to results table
mean_values_crosssect <-
  joineddat_wide %>%
  ungroup() %>%
  select(Vaccine, contains("_VZV"), contains("_gE")) %>%
  group_by(Vaccine) %>%
  mutate_at(2:ncol(.), .funs = function(x) { log10(x + 1) } ) %>%
  summarize_all(summaryMeanSD, digits = 2, na.rm = T) %>%
  pivot_longer(cols = 2:ncol(.), values_to = "Mean (SD)") %>%
  separate(col = name, into = c("Marker", "Time", "Stim"), sep = "_") %>%
  pivot_wider(names_from = Vaccine, values_from = `Mean (SD)`,
              names_prefix = "Mean (SD)")

# final results
results_core %>%
  mutate(Time = as.numeric(Time) %>% as.character()) %>%
  left_join(mean_values_crosssect, .,
            by = c("Stim", "Marker", "Time")) %>%
  mutate(Time = gsub("999", "Peak", Time)) %>%
  arrange(Stim, Marker != "IL-2", Marker != "IFN-γ", Time) %>%
  mutate_all(., function(x) { x[is.na(x)] <- ""; x }) %>%
  View(., title = "Cross-Sectional Results (Suppl. Table 2)")




# run above, but covariate-adjusted ============================================

# e.g., age-adjusted
results_ageadj <-
  pmap(modeling_matrix %>% 
         select(c("Outcome", "Baseline")) %>%
         mutate(covariates = " + AgeGrp"),
       run_models) %>%
  bind_rows %>%
  bind_cols(modeling_matrix, .) %>%
  clean_results_df

# sex-adj 
results_sexadj <-
  pmap(modeling_matrix %>% 
         select(c("Outcome", "Baseline")) %>%
         mutate(covariates = " + Sex"),
       run_models) %>%
  bind_rows %>%
  bind_cols(modeling_matrix, .) %>%
  clean_results_df

# example of comparing two sets of results:
# highly correlated effect estimates and p-values
# so robust analysis to age effect
plot_grid(
  ggplot(data = NULL,
         aes(x = results_core$Estimate %>% as.numeric,
             y = results_sexadj$Estimate %>% as.numeric)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, lty = 3) +
    xlab("Core Results") + ylab("Age-Adjusted") +
    theme_classic(),
  ggplot(data = NULL,
       aes(x = results_core$P %>% as.numeric %>% log10,
           y = results_sexadj$P %>% as.numeric  %>% log10)) +
  geom_point() +
    geom_abline(intercept = 0, slope = 1, lty = 3) +
  xlab("Core Results") + ylab("Age-Adjusted") +
    theme_classic()
)




# look at *group* ==============================================================
# group (A, B, C, D); prior zoster vaccine info vs. binary vaccine ID (ZVL, RZV)

# run group models
results_grpinteract <-
  pmap(modeling_matrix%>% 
         select(c("Outcome", "Baseline")),
       
       function(Outcome, Baseline) {
         base_formula <-
           paste0("log10(`", Outcome, "` + 1) ~ Group + log10(`",
                  Baseline, "` + 1)") %>% as.formula
         
         tryCatch(
           anova(
             lm(base_formula, joineddat_wide),
             lm(update(base_formula, . ~ . - Group), joineddat_wide)) %>%
             .$`Pr(>F)` %>% .[2]
           ,
           error = function(e) { rep(NA_real_, 1) })
       }) %>%
  unlist %>%
  bind_cols(modeling_matrix, ANOVA_F = .)
results_grpinteract$ANOVA_FDR <- p.adjust(results_grpinteract$ANOVA_F, "fdr")
results_grpinteract$GroupEffect <- if_else(results_grpinteract$ANOVA_F < 0.05, "*", "")

# not presented in final manuscript
results_grpinteract %>%
  rowwise() %>%
  transmute(Time, Stim, Marker,
            P = ANOVA_F,
            FDR = ANOVA_FDR,
            Signif = GroupEffect) %>%
  mutate_if(is.numeric, formatC, digits = 2, format = "g", flag = "#") %>%
  View(title = "Differences by Group (ABCD)")




# look at vac-by-zoster interaction ============================================
# significant interaction would imply prior zoster vaccine modifies associations
# labeled "alt" because similar idea as looking at group (A, B, C, D)

# run models for prior vac x vaccine ID interaction
results_interact_ALT <-
  pmap(modeling_matrix %>% 
         select(c("Outcome", "Baseline")),
       
       function(Outcome, Baseline) {
         
         base_formula <-
           paste0("log10(`", Outcome, "` + 1) ~ Vaccine*PriorZosterVac + log10(`",
                  Baseline, "` + 1)") %>% as.formula
         
         tryCatch(
           anova(
             lm(base_formula, joineddat_wide),
             lm(update(base_formula, . ~ . - PriorZosterVac - Vaccine:PriorZosterVac),
                joineddat_wide)) %>%
             .$`Pr(>F)` %>% .[2]
           ,
           error = function(e) { rep(NA_real_, 1) })
       }
       
       ) %>%
  unlist %>%
  bind_cols(modeling_matrix, ANOVA_F = .)
results_interact_ALT$ANOVA_FDR <- p.adjust(results_interact_ALT$ANOVA_F, "fdr")
results_interact_ALT$GroupEffect <- if_else(results_interact_ALT$ANOVA_FDR < 0.05, "*", "")

# results for age-by-vaccine interaction
results_interact_ALT %>%
  rowwise() %>%
  transmute(Time, Stim, Marker,
            P = ANOVA_F,
            FDR = ANOVA_FDR,
            Signif = GroupEffect) %>%
  mutate_if(is.numeric, formatC, digits = 2, format = "g", flag = "#") %>%
  View(title = "Prior Vac-by-Vaccine Interaction")



# agegrp-by-zoster interaction =================================================
# significant interaction would imply age group (A or B) modifies associations
# checking since age a priori known impact on vaccine outcomes

# run models for agegroup interaction
results_interact_ALTAGE <-
  pmap(modeling_matrix %>% 
         select(c("Outcome", "Baseline")),
       
       function(Outcome, Baseline) {
         
         base_formula <-
           paste0("log10(`", Outcome, "` + 1) ~ Vaccine*AgeGrp + log10(`",
                  Baseline, "` + 1)") %>% as.formula
         
         tryCatch(
           anova(
             lm(base_formula, joineddat_wide),
             lm(update(base_formula, . ~ . - AgeGrp - Vaccine:AgeGrp),
                joineddat_wide)) %>%
             .$`Pr(>F)` %>% .[2]
           ,
           error = function(e) { rep(NA_real_, 1) })
       }
       
  ) %>%
  unlist %>%
  bind_cols(modeling_matrix, ANOVA_F = .)
results_interact_ALTAGE$ANOVA_FDR <- p.adjust(results_interact_ALT$ANOVA_F, "fdr")
results_interact_ALTAGE$GroupEffect <- if_else(results_interact_ALT$ANOVA_FDR < 0.05, "*", "")

# results for agegroup-by-vaccine interaction
results_interact_ALTAGE %>%
  rowwise() %>%
  transmute(Time, Stim, Marker,
            P = ANOVA_F,
            FDR = ANOVA_FDR,
            Signif = GroupEffect) %>%
  mutate_if(is.numeric, formatC, digits = 2, format = "g", flag = "#") %>%
  View(title = "Age-by-Vaccine Interaction")



