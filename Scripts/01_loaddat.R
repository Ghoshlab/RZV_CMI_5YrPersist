# ==============================================================================
# 01_loaddat.R
# import fluorospot & metadata from excel files
# ==============================================================================




# fluorospot data (outcome) ====================================================

# count of cells positive out of mil PBMC
# 2,216 x 8
fluoro <-
  read_excel("./Data/2020.11.05 GSK Fluorospot Data (submitted to Choo 11.9.20).xlsx") %>% 
  transmute(PID, Group, Month, Stim, # nicer names
            AgeGroup = Age,
            `IFN-γ` = `IFNg (Spots/1 million PBMC)` %>% as.numeric,
            `IL-2` = `IL2 (Spots/1 million PBMC)` %>% as.numeric,
            DP = `Double Positive (Spots/1 million PBMC)` %>% as.numeric)

# interpreting groups
# design matrix from "GSK protocol consensus.5.20.13"
design <-
  tribble(~Arm, ~PriorZosterVac, ~StudyVac1, ~StudyVac2,
          "A", 0, "Live", "none",
          "B", 0, "gE", "gE",
          "C", 1, "Live", "none",
          "D", 1, "gE", "gE") %>%
  
  # alternative names for vaccines:
  # Hu/su or gE --> RZV   &&   live or VZ --> ZVL
  mutate_all(.funs = function(x) {
    x %>% gsub("Live", "ZVL", .) %>% gsub("gE", "RZV", .)
  })




# participant metadata (age, ancestry, sex) ====================================

metadat <-
  map(1:3, function(sheetnum) {
    read_excel("./Data/Randomization Enrollment Log w Demog.xlsx",
               sheet = sheetnum) %>%
      select(1:8) %>%
      filter(Group %in% c("I", "II", "III"))
  }) %>%
  bind_rows() %>%
  mutate(PID =
            gsub("III", "3", `Subj ID`) %>%
            gsub("II", "2", .) %>%
            gsub("I", "1", .)) %>%
  rename(Cohort = Group) %>%
  set_names(., names(.) %>% gsub(" ", "", .))




# join fluoro with metadata ====================================================

# one subject (3-62D) is in file, but doesn't have any non-missing fluoro vals
fluoro %>%
  .[apply(., 1,
          function(x) { is.na(x) %>% any() }), ]

setdiff(metadat$PID, fluoro$PID) # one partic (1-05B) in metadat but not flouro 
setdiff(fluoro$PID,metadat$PID) # no partic w/ fluoro dat but no metadata

# "joineddat" tibble with both
joineddat <-
  fluoro %>% 
  filter(PID != "3-62D") %>%
  left_join(design, by = c("Group" = "Arm")) %>%
  left_join(metadat, by = "PID") %>%
  group_by(PID) 

# since peak of RZV and ZVL are at different times,
# also define Th1 responses at a "peak" time (dummy value Month = 999)
peak_responses <-
  joineddat %>%
  filter( (StudyVac1 == "RZV" & Stim == "gE" & Month == 3) |
            (StudyVac1 == "ZVL" & Stim == "gE" & Month == 1) |
            (StudyVac1 == "RZV" & Stim == "VZV" & Month == 3) |
            (StudyVac1 == "ZVL" & Stim == "VZV" & Month == 1) ) %>%
  mutate(Month = 999)

# add "Time" variable for interpretable month labels
joineddat <-
  joineddat %>%
  bind_rows(peak_responses) %>%
  mutate(Vaccine = StudyVac1) %>%
  mutate(Time = Month %>% as.factor %>%
           fct_recode(`0` = "0", `1m` = "1", `3m` = "3",
                      `1yr` = "12", `2yr` = "24", `4yr` = "48", `5yr` = "60",
                      Peak = "999") )

# swapping between month and descriptive times for plots
# time points available: 0, 1, 3, 12, 24, 48, 60
labels_month2time <-
  joineddat$Month %>% unique %>%
  set_names(joineddat$Time %>% unique)
labels_time2month <-
  names(labels_month2time) %>% set_names(labels_month2time)




# check sample sizes / missingness rates =======================================

# sample sizes by group x time
joineddat %>%
  filter(!is.na(`IFN-γ`) & !is.na(`IL-2`) & !is.na(DP)) %>%
  group_by(Vaccine, Time) %>%
  filter(!duplicated(PID)) %>%
  tally() %>%
  pivot_wider(Vaccine,
              names_from = c(Time),
              values_from = n)

# checking percent missing rates, by arm (A, B, C, D)
# will show NA in a column if column has nothing missing
joineddat %>%
  mutate(Available = !is.na(`IFN-γ`) & !is.na(`IL-2`) & !is.na(DP)) %>%
  group_by(Month) %>%
  group_split() %>%
  map( ~ .x %>%
         group_by(Available, Group) %>%
         tally() %>%
         pivot_wider(id_cols = c("Available"),
                     names_from = "Group", values_from = n) %>%
         apply(., 2, function(x) { x/sum(x)} )
  ) %>%
  set_names(labels_month2time)

# checking percent missing rates, by age group
joineddat %>%
  mutate(Available = !is.na(`IFN-γ`) & !is.na(`IL-2`) & !is.na(DP)) %>%
  bind_rows(., mutate(., AgeGroup = "All")) %>%
  group_by(Month) %>%
  group_split() %>%
  map( ~ .x %>%
         group_by(Available, AgeGroup) %>%
         tally() %>%
         pivot_wider(id_cols = c("Available"),
                     names_from = "AgeGroup", values_from = n) %>%
         apply(., 2, function(x) { x/sum(x)} )
  ) %>%
  set_names(labels_month2time)




# final dataset ================================================================

joineddat_wide <-
  joineddat %>%
  pivot_wider(id_cols = c("PID", "Group", "Cohort", "Age", "AgeGrp",
                          "Race", "Ethnicity", "Sex",
                          "Vaccine", "PriorZosterVac"),
              names_from = c("Month", "Stim"),
              values_from = c("IFN-γ", "IL-2", "DP"))
  


