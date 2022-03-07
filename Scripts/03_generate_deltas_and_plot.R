# ==============================================================================
# 03_generate_deltas_and_plot.R
# calculate each participant's fluoro value at time T, minus baseline
# & plot deltas over time
# ==============================================================================




# generate each set of deltas ==================================================
# based on model matrix (script 02)
# store in "joineddat_wide_BLadj" (wide fmt) and "joineddat_BLadj" (long fmt)

# function to loop through each outcome
generate_deltas <-
  function(Outcome, Baseline) {
    joineddat_wide_BLadj <<-
      joineddat_wide_BLadj %>%
      rowwise() %>%
      mutate(!!sym(paste0("Δ", Outcome)) := !!sym(Outcome) - !!sym(Baseline))
    ""
  }

# run above function
joineddat_wide_BLadj <- joineddat_wide
pmap(
  modeling_matrix_add0 %>%
    select(Outcome, Baseline),
  generate_deltas
)

# covert wide output from above into long format for plotting
joineddat_BLadj <-
  joineddat_wide_BLadj %>%
  select(1:10, contains("Δ")) %>%
  pivot_longer(
    cols = names(joineddat_wide_BLadj) %>% .[grepl("Δ", .)],
    names_to = c("Measurement", "Month", "Stim"), names_sep = "_",
    names_transform = list(Month = as.double)
  )

# also for plotting, calculate signed-log10 value
joineddat_BLadj <-
  joineddat_BLadj %>%
  mutate(log10val = sign(value) * log10(abs(value + 1))) %>%
  mutate(Measurement = as.factor(Measurement) %>% fct_rev())




# plot mean DELTA (log-scale) over time ========================================
# stimname = gE or VZV

plot_delta_by_time_log10 <-
  function(stimname, asterisklocation = 2.75) {
    filtered_dat <-
      joineddat_BLadj %>%
      filter(Stim == stimname & !(Month %in% c("999", 999))) %>%
      mutate(Month = as.numeric(Month)) %>%
      arrange(Month) %>%
      mutate(Month = as.factor(Month) %>% fct_inorder())
    ggplot(
      filtered_dat,
      aes(Month, y = log10val)
    ) +
      geom_hline(aes(yintercept = 0)) +
      stat_summary(
        geom = "line", color = "#888888", show.legend = F,
        fun = median, aes(group = Vaccine, lty = Vaccine),
      ) +
      geom_boxplot2(aes(fill = Vaccine),
        notch = F, outlier.size = 0.1, outlier.shape = NA,
        position = position_dodge(preserve = "single"),
        width.errorbar = 0, alpha = 1, width = 0.4
      ) +
      geom_text(
        data = results_core %>%
          mutate(Measurement = paste0("Δ", Marker) %>% as.factor() %>% fct_rev()) %>%
          mutate(Month = as.numeric(Time)) %>%
          filter(!(Month %in% c(999)) & Stim == stimname) %>%
          arrange(Month) %>%
          mutate(Month = as.factor(Month) %>% fct_inorder()),
        aes(Month, label = Signif), y = asterisklocation, size = 5
      ) +
      facet_grid(Measurement ~ ., switch = "y", scales = "free") +
      theme_linedraw() +
      theme(
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      scale_linetype_manual(values = c(5, 5)) +
      scale_fill_manual(name = "Vaccine", values = palette_vaccine) +
      scale_x_discrete("Time",
                       labels = labels_time2month[levels(filtered_dat$Month)]) +
      scale_y_continuous(
        bquote(atop("Change in SFC/10"^6 ~
                      " from Baseline", "(log"[10] ~ "scales)")),
        breaks = seq(-2, 3), labels = c(-100, -10, 0, 10, 100, 1000)
      ) +
      ggtitle(paste0(stimname, "-Stimulation"))
  }




# run for each stimname ========================================================
# 800 x 600

plot_grid(
  plot_delta_by_time_log10("gE"),
  plot_delta_by_time_log10("VZV"),
  labels = c("A", "B")
)

ggsave("Fig1.png", width = 8, height = 6)

