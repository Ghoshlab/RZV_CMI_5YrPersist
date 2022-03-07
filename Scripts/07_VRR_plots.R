# ==============================================================================
# 07_VRR_plots.R
# plot VRRs over time
# ==============================================================================




# plot VRR over time ===========================================================

# apply plotting function concept from script 03 to VRRs calculated in script 06
# stimname = gE or VZV
# asteriskposition = putting significance asterisk close to line on plot

plot_vrrs <-
  function(stimname, asteriskposition) {

  ggplot(VR_dat_summary %>%
           filter(Stim == stimname & Month != 999),
         aes(factor(Month), y = VRR)) +
    stat_summary(geom = "line", fun = median,
                 aes(group = Vaccine, lty = Vaccine, color = Vaccine)) +
    geom_point(size = 3, aes(shape = Vaccine)) +
    geom_text(data = delta_vrr_results %>%
                filter(Stim == stimname  & Month != 999),
              aes(as.factor(Month), label = Signif),
              y = asteriskposition, size = 5) +
    facet_grid(Measurement ~ ., switch = "y", scales = "fixed") +
    theme_linedraw() +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key.width = unit(0.3, "inches")) +
    scale_linetype_manual(values = c(RZV = 6, ZVL = 1)) +
    scale_shape_manual(values = c(ZVL = "Z", RZV = "R")) +
    scale_color_manual(values = palette_vaccine) +
    scale_x_discrete("Time", labels = labels_time2month) +
    scale_y_continuous("Vaccine Response Rate\n(% Recipients with Time t Response / Baseline ≥ 2)",
                       limits = c(0, 100),
                       breaks = c(0.0, 0.25, 0.50, 0.75, 1.0) %>% `*`(., 100),
                       labels = c(0.0, 0.25, 0.50, 0.75, 1.0) %>% `*`(., 100) %>% paste0(., "%"),
                       expand = c(0.05, 0.15)) +
    ggtitle(paste0("VRR, ", stimname, "-Stimulation"))
}




# run for each stimname ========================================================
# 800 x 600

plot_grid(
  plot_vrrs("gE", 1.01*100),
  plot_vrrs("VZV", 0.8*100),
  labels = c("A", "B"))
ggsave("Fig2.png", width = 8, height = 6)



