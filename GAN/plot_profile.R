library(ggplot2)

profile_raw <- read.csv(
  file = "/home/emilie/Documents/Fil_rouge/Code/Rproj/input_profiles/M1_pl1.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = ";"
)
profile_raw <- read_raw_profile(profile_raw)
number_of_dyes <- 4

plot_max_rfu <- max(as.numeric(unlist(profile_raw[, 2:(number_of_dyes + 1)])))
plot_min_rfu <- min(as.numeric(unlist(profile_raw[, 2:(number_of_dyes + 1)])))

plot_list <- vector("list", number_of_dyes)

for (dye in 1:number_of_dyes) {
  
  col_name <- paste0("channel_", dye)
  
  p <- ggplot(profile_raw, aes(x = index)) +
    geom_line(
      aes(y = .data[[col_name]]),
      color = "black",
      linewidth = 0.2,
      alpha = 0.7
    ) +
    coord_cartesian(ylim = c(plot_min_rfu, plot_max_rfu)) +
    ylab("RFU") +
    xlab("") +
    ggtitle(paste0("Channel ", dye)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 8),
      axis.text.x = element_text(size = 5)
    )
  
  plot_list[[dye]] <- p
}
plot_list[[1]]
plot_list[[2]]
plot_list[[3]]
plot_list[[4]]