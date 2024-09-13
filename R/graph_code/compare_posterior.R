library(ggplot2)
library(showtext)
library(scales)

font_add_google(name = "Spectral", family = "Spectral")
showtext_auto()

# Define the color palette
color_palette <- c(
  jet = "#131516",
  accent_blue = "#0072B2",
  accent_green = "#009E73",
  coral_red = "#FF6347",
  sky_blue = "#87CEEB",
  graphite_gray = "#474A51"
)

theme_custom <- function(base_size = 12, base_family = "Spectral") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Text elements
      text = element_text(color = color_palette["jet"]),
      axis.title = element_text(size = rel(1.2)),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(color = color_palette["jet"]),
      
      # Plot title and subtitle
      plot.title = element_text(size = rel(1.5), hjust = 0.05),
      plot.subtitle = element_text(size = rel(1.2), hjust = 0.05),
      
      # Legend appearance
      legend.title = element_blank(),
      legend.text = element_text(size = rel(1)),
      legend.background = element_blank(),
      legend.key = element_rect(fill = "white", colour = color_palette["jet"]),
      legend.position = "bottom",
      
      # Panel background
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      
      plot.background = element_blank(),
      
      # Plot margins
      # plot.margin = margin(15, 15, 15, 15)
    )
}



compare_posteriors <- function(..., data_pts = NULL, vars = NULL, dodge_width = 0.5) {
  dots <- rlang::dots_list(..., .named = TRUE)

  draws <- lapply(dots, function(x) {
    if (class(x)[1] == "stanreg") {
      posterior::subset_draws(posterior::as_draws(x$stanfit),
                              variable = names(fixef(x))
      )
    } else if (class(x)[1] == "brmsfit") {
      brm_draws <- posterior::subset_draws(posterior::as_draws(x$fit),
                                           variable = paste0("b_", rownames(fixef(x)))
      )
      posterior::variables(brm_draws) <- stringr::str_split(posterior::variables(brm_draws), "_", simplify = T)[, 2]
      posterior::rename_variables(brm_draws, `(Intercept)` = Intercept)
    } else if (class(x)[1] == "CmdStanMCMC") {
       draws_out <- x$draws(variables = vars, format = "matrix") 
       x |>
         spread_draws(mu_study[i]) |>
         bind_rows(mu_draws) 
    } else {
      stop(paste0(class(x)[1], " objects not supported."))
    }
  })

  intervals <- lapply(draws, bayesplot::mcmc_intervals_data)

  if (!is.null(data_pts)) {
    data_names <- names(intervals[[1]])
    levels_ <- levels(unlist(data_pts[, 1]))
    
    data_pts <- tibble(
      model = "data",
      parameter = meta_dat[, "i"],
      outer_width = 0,
      inner_width = 0,
      point_est = "data",
      ll = 0,
      l = 0,
      m = meta_dat[, "yi"],
      h = 0,
      hh = 0,
      overall = "no"
    )
    intervals <- lapply(intervals, function(x) {
      x$parameter <- factor(x$parameter, levels = levels_)
      x
    })

  }


  combined <- dplyr::bind_rows(intervals, .id = "model") |>
    mutate(overall = if_else(parameter == "intercept", paste0(model, "_yes"), paste0(model, "_no")))

  combined |>
  ggplot(aes(x = m, y = parameter, group = model, alpha = model)) +
    geom_linerange(aes(xmin = l, xmax = h, color = overall), linewidth = 1.3, position = position_dodge(dodge_width)) +
    geom_linerange(aes(xmin = ll, xmax = hh, color = overall), linewidth = 0.6, position = position_dodge(dodge_width)) +
    geom_point(data = data_pts, aes(x = m, y = parameter)) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "#474A51") +
    theme_custom()   #return(combined)
}

plot_2lvl <- function(x, data) {
  x |>
  ggplot(aes(y = mu_study, x = i, group = i)) +
  stat_halfeye(
    aes(fill = after_stat(level)),
    shape = 124,
    slab_linewidth = 0.01,
    interval_color = "gray80",
    slab_color = "gray80",
    point_colour = "#474A51",
    alpha = 1,
    size = 0.01
  ) +
  scale_fill_manual(values = c("gray75", "#009E73")) +
  geom_point(
    data = data,
    colour = "#131516",
    size = 2,
    aes(x = i, y = yi, group = i)
  ) +
  geom_hline(yintercept = 0,
             color = "#474A51",
             linetype = "dashed") +
  theme_minimal() +
  coord_flip() +
  labs(subtitle = "Meta Analysis ", x = NULL, y = NULL) +
  theme(
    text = element_text(family = "Spectral"),
    axis.title.y = element_blank(),
    plot.margin = unit(c(1, 1, 1, -20), "cm"),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = c("Overall", data$author_year)[as.integer(levels(data$i)) + 1]) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
}

# p <- compare_posteriors(mod_two_level_cp_out, mod_three_level_cp_out, data_pts = plot_dat, vars = c("mu", "mu_study"), dodge_width = 1.8) 
# group_colors <- c(data = "black" , mod_three_level_cp_out_yes ="#FF6347",mod_three_level_cp_out_no ="#FF6347", mod_two_level_cp_out_yes ="#0072B2",
#                   mod_two_level_cp_out_no ="#0072B2")
# 
# group_alpha <- c(data = 1 , mod_three_level_cp_out = 0.7, mod_two_level_cp_out = 0.8)
# group_sizes <- c(data = 3 , mod_three_level_cp_out = 0.5, mod_two_level_cp_out = .5)
# showtext_opts(dpi = 500)
# p + scale_alpha_manual(values = group_alpha, guide = "none") +
#   scale_size_manual(values = group_sizes, guide = "none") +
#   labs(subtitle = "Meta Analysis ",
#        x = NULL, 
#        y = NULL) +
#   annotate("text", x= 0.4, y=1, label= "Overall", family = "Spectral", size = 2.1) +
#   # scale_y_discrete(labels=function(l) parse(text=gsub("^mu$", "Overall", l))) + 
#   # scale_y_discrete(labels=function(l) parse(text=gsub("^.*\\[([0-9]+)\\].*$", "\\1", l))) +
#   scale_color_manual(values = group_colors, 
#                      labels = c(Data = 'Data', 
#                                 mod_three_level_cp_out = 'Three Level Model', 
#                                 mod_two_level_cp_out = 'Two Level Model')) +
#   theme(legend.position = "none") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank()) 
# ggsave("meta.png", width = 19, height = 14, units = "cm")
meta_dat[, "i"]
