# Load required libraries
library(ggplot2)
library(showtext)
library(ggbrace)

# Add Google font 'Spectral'
font_add_google(name = "Spectral", family = "Spectral")
showtext_auto()

# Create data for the curve (-x^3) and the tangent line (-3x^2)
curve_data <- data.frame(
  x = seq(0.05, 0.42, length.out = 200),
  y = -(seq(0.05, 0.42, length.out = 200))^3 # Curve -x^3
)

# Derivative at x = 0.2
x_intersect = 0.2
slope = -3 * x_intersect^2 # Slope of the tangent at x = 0.2
intercept = -x_intersect^3 - slope * x_intersect # y = -x^3 at x = 0.2, line: y = slope * x + intercept

# Tangent line data at the point x = 0.2
tangent_data <- data.frame(
  x = seq(0.05, 0.42, length.out = 200),
  y = slope * seq(0.05, 0.42, length.out = 200) + intercept # Tangent line at x = 0.2
)

# Intersection point (current)
current_point <- data.frame(
  x = x_intersect,
  y = -x_intersect^3,
  label = "Particle"
)

# Additional point on the derivative line at x = 0.3
additional_point <- data.frame(
  x = c(0.3, 0.4),
  y = slope * c(0.3, 0.4) + intercept
)

# Plot
showtext_opts(dpi = 800)
text_size = 3
line_size <- 1.5

p <- ggplot() +
  geom_line(data = curve_data, aes(x = x, y = y), color = "#131516", size = line_size) + # Curve -x^3
  geom_line(data = tangent_data, aes(x = x, y = y), color = "#0072B2", size = line_size) + # Tangent line -3x^2 at x = 0.2
  geom_point(data = current_point, aes(x = x, y = y), shape = 16, color = "#009E73", size = 5) + # Intersection point
  geom_text(data = current_point, aes(x = x, y = y, label = label), vjust = -1, hjust = -0.1, family = "Spectral", size = text_size) +
  geom_point(data = additional_point, aes(x = x, y = y), shape = 16, color = "#58212b", size = 3) +
  theme_minimal(base_family = "Spectral") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(legend.position = "none")  +
  labs(x = NULL, y = NULL) 
p <- p + stat_brace(data = data.frame(x = c(0.3,0.4), y = slope * c(0.3, 0.355) + intercept, text = "Small delta"), aes(x = x, y = y), rotate = 270, size = 0.75) +
  stat_bracetext(data = data.frame(x = c(0.3,0.4), y = slope * c(0.3, 0.355) + intercept, text = "Small delta"), aes(x = x, y = y, label = text), family = "Spectral", rotate = 270, size = text_size) +
  stat_brace(data = data.frame(x = c(0.4,0.5), y = slope * c(0.4, 0.66) + intercept, text = "Large delta"), aes(x = x, y = y), rotate = 270, size = 0.75) +
  stat_bracetext(data = data.frame(x = c(0.4,0.5), y = slope * c(0.4, 0.66) + intercept, text = "Large delta"), aes(x = x, y = y, label = text), family = "Spectral", rotate = 270, , size = text_size) 
p  
# Save plot as SVG
ggsave("log_density_plot.png", width = 8, height = 7)
