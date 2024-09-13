params <- function(k, p) {
  return( c(alpha = (4 * k * (1 - p)) / (2 - p)^2, beta = (2 - p) / (2 * p)))
}

int_end <- function(k, p, r) {
  out <- (r + 0.5 + (k * (1 - p))) / (2 - p)
  out
}

k <- 0.1
p <- 0.02
r <- 0

d <- params(k = k, p = p)
x <- int_end(k = k, p = p, r = r)

integrate(function(x) dgamma(x, shape = d[1], scale = d[2]), 0, upper = x + 0.12 )
pgamma(x + 0.305, shape = d[1], scale = d[2])
dnbinom(r, size = k, p = p)

X <- ((2 * p) * (2 * r + 1 + (2 * k *(1 - p)) / (2 - p))) / (2 - p)

library(expint)
gammainc(p + 1, k)


bs_lpdf <- function(x, mu, sigma, T) {
  a <- -log(2 * sigma) - 0.5 * log(2 * pi)
  b <- -0.5 * ( (T - x * mu) / (sigma * sqrt(x)) )^2
  c <- log(T + x * mu) - (3/2) * log(x)

  return( a + b + c)
  }
bs_lpdf2 <- function(x, beta, alpha) {
  a <- 0.5 * log(beta) - log(2 * alpha) - 0.5 * log(2 * pi)
  b <- -beta * (1 - x / beta)^2 / (2 * alpha^2 * x)
  c <- log1p(x/beta) - (3/2) * log(x)
  
  return(exp(a+b+c))
}


mu <- 1
sigma <- 5
T <- 15
bs_lpdf(x, mu, sigma, T)

alpha <- sigma / sqrt(mu * T)
beta <- T / mu

alpha <- 1.2
beta <- 6

bs_lpdf2(x, beta, alpha)


# Generate x values from 0.1 to 10, avoiding zero because of division by x in the density formula
x_values <- seq(0.1, 30, by = 0.1)
densities <- sapply(x_values, bs_lpdf2, beta = beta, alpha = alpha)

# Create a data frame for ggplot
data_for_plot <- data.frame(x = x_values, Density = densities)

# Plotting
ggplot(data_for_plot, aes(x = x, y = Density)) +
  geom_line() +  # Plot a line
  labs(title = "Density Function", x = "x values", y = "Density") +
  theme_minimal()  # Use a minimal theme for clarity

library(cmdstanr)
mod_bs <- cmdstan_model("./stan/fit_bs.stan")

mod_bs$sample(
  data = list(N = stan_data$L,
              y = stan_data$length),
  parallel_chains = 4
)

integrate(bs_lpdf2, lower = 0, upper=3.5, beta = 6, alpha = 1.2)

library(pracma)

beta <- 6
alpha <- 1.2

0.5 * alpha * sqrt(beta) / ((sqrt(beta / alpha^2)) )
alpha^2 * 0.5

a <- 0.5/(sqrt(beta / alpha^2)) * alpha * sqrt(beta)
b <- exp(1/alpha^2 - sqrt(1/(alpha^2 * beta)) * sqrt(beta/alpha^2)) * sqrt(1/(alpha^2 * beta))
inner_a <- beta * sqrt(1/(alpha^2 * beta)) + sqrt(beta/alpha^2) * erf( (x * sqrt(1 / (alpha^2 * beta)) - (sqrt(beta/alpha^2))) / sqrt(2 * x) ) 

inner_b <- (sqrt(beta /alpha^2) - beta * sqrt(1/(alpha^2 * beta)))
c <- -exp(2 * sqrt(1/(alpha^2 *beta))) * sqrt(beta/alpha^2) + 
  exp(2 * sqrt(1/(alpha^2 * beta) * sqrt(beta/alpha^2))) * erf(sqrt(beta/alpha^2) + x * sqrt(1/(alpha^2 * beta)) / (sqrt(2) * sqrt(x))) + 1

a * b * (inner_a  + inner_b * c)

a_times_b <- alpha/sqrt(beta)

inner_a2 <- sqrt(beta) / alpha + (sqrt(beta) /alpha) * erf( (x * 1 / (alpha * sqrt(beta)) - (sqrt(beta)/ alpha)) / sqrt(2 * x) )


0.5 + 0.5 * erf( (x * 1 / (alpha * sqrt(beta)) - (sqrt(beta)/ alpha)) / sqrt(2 * x) )

0.5 + 0.5 * erf( (x - beta) / (sqrt(2 * beta * x) * alpha))


