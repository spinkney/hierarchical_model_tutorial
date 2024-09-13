library(cmdstanr)
library(rstanarm)
library(tidybayes)
library(ggplot2)
library(bayesplot)
library(tidybayes)
library(ggdist)
library(broom.mixed)
library(emmeans)

hdp <- read.csv("https://stats.idre.ucla.edu/stat/data/hdp.csv")
hdp <- within(hdp, {
  Married <- factor(Married, levels = 0:1, labels = c("no", "yes"))
  DID <- factor(DID)
  HID <- factor(HID)
  
  CancerStage <- factor(CancerStage)
})

remission_by_hospital <- hdp |>
  group_by(HID) |>
  summarize(remissions = sum(remission),
            total_patients = n(),
            lawsuits = sum(unique(Lawsuits)),
            num_doctors = n_distinct(DID)) |>
  filter(remissions != 0, total_patients != 0) |>
  mutate(empirical_prob = remissions / total_patients,
         num_doctors = factor(num_doctors, ordered = TRUE)
     #    lawsuits_per_doc = lawsuits / num_doctors,
      #   patients_per_doc = total_patients / num_doctors
     )

remission_by_hospital |>
 reframe(unique(num_doctors))

total_avg <- remission_by_hospital |>
  summarize(total_avg = sum(remissions) / sum(total_patients))

N <- nrow(remission_by_hospital)
# complete pooling
# remission ~ alpha_all

rstanarm

mod_complete_pooling <- cmdstan_model("./stan/complete_pooling_binom.stan")
mod_complete_pooling_out <- mod_complete_pooling$sample(
  data = list(N = N,
              y = remission_by_hospital$remissions,
              K = remission_by_hospital$total_patients,
              J = 2,
              X = as.matrix(select(remission_by_hospital, lawsuits, num_doctors))),
  parallel_chains = 4
)

mod_out <- brm(bf(remissions | trials(total_patients) ~ 1 + num_doctors),
                      family = binomial(link = "logit"),
                      data = remission_by_hospital,
               backend = "cmdstanr",
               cores = 4,
               chains = 4)

inv_logit_scaled(1.05 + .24 * log(200) - 1.01 * log(500))
mod_out
emmeans(mod_out, ~ num_doctors, epred=TRUE, re_formula=NULL)
plot(conditional_effects(mod_out))
pooling <- 
  mod_complete_pooling_out$summary(variables = "phi",
  quantiles = ~ posterior::quantile2(., probs = c(0.1, 0.5, 0.9))) |>
  tibble(new_col = 1:N) |>
  select(last_col(), everything(), -variable)

# no pooling
# estimate separate regressions for each group
# remission ~ alpha_k
mod_no_pooling <- cmdstan_model("./stan/no_pooling_bern.stan")

mod_no_pooling_out <- mod_no_pooling$sample(
  data = list(N = N,
              y = remission_by_hospital$remissions,
              K = remission_by_hospital$total_patients),
  parallel_chains = 4
)
mod_no_pooling_out$summary("phi")

no_pooling <- 
  mod_no_pooling_out$summary(variables = "phi",
                            quantiles = ~ posterior::quantile2(., probs = c(0.1, 0.5, 0.9))) |>
  tibble(new_col = 1:N) |>
  select(last_col(), everything(), -variable)

# mod_no_pooling_out |>
# spread_draws(phi[N]) |>
#   filter(N <= 10) |>
#   ggplot(aes(y = as.factor(N), x = phi)) +
#   stat_halfeye(.width = c(.90, .5))


# partial pooling non-centered
# remissions ~ global_alpha + alpha_k
mod_partial_pooling <- cmdstan_model("./stan/partial_pooling_binom.stan")

mod_partial_pooling_out <- mod_partial_pooling$sample(
  data = list(N = N,
              y = remission_by_hospital$remissions,
              K = remission_by_hospital$total_patients),
  parallel_chains = 4,
  seed = 349823
)

mod_partial_pooling_out$summary("phi")

partial_pooling_nc <- 
  mod_partial_pooling_out$summary(variables = "phi",
                             quantiles = ~ posterior::quantile2(., probs = c(0.1, 0.5, 0.9))) |>
  tibble(new_col = 1:N) |>
  select(last_col(), everything(), -variable)

# partial pooling centered
# remissions ~ global_alpha + alpha_k
mod_partial_pooling_center <- cmdstan_model("./stan/partial_pooling_binom_centered.stan")

mod_partial_pooling_center_out <- mod_partial_pooling_center$sample(
  data = list(N = N,
              y = remission_by_hospital$remissions,
              K = remission_by_hospital$total_patients),
  parallel_chains = 4,
  seed = 349823
)

mod_partial_pooling_center_out$summary("phi")

partial_pooling_c <- 
  mod_partial_pooling_center_out$summary(variables = "phi",
                                  quantiles = ~ posterior::quantile2(., probs = c(0.1, 0.5, 0.9))) |>
  tibble(new_col = 1:N) |>
  select(last_col(), everything(), -variable)


# partial pooling pcp
# remissions ~ global_alpha + alpha_k
mod_partial_pooling_pcp <- cmdstan_model("./stan/partial_pooling_binom_pcp.stan")

# when w = 0 => centered
# when w = 1 => non-centered

# empirical variance binomial n * p * (1 - p) 
var_y <- (remission_by_hospital$doctors/remission_by_hospital$total_patients) * remission_by_hospital$empirical_prob *(1 - remission_by_hospital$empirical_prob)

w <- (1 + 1/sigma)^-1
w
mod_partial_pooling_pcp_out <- mod_partial_pooling_pcp$sample(
  data = list(N = N,
              y = remission_by_hospital$remissions,
              K = remission_by_hospital$total_patients, 
              w = w),
  parallel_chains = 4,
  seed = 349823
)

mod_partial_pooling_pcp_out$summary("phi")

partial_pooling_pcp <- 
  mod_partial_pooling_pcp_out$summary(variables = "phi",
                                  quantiles = ~ posterior::quantile2(., probs = c(0.1, 0.5, 0.9))) |>
  tibble(new_col = 1:N) |>
  select(last_col(), everything(), -variable)

## plot

library(ggplot2)
models <- c("complete pooling", "no pooling", "partial pooling partially noncentered", "partial pooling centered", "partial pooling nc")
estimates <- rbind(pooling, no_pooling, partial_pooling_pcp, partial_pooling_c, partial_pooling_nc) |>
  mutate(new_col = new_col / max(new_col))
colnames(estimates) <- c("x", "lb", "median", "ub")

plotdata <- data.frame(estimates, 
                       observed = rep(remission_by_hospital$empirical_prob, times = length(models)), 
                       model = rep(models, each = N), 
                       row.names = NULL) |>
  group_by(model) |>
  arrange(observed, .by_group = TRUE)

ggplot(plotdata, aes(x = observed, y = median, ymin = lb, ymax = ub)) +
  geom_hline(yintercept = total_avg$total_avg, color = "lightpink", linewidth = 0.75) +
  geom_abline(intercept = 0, slope = 1, color = "skyblue") + 
  geom_linerange(color = "gray60", size = 0.75) + 
  geom_point(size = 2.5, shape = 21, fill = "gray30", color = "white", stroke = 0.2) + 
  facet_grid(. ~ model) +
  coord_fixed() +
#  scale_x_discrete() +
  labs(x = "Observed Remissions / Number of Patients", y = "Predicted Remission") +
  ggtitle("Posterior Medians and 80% Intervals") +
  theme_tidybayes()

library(brms)

stancode(remission ~ (1 +  IL6 + CRP + CancerStage + LengthofStay + Experience || HID),
         family = bernoulli(link="logit"),
         data = hdp)

test <- make_standata(remission ~ Age + LengthofStay + FamilyHx + IL6 + CRP + CancerStage +  
                        Experience + (1 + LengthofStay | DID) + (1 || HID),
              family = bernoulli(link="logit"),
              data = hdp)

prior <- set_prior("multi_normal(")
fit_brms <- brm(remission ~ Age + LengthofStay + FamilyHx + IL6 + CRP + CancerStage +  
                     + (1 + LengthofStay | DID) + (1 || HID),
                family = bernoulli(link="logit"),
    data = hdp,
    iter = 400,
    backend = "cmdstanr",
    cores = 4,
    chains = 4,
    seed = 234523)

get_variables(fit_brms)

fit_brms |> 
  broom.mixed::tidy(effects = c("ran_pars"), conf.level = 0.8)
fit_brms |> 
  broom.mixed::tidy(effects = c("fixed"), conf.level = 0.8)

newdata <- data.frame(y = rnorm(500, mean = 10))
fit <- update(compiled, newdata = newdata, 
              chains = 1, iter = 1000, warmup = 500)
# dv ~ bv * wv + (1+wv|group)
# where bv is one or more between-group predictors and wv is one or more within-group predictors
mod_out <- stan_glmer(remission ~ Age + FamilyHx + LengthofStay +
                        IL6 + CRP + CancerStage + Experience +
                        (1 + LengthofStay | DID) + (1 | HID),
           family = binomial(link = "logit"),
           data = hdp,
           iter = 500,
           QR = TRUE,
           algorithm = "sampling",
           cores = 4
           )
summary(mod_out)
rstan::get_stanmodel(mod_out$stanfit)

stancode(remission ~ Age + FamilyHx + LengthofStay +
           IL6 + CRP + CancerStage + Experience +
           (1 + LengthofStay | DID) + (1 | HID),
         family = bernoulli(link="logit"),
         data = hdp)

# IL6: Continuous, interleukin 6, a proinflammatory cytokine commonly examined as an indicator of inflammation, cannot be lower than zero.
# CRP: Continuous, C-reactive protein, a protein in the blood also used as an indicator of inflammation. It is also impacted by BMI.

fit_brms <- brm(remission ~ Age + FamilyHx + LengthofStay +
                  IL6 + CRP + CancerStage + Experience +
                  (1 + LengthofStay | DID) + (1 | HID),
                family = bernoulli(link="logit"),
                data = hdp,
                iter = 400,
                backend = "cmdstanr",
                cores = 4,
                chains = 4,
                seed = 234523)
hdp


###

ss_quantile <- function(ss, N, q) {
  result <- rep(NA, N);
  for (n in 1:N) {
    result[n] <- sort(ss$theta[,n])[M * q];
  }
  return(result);
}

theta_10_pool <- ss_quantile(ss_pool, N, 0.1);
theta_50_pool <- ss_quantile(ss_pool, N, 0.5);
theta_90_pool <- ss_quantile(ss_pool, N, 0.9);

theta_10_no_pool <- ss_quantile(ss_no_pool, N, 0.1);
theta_50_no_pool <- ss_quantile(ss_no_pool, N, 0.5);
theta_90_no_pool <- ss_quantile(ss_no_pool, N, 0.9);

theta_10_hier <- ss_quantile(ss_hier, N, 0.1);
theta_50_hier <- ss_quantile(ss_hier, N, 0.5);
theta_90_hier <- ss_quantile(ss_hier, N, 0.9);

theta_10_hier_logit <- ss_quantile(ss_hier_logit, N, 0.1);
theta_50_hier_logit <- ss_quantile(ss_hier_logit, N, 0.5);
theta_90_hier_logit <- ss_quantile(ss_hier_logit, N, 0.9);

pop_mean <- sum(y) / sum(K);

df_plot2 <- data.frame(x = rep(y / K, 4),
                       y = c(theta_50_pool, theta_50_no_pool,
                             theta_50_hier, theta_50_hier_logit),
                       model = c(rep("complete pooling", N),
                                 rep("no pooling", N),
                                 rep("partial pooling", N),
                                 rep("partial pooling (log odds)", N)));

plot_bda3_fig_5_4 <-
  ggplot(df_plot2, aes(x=x, y=y)) +
  geom_hline(aes(yintercept=pop_mean), colour="lightpink") +
  geom_abline(intercept=0, slope=1, colour="skyblue") +
  facet_grid(. ~ model) +
  geom_errorbar(aes(ymin=c(theta_10_pool, theta_10_no_pool,
                           theta_10_hier, theta_10_hier_logit),
                    ymax=c(theta_90_pool, theta_90_no_pool,
                           theta_90_hier, theta_90_hier_logit)),
                width=0.005, colour="gray60") +
  geom_point(colour="gray30", size=0.75) +
  coord_fixed() +
  scale_x_continuous(breaks = c(0.2, 0.3, 0.4)) +
  xlab("observed rate, y[n] / K[n]") +
  ylab("chance of success, theta[n]") +
  ggtitle("Posterior Medians and 80% intervals\n(red line: population mean;  blue line: MLE)")
plot_bda3_fig_5_4;

