library(cmdstanr)
library(data.table)
library(tidybayes)
library(ggplot2)
library(showtext)
library(metafor)
source("./R/graph_code/compare_posterior.r")

# Results from 48 studies on the effectiveness of 
# school-based writing-to-learn interventions on academic achievement.

# id	      numeric	  study number
# author	  character	study author(s)
# year	    numeric	  publication year
# grade	    numeric	  grade level (1 = elementary; 2 = middle; 3 = high-school; 4 = college)
# length	  numeric	  treatment length (in weeks)
# minutes	  numeric	  minutes per assignment
# wic	      numeric	  writing tasks were completed in class (0 = no; 1 = yes)
# feedback	numeric	  feedback on writing was provided (0 = no; 1 = yes)
# info	    numeric	  writing contained informational components (0 = no; 1 = yes)
# pers	    numeric	  writing contained personal components (0 = no; 1 = yes)
# imag	    numeric	  writing contained imaginative components (0 = no; 1 = yes)
# meta	    numeric	  prompts for metacognitive reflection (0 = no; 1 = yes)
# subject	  character	subject matter
# ni	      numeric	  total sample size of the study
# yi	      numeric	  standardized mean difference
# vi	      numeric	  corresponding sampling variance

### data prep
meta_dat <- dat.bangertdrowns2004 |>
  mutate(
    id = forcats::as_factor(id),
    ni100 = ni / 100,
    author_year = paste0(author, " ", year),
    subject = factor(
    case_when(
      subject %in% c("Algebra", "Calculus", "Math", "Math in Science", "Statistics") ~ "Math",
      subject %in% c(
        "Chemistry",
        "Comp Science and Chemistry",
        "Biology",
        "Earth Science",
        "Natural Resources",
        "Nursing",
        "Science"
      ) ~ "Sci",
      TRUE ~ "Soc"  # Default case if none of the above matches
    )
  ))
res1 <- rma(yi, vi, data=meta_dat, test="knha")

summary_dat <-
  meta_dat |>
  group_by(subject) |>
  summarize(Mean = mean(yi),
            Std_Dev = sd(yi),
            Min = min(yi),
            Max = max(yi))

summary_dat

stan_data_fun <- \(x) {
  x_dt <- data.table::data.table(x)

    list(
      N = x_dt[, .N],
      y = x_dt[, yi],
      study_sigma = x_dt[, sqrt(vi)],
      D = x_dt[, length(unique(subject))],
      G = x_dt[, length(unique(grade))],
      subject_idx =
        x_dt[, subject_grp := .GRP, subject][, subject_grp],
      grade_idx =
        x_dt[, grade_grp := .GRP, grade][, grade_grp],
      L = x_dt[!is.na(length), .N],
      length = x_dt[!is.na(length), length],
      missing_length_idx = x_dt[, which(is.na(length))],
      complete_length_idx = x_dt[, which(!is.na(length))],
      ni = x_dt[, ni100]
      )
}

stan_data <- stan_data_fun(meta_dat)

meta_dat <- meta_dat |>
  mutate(subject_id = stan_data$subject_idx)

# stan_data <- stan_data_fun(meta_dat)

####################################################################################
####################################################################################

### stan model 1: 2 level centered parameterization
## question: the code only estimates one 'sigma' instead of N
##  what is the interpretation of sigma mean here and 
##  what would it mean if we had N sigma's?

mod_two_level_cp <- cmdstan_model("./stan/2_meta_two_level_cp.stan")
cat(mod_two_level_cp$code(), sep = "\n")

mod_two_level_cp_out <- mod_two_level_cp$sample(
  data = stan_data,
  parallel_chains = 4,
  seed = 1231290
)
mod_two_level_cp_out$summary(c("intercept", "sigma"))

res_2lvl <- rma(yi, vi, data= meta_dat, test="knha")
res_2lvl # look at tau value for sigma and estimate value for mu

mu_draws <- mod_two_level_cp_out |>
       spread_draws(intercept, mu_study[i]) |>
       mutate(mu_study = intercept,
              i = 0) |>
  select(i, mu_study, .chain, .iteration, .draw)

### Plot 
meta_dat <- meta_dat |>
  mutate(i = as.factor(id)) |>
  mutate(i = forcats::fct_reorder(i, yi)) |>
  mutate(i = forcats::fct_expand(i, "0", after = 0))

plot_2lvl_dat <- mod_two_level_cp_out |>
  spread_draws(mu_study[i]) |>
  bind_rows(mu_draws) 

plot_2lvl_dat$i <- factor(plot_2lvl_dat$i, levels(meta_dat$i))

plot_2lvl(plot_2lvl_dat, meta_dat)

####################################################################################
####################################################################################

### stan model 2: Add study size regressor
## question: What other regressors/predictors would be good to add?
##   Are there any issues with adding these predictors?
mod_two_level_ncp_reg <- cmdstan_model("./stan/2_meta_two_level_ncp_reg.stan")

stan_data$X <- stan_data$ni 

mod_two_level_ncp_reg_out <- mod_two_level_ncp_reg$sample(
  data = stan_data,
  parallel_chains = 4,
  seed = 1231290
)
mod_two_level_ncp_reg_out$summary(c("intercept", "beta_mu",  "beta_sigma", "intercept_sigma", "mu_est"))

# compare to frequentist
res_2lvl_loc_scale <- rma(yi, vi, mods = ~ ni100 , scale = ~ ni100, data= meta_dat, test="knha")
summary(res_2lvl_loc_scale)

mod_two_level_cp_reg_out$summary(c("intercept_est", 
                                   "beta_mu", 
                                   "intercept_sigma_est",
                                   "beta_sigma"))

####################################################################################
####################################################################################

### stan model 3: Three-level meta-analysis with subject
## Three-level meta-analysis with district
## question: Why didn't we put a regessor on the 'subject'-level std.dev?

mod_three_level_ncp_reg <- cmdstan_model("./stan/2_meta_three_level_ncp_reg.stan")

mod_three_level_ncp_reg_out <- mod_three_level_ncp_reg$sample(
  data = stan_data,
  parallel_chains = 4,
  adapt_delta = 0.95,
  seed = 1231290
)

mod_three_level_ncp_reg_out$summary(c("intercept_sigma_study","beta_sigma",
                                  "sigma_subject", 
                                 "intercept", "beta_mu"))

mod_three_level_ncp_reg_out$summary("mu_subject")

mod_three_level_ncp_reg$summary(c("mu_est",
                                   "intercept_est", 
                                   "beta_mu", 
                                   "intercept_sigma_est",
                                   "beta_sigma"))
res_3lvl <- rma.mv(yi, 
                   vi, 
                   random = ~ 1 | subject/id,
                   data = meta_dat)

round(res_3lvl$sigma2[1] / sum(res_3lvl$sigma2), 3)

rma.mv(yi, 
       vi, 
       random = ~ id | subject,
       data = meta_dat)
