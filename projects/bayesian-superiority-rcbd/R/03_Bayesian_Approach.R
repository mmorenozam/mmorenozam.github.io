library(brms)

# 1. We set the general parameters for the simulation

set.seed(76)

n_trials     <- 3
n_blocks     <- 3
n_treatments <- 4

# Overall mean
mu <- 100

# Fixed treatment effects
trt_effects <- c(
  Treatment_1 = 0,
  Treatment_2 = 0,
  Treatment_3 = 3,
  Treatment_4 = 5
)

# Standard deviations of random effects
trial_sd      <- 10
block_sd      <- 2
gxe_sd        <- 3
residual_sd   <- 4


# 2. Experimental layout 

design_data <- expand.grid(
  Trial     = factor(1:n_trials),
  Block     = factor(1:n_blocks),
  Treatment = factor(names(trt_effects), levels = names(trt_effects))
)

# Block nested within Trial
design_data$Trial_Block <-
  interaction(design_data$Trial,
              design_data$Block,
              sep = ":")

# Treatment × Trial interaction
design_data$Treatment_Trial <-
  interaction(design_data$Treatment,
              design_data$Trial,
              sep = ":")


# 3. Simulate random effects 

# Trial effects
trial_effects <-
  rnorm(nlevels(design_data$Trial),
        mean = 0,
        sd = trial_sd)

names(trial_effects) <- levels(design_data$Trial)

# Trial:Block effects
block_effects <-
  rnorm(nlevels(design_data$Trial_Block),
        mean = 0,
        sd = block_sd)

names(block_effects) <- levels(design_data$Trial_Block)

# Treatment × Trial interaction effects
gxe_effects <-
  rnorm(nlevels(design_data$Treatment_Trial),
        mean = 0,
        sd = gxe_sd)

names(gxe_effects) <- levels(design_data$Treatment_Trial)

# Residual errors
residual_errors <-
  rnorm(nrow(design_data),
        mean = 0,
        sd = residual_sd)


# 4. Generate response 

design_data$y <-
  mu +
  trt_effects[design_data$Treatment] +
  trial_effects[design_data$Trial] +
  block_effects[design_data$Trial_Block] +
  gxe_effects[design_data$Treatment_Trial] +
  residual_errors

priors <- c(
  prior(normal(100, 20), class = "Intercept"),
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 0, 10), class = "sd"),
  prior(student_t(3, 0, 10), class = "sigma")
)

# bayes_mod <- brm(
#   y ~ Treatment + (1 | Trial) + (1 | Trial:Block) + (1 | Treatment:Trial),
#   data = design_data,
#   prior = priors,
#   chains = 4, iter = 6000, warmup = 3000, cores = 4,
#   seed = 1410,
#   sample_prior = "yes",
#   control = list(adapt_delta = 0.99, max_treedepth = 12),
#   refresh = 1000
# )

# saveRDS(bayes_mod, file = here::here("projects/bayesian-superiority-rcbd/bayes_model/bayes_mod.rds"))

bayes_mod <- readRDS(here::here("projects/bayesian-superiority-rcbd/bayes_model/bayes_mod.rds"))

bayes_mod2 <- update(
  bayes_mod,
  newdata = design_data,
  prior = priors,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  chains = 4, iter = 8000, warmup = 4000, cores = 4,
  recompile = FALSE
)

saveRDS(bayes_mod2, file = here::here("projects/bayesian-superiority-rcbd/bayes_model/bayes_mod.rds"))

draws <- as_draws_df(bayes_mod2)
poc_threshold <- 3

summarize_treatment <- function(delta, label, prior_vector) {
  prior_odds <- (sum(prior_vector > poc_threshold)/length(prior_vector))/(sum(prior_vector < poc_threshold)/length(prior_vector))
  ha <- hypothesis(bayes_mod, paste0("Treatment", label, " > ", poc_threshold))
  posterior_odds <- ha$hypothesis$Evid.Ratio
  data.frame(
    Treatment = label,
    `P(delta > 3%, meets POC)` = round(mean(delta >= poc_threshold), 3),
    `Bayes Factor (BF10)` = round(posterior_odds/prior_odds, 2),
    check.names = FALSE
  )
}

bf_table <- bind_rows(
  summarize_treatment(draws$b_TreatmentTreatment_2, "Treatment_2", draws$prior_b),
  summarize_treatment(draws$b_TreatmentTreatment_3, "Treatment_3", draws$prior_b),
  summarize_treatment(draws$b_TreatmentTreatment_4, "Treatment_4", draws$prior_b)
)

bf_table
