library(lme4)
library(lmerTest)
library(dplyr)
library(emmeans)

simulate_and_test <- function(n_trials, 
                              n_blocks = 3, 
                              n_treatments = 4,
                              mu = 100, 
                              trt_effects = c(
                                Treatment_1 = 0,
                                Treatment_2 = 0,
                                Treatment_3 = 3,
                                Treatment_4 = 5
                              ),
                              trial_sd = 10, 
                              block_sd = 2.5,
                              gxe_sd = 2, 
                              residual_sd = 5) {
  
  design_data <- expand.grid(
    Trial     = factor(1:n_trials),
    Block     = factor(1:n_blocks),
    Treatment = factor(names(trt_effects), levels = names(trt_effects))
  )
  
  design_data$Trial_Block <-
    interaction(design_data$Trial,
                design_data$Block,
                sep = ":")
  
  design_data$Treatment_Trial <-
    interaction(design_data$Treatment,
                design_data$Trial,
                sep = ":")
  
  trial_effects <-
    rnorm(nlevels(design_data$Trial),
          mean = 0,
          sd = trial_sd)
  
  names(trial_effects) <- levels(design_data$Trial)
  
  block_effects <-
    rnorm(nlevels(design_data$Trial_Block),
          mean = 0,
          sd = block_sd)
  
  names(block_effects) <- levels(design_data$Trial_Block)
  
  gxe_effects <-
    rnorm(nlevels(design_data$Treatment_Trial),
          mean = 0,
          sd = gxe_sd)
  
  names(gxe_effects) <- levels(design_data$Treatment_Trial)
  
  residual_errors <-
    rnorm(nrow(design_data),
          mean = 0,
          sd = residual_sd)
  
  design_data$y <-
    mu +
    trt_effects[design_data$Treatment] +
    trial_effects[design_data$Trial] +
    block_effects[design_data$Trial_Block] +
    gxe_effects[design_data$Treatment_Trial] +
    residual_errors
  
  tryCatch({
    mod <- suppressMessages(suppressWarnings(
      lmer(y ~ Treatment + (1|Trial) + (1|Trial:Block) + (1|Treatment:Trial),
           data = design_data)
    ))
    em <- emmeans(mod, specs = ~ Treatment)
    ct <- as.data.frame(contrast(em, method = "trt.vs.ctrl", ref = "Treatment_1", adjust = "dunnett"))
    ct$p.value[ct$contrast == "Treatment_4 - Treatment_1"]
  }, error = function(e) NA, warning = function(w) NA)
}

set.seed(624)
n_trials_grid <- c(3, 5, 10, 20, 30, 50, 70, 90)
n_reps <- 250

power_results <- lapply(n_trials_grid, function(nt) {
  pvals <- unlist(replicate(n_reps, simulate_and_test(n_trials = nt)))
  data.frame(n_trials = nt, power = mean(pvals < 0.05, na.rm = TRUE))
})
power_df <- do.call(rbind, power_results)

write.csv(power_df, here::here("projects/bayesian-superiority-rcbd/data/power_df.csv"), row.names = F)



