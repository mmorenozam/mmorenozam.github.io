library(ggplot2)
library(dplyr)
library(lmerTest)
library(emmeans)

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

p1 <- ggplot(design_data,
             aes(Treatment, y, color = Trial)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(
    aes(group = Treatment),
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "black"
  ) +
  theme_bw() +
  labs(y = "Observed mean of treatment per trial\n(M-TON/ha)",
       title = "Overall yield obtained per treatment/trial",
       caption = "Black segment represents the overall mean of the treatment across all trials")+
  theme(plot.caption = element_text(hjust = 0))

p1

# 5. Fit the Trial Series Model with lme4
# Syntax: (1 | Trial) treats Trial as random
# Syntax: (1 | Trial:Block) handles nesting explicitly
# Syntax: (1 | Treatment:Trial) models the GxE interaction
meta_model <- lmer(y ~ Treatment + (1 | Trial) + (1 | Trial:Block) + (1 | Treatment:Trial), 
                   data = design_data)

summary(meta_model)

emmeans(meta_model, specs = ~Treatment)

contrast(emmeans(meta_model, specs = ~ Treatment), method = "trt.vs.ctrl", ref = "Treatment_1", adjust = "dunnett")
