# Set the following variables to the desired values:
OUTPUT_DIR <- 'results'  
STATES <- c(
    'Tamil Nadu', 'Andhra Pradesh', 'Odisha', 'West Bengal'
)

# Load libraries and local scripts
library(tidyverse)
library(INLA)

inla.setOption(scale.model.default = TRUE)

source("R/inla_formula.R")
source("R/train.R")
source("R/models.R")
source("R/model_average.R")
source("R/postprocessing.R")

cases <- read_csv("data/data.csv") %>%
    filter(state %in% STATES) %>% 
    prep_train_data

# Set up a data frame with models to fit.
train_grid <- list(NULL, NULL)
train_grid <- cases %>%
    group_by(state) %>%
    group_modify(function(d, g) {
        model_fmls()
    }) %>%
    ungroup %>%
    arrange(state, model_name)


# At this point train_grid contains a row for each state, and model
# to be trained. Next, train the models. 
fits <- train_grid %>%
    train_many(data = cases, num_posterior_draws = 1000)

# Combine posterior draws to obtain draws from the model average
results <- combine_draws(
    fits, num_posterior_draws = 16000, with_replacement = FALSE)

results <- results %>%
    calculate_excess_cases

# Summarise results
result_summaries <- summarise_results(results)
# samples_summary <- extract_all_samples(results)

# Save results
dir.create(OUTPUT_DIR, showWarnings = FALSE)
write_csv(result_summaries, file.path(OUTPUT_DIR, 'result_summaries.csv'))
# write_csv(samples_summary, file.path(OUTPUT_DIR, 'samples.csv'))