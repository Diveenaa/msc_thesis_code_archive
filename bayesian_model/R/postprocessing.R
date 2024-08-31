calculate_excess_cases <- function(result_grid) {
    monthly_excess <- result_grid %>%
        group_by(state) %>%
        group_modify(function(d, g) {
            oos_mask <- d$data[[1]]$out_of_sample

            excess_case_draws <-
                d$data[[1]][oos_mask, ]$cases - d$case_draws[[1]][oos_mask, ]
            excess_rate_draws <-
                excess_case_draws / d$data[[1]][oos_mask, ]$population
            relative_increase_draws <-
                excess_case_draws / d$case_draws[[1]][oos_mask, ]
            bind_cols(d, tibble(
                excess_case_draws = list(excess_case_draws),
                excess_rate_draws = list(excess_rate_draws),
                relative_increase_draws = list(relative_increase_draws)
            ))
        }) %>%
        ungroup
    all_month_excess <- monthly_excess %>%
        group_by(state) %>%
        group_modify(function(d, g) {
            oos_mask <- d$data[[1]]$out_of_sample
            total_excess_cases <- colSums(d$excess_case_draws[[1]])
            total_excess_rates <-
                total_excess_cases / mean(d$data[[1]][oos_mask, ]$population)
            total_relative_increase <-
                total_excess_cases / colSums(d$case_draws[[1]][oos_mask, ])
            bind_cols(d, tibble(
                all_month_excess_case_draws = list(total_excess_cases),
                all_month_excess_rate_draws = list(total_excess_rates),
                all_month_relative_increase_draws = list(total_relative_increase)
            ))
        }) %>%
        ungroup
    all_month_excess
}


summarise_results <- function(
    result_grid, quantiles = seq(0.025, 0.975, .025)
) {
    result_grid %>%
        group_by(state) %>%
        group_modify(function(d, g) {
            oos_mask <- d$data[[1]]$out_of_sample
            months <- as.character(d$data[[1]][oos_mask, ]$month)

            draws_summary <- function(draws) {
                if (is.matrix(draws)) {
                    qs <- matrixStats::rowQuantiles(draws, probs = quantiles)
                } else {
                    stopifnot(is.vector(draws))
                    qs <- quantile(draws, probs = quantiles) %>% t
                }
                qs %>%
                    as_tibble %>%
                    set_names(paste0('q', sprintf(
                        '%03.0f', parse_number(sub('%$', '', names(.))) * 10)))
            }

            cases_summary <- bind_rows(
                d$case_draws[[1]][oos_mask, ] %>%
                    draws_summary %>%
                    mutate(month = months),
                colSums(d$case_draws[[1]][oos_mask, ]) %>%
                    draws_summary %>%
                    mutate(month = 'all')
            ) %>% mutate(quantity = 'Predicted cases')

            excess_cases_summary <- bind_rows(
                d$excess_case_draws[[1]] %>%
                    draws_summary %>%
                    mutate(month = months),
                d$all_month_excess_case_draws[[1]] %>%
                    draws_summary %>%
                    mutate(month = 'all')
            ) %>% mutate(quantity = 'Excess cases')

            excess_rates_summary <- bind_rows(
                d$excess_rate_draws[[1]] %>%
                    draws_summary %>%
                    mutate(month = months),
                d$all_month_excess_rate_draws[[1]] %>%
                    draws_summary %>%
                    mutate(month = 'all')
            ) %>%
                mutate_at(vars(matches('^q\\d{3}$')), ~.x * 1e5) %>%
                mutate(quantity = 'Excess cases per 100,000')

            relative_increase_summary <- bind_rows(
                d$relative_increase_draws[[1]] %>%
                    draws_summary %>%
                    mutate(month = months),
                d$all_month_relative_increase_draws[[1]] %>%
                    draws_summary %>%
                    mutate(month = 'all')
            ) %>% mutate(quantity = 'Relative change in cases')

            bind_rows(
                cases_summary, excess_cases_summary, excess_rates_summary,
                relative_increase_summary
            )
        }) %>%
        ungroup %>%
        select(quantity, state, month, month, matches('^q\\d{3}$'))
}

extract_all_samples <- function(result_grid) {
  result_grid %>%
    group_by(state) %>%
    group_modify(function(d, g) {
      oos_mask <- d$data[[1]]$out_of_sample
      months <- as.character(d$data[[1]][oos_mask, ]$month)
      
      # Extract raw posterior draws for each month
      cases_samples <- bind_rows(
        as_tibble(d$case_draws[[1]][oos_mask, ]) %>%
          mutate(month = months),
        as_tibble(matrix(colSums(d$case_draws[[1]][oos_mask, ]), nrow = 1)) %>%
          mutate(month = 'all')
      ) %>%
        mutate(quantity = 'Predicted cases')
      
      # Extract raw posterior draws for excess cases
      excess_cases_samples <- bind_rows(
        as_tibble(d$excess_case_draws[[1]]) %>%
          mutate(month = months),
        as_tibble(matrix(d$all_month_excess_case_draws[[1]], nrow = 1)) %>%
          mutate(month = 'all')
      ) %>%
        mutate(quantity = 'Excess cases')
      
      # Extract raw posterior draws for excess rates
      excess_rates_samples <- bind_rows(
        as_tibble(d$excess_rate_draws[[1]]) %>%
          mutate(month = months),
        as_tibble(matrix(d$all_month_excess_rate_draws[[1]], nrow = 1)) %>%
          mutate(month = 'all')
      ) %>%
        mutate(quantity = 'Excess cases per 100,000')
      
      # Extract raw posterior draws for relative increases
      relative_increase_samples <- bind_rows(
        as_tibble(d$relative_increase_draws[[1]]) %>%
          mutate(month = months),
        as_tibble(matrix(d$all_month_relative_increase_draws[[1]], nrow = 1)) %>%
          mutate(month = 'all')
      ) %>%
        mutate(quantity = 'Relative change in cases')
      
      # Combine all the samples
      bind_rows(
        cases_samples, excess_cases_samples, excess_rates_samples,
        relative_increase_samples
      )
    }) %>%
    ungroup %>%
    select(quantity, state, month, everything())  # Return all posterior draws
}