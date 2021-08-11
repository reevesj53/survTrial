# Set globalVariables to minimize check notes

# trial_sim
utils::globalVariables(c("death", "eventt", "subjid", "totalt", "totalmax"))

# calc_km
utils::globalVariables(c('rx', 'data', 'kmfit', 'median', 'sim_low', 'sim_high'))

# calc_hr
utils::globalVariables(c('estimate', 'HR'))

# plot_sim
utils::globalVariables(c('description', 'KM_median', 'treatment'))

