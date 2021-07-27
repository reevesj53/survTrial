# Set globalVariables to minimize check notes

# trial_sim
utils::globalVariables(c("death", "eventt", "subjid", "totalt"))

# calc_km
utils::globalVariables(c('kmfit', 'rx', 'data', 'km', 'time',
                         'surv', 'median', 'n', 'int_low', 'int_high'))
