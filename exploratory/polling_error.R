library(tidyverse)
library(rstanarm)

# senate races
polls = read_csv("data/raw_polls.csv") %>%
    filter(type_detail == "Sen-G") %>%
    mutate(wgt = 100^2 * samplesize / (cand1_pct*(100-cand1_pct)),
           dem_poll = cand1_pct / (cand1_pct + cand2_pct),
           dem_act = cand1_actual / (cand1_actual + cand2_actual)) %>%
    select(year, state=location, samplesize, wgt, dem_poll, dem_act) %>%
    group_by(year, state) %>%
    summarize(dem_poll = weighted.mean(dem_poll, wgt),
              dem_act = mean(dem_act),
              error = qlogis(dem_poll) - qlogis(dem_act))

m = stan_lmer(error ~ (1|year), data=polls, chains=1)
m_sds = as.data.frame(VarCorr(m, sigma=1))$sdcor

errors = list(prior_natl_poll_error = m_sds[1],
              prior_race_poll_error = m_sds[2])

# generic ballot
polls = read_csv("data/raw_polls.csv") %>%
    filter(type_detail == "House-G", location=="US") %>%
    mutate(wgt = 100^2 * samplesize / (cand1_pct*(100-cand1_pct)),
           dem_poll = cand1_pct / (cand1_pct + cand2_pct),
           dem_act = cand1_actual / (cand1_actual + cand2_actual)) %>%
    select(year, samplesize, wgt, dem_poll, dem_act) %>%
    group_by(year) %>%
    summarize(dem_poll = weighted.mean(dem_poll, wgt),
              dem_act = mean(dem_act),
              error = qlogis(dem_poll) - qlogis(dem_act))

errors$prior_generic_poll_bias = mean(polls$error) / 2
errors$prior_generic_poll_error = sd(polls$error)

write_rds(errors, "output/polling_error.rdata")
