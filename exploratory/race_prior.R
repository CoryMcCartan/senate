library(tidyverse)
library(rstanarm)

############################
# LOAD AND PREP DATA
############################

control = read_csv("data/party_control.csv") %>%
    select(year=elect_year, senate=dem_senate, pres=dem_pres) %>%
    mutate(senate = 2*senate - 1,
           pres = 2*pres - 1) %>%
    bind_rows(tibble(year=2020, senate=-1, pres=-1))
pres = read_csv("data/state_data_combined.csv") %>%
    mutate(dem = dem/(dem+gop)) %>%
    select(year, abbr, pres_dem=dem, pres_natl=natl_dem) %>%
    mutate(pres_dem = qlogis(pres_dem),
           pres_natl = qlogis(pres_natl))
generic = read_csv("data/generic_ballot.csv") %>%
    select(year, generic = logit_pop)

d_2020 = read_csv("data/2020.csv") %>%
    mutate(midterm = as.logical(midterm))

raw = read_csv("data/returns.csv")
d = raw %>%
    mutate(party = if_else(party == "democratic-farmer-labor", "democrat", party)) %>%
    filter(party %in% c("democrat", "republican"), !special) %>%
    mutate(candidate = str_remove_all(candidate, "Jr") %>%
               str_remove_all("Sr") %>%
               str_remove_all("III") %>%
               str_to_lower("") %>%
               str_remove_all("[^a-z ]") %>%
               str_replace_all(" \\w ", " ")) %>%
    select(year, abbr=state_po, party, candidate, votes=candidatevotes) %>%
    group_by(year, abbr, party) %>%
    summarize(candidate = candidate[which.max(votes)],
              votes = sum(votes)) %>%
    mutate(cycle = 1 + (year %% 6)/2) %>%
    filter(n() == 2) %>%
    group_by(cycle, abbr, party) %>%
    mutate(incumbent = candidate == lag(candidate, default="")) %>%
    ungroup %>%
    select(year, cycle, abbr, party, votes, incumbent) %>%
    mutate(party = str_sub(party, 1, 1)) %>%
    pivot_wider(names_from=party, values_from=c(votes, incumbent)) %>%
    mutate(dem = log(votes_d / votes_r),
           inc = incumbent_d - incumbent_r,
           midterm = year %% 4 == 2) %>%
    select(-starts_with("votes_"), -starts_with("incumbent_")) %>%
    bind_rows(d_2020) %>%
    left_join(control, by="year") %>%
    full_join(pres, by=c("abbr", "year")) %>%
    left_join(generic, by="year") %>%
    arrange(abbr, year) %>%
    fill(pres_dem, pres_natl, .direction="down") %>%
    arrange(year, abbr) %>%
    group_by(abbr) %>%
    mutate(pres_avg = 0.5*lag(pres_dem - pres_natl, 1) + 0.5*lag(pres_dem - pres_natl, 3),
           pres_level = lag(pres_dem)) %>%
    group_by(cycle, abbr) %>%
    mutate(level = dem - generic,
           ref_level = if_else(inc != 0, lag(dem), pres_level) - lag(generic)) %>%
    filter(year >= 1982)

# manually fill in special elections and fix errors
d$inc[d$year==2014 & d$abbr=="KY"] = -1
d$inc[d$year==2020 & d$abbr=="KY"] = -1
d$ref_level[d$year==2020 & d$abbr=="AZ"] = -0.1733
d$ref_level[d$year==2020 & d$abbr=="GA-S"] = -0.1523
d$pres_avg[d$year==2020 & d$abbr=="AZ"] = -0.2505
d$pres_avg[d$year==2020 & d$abbr=="GA-S"] = -0.331

write_csv(d, "data/senate_combined.csv")

############################
# FIT MODEL
############################

d = filter(d, year != 2020) %>%
    drop_na

m = stan_lmer(level ~ ref_level + pres_avg +  senate*inc + pres*inc*midterm + (1|year),
             data=d, prior=cauchy(), chains=2)

