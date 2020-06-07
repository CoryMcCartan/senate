#!/usr/bin/env Rscript

#####################################
#   U.S. SENATE MODEL 2020          #
#   CORY McCARTAN                   #
#   (c) 2020                        #
#####################################

library(optparse)

###
### Parse options
###

option_list = list(
    make_option("--dry", action="store_true", default=F,
                help="Dry run, results not saved."),
    make_option("--date", type="character", default=as.character(Sys.Date()),
                help="The date to estimate from."),
    make_option("--iter", type="integer", default=1500,
                help="Number of MCMC iterations for voter intent estimation,
                      not including warmup iterations."),
    make_option("--chains", type="integer", default=2,
                help="Number of MCMC chains for voter intent estimation."),
    make_option("--recompile", action="store_true", default=F,
                help="Force recompile of Stan models."),
    make_option("--model_dir", type="character", default="stan",
                help="The directory in which the models are stored"),
    make_option("--output_file", type="character", default="docs/estimate.json",
                help="The file to save estimates to."),
    make_option("--sims_file", type="character", default="docs/sims.json",
                help="The file to save race simulations to."),
    make_option("--history_file", type="character", default="docs/history.csv",
                help="The file to save model history to.")
)
opt = parse_args(OptionParser(option_list=option_list,
                              description="Forecast the 2020 U.S. Senate elections."))
opt$iter = 3*ceiling(opt$iter/3)

suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(forcats))
suppressMessages(library(readr))
suppressMessages(library(lubridate))
suppressMessages(library(fredr))
suppressMessages(library(cli))
suppressMessages(library(glue))
suppressMessages(library(tidybayes))
suppressMessages(library(jsonlite))

election_day = as.Date("2020-11-03")
start_date = ymd("2020-01-15")
from_date = as.Date(opt$date)
n_days = ceiling(as.numeric(election_day - start_date) / 3) + 1
n_weeks = ceiling(as.numeric(election_day - start_date) / 21) + 1
wnum = (0:(n_days-1)) / 7
week_frac = if_else(floor(wnum)==wnum, 0, wnum - floor(wnum))
to_date = function(day)  election_day - 3*(n_days - day)

cur_dem = 47
cur_gop = 53
fixed_dem = 35
fixed_gop = 30

polls_model_path = file.path(opt$model_dir, "polls")
natl_model_path = file.path(opt$model_dir, "natl-prior")
race_model_path = file.path(opt$model_dir, "race-prior")

# past race results and other covariates
natl_d = suppressMessages(read_csv("data/generic_ballot.csv")) %>%
    mutate(lag_pop = lag(logit_pop)) %>%
    drop_na
race_d = suppressMessages(read_csv("data/senate_combined.csv")) %>%
    select(year, abbr, level, ref_level, pres_avg, senate, inc, pres, midterm) %>%
    filter(year==2020 | !is.na(level))
race_prior_d2020 = filter(race_d, year==2020) %>%
    select(-level) %>%
    arrange(abbr)

state_abbr = read_rds("data/state_data_2016.rdata") %>%
    select(state, abbr) %>%
    filter(!str_detect(state, "CD-")) %>%
    bind_rows(tibble(state="Georgia-S", abbr="GA-S")) %>%
    semi_join(race_prior_d2020, by="abbr") %>%
    arrange(abbr)



###
### Download new polling data
###
cli_h1("Downloading data")
source("model/get_data.R")

# presidential approval, general election polls, and distribution of past polling errors
pres_appr = get_approval() # mean pct. pt. gap appr-disappr
cli_alert_success("Presidential approval polling downloaded.")

old_polls = suppressMessages(read_csv("docs/polls.csv", col_types="cilcd"))
polls_d = get_elec_polls()
polls_d %>%
    select(date, race, national, firm, dem) %>%
    write_csv("docs/polls.csv")
if (from_date == Sys.Date() &&
        isTRUE(all.equal(old_polls, select(polls_d, date, race, national, firm, dem)))) {
    cli_alert_warning("No new polls.")
    system("osascript -e 'display notification \"No new polls.\" with title \"Senate Model\"'")
    system("osascript -e beep"); system("osascript -e beep")
    Sys.sleep(10)
}
cli_alert_success("{nrow(polls_d)} election poll{?s} downloaded.")

poll_errors = read_rds("output/polling_error.rdata")

# economy
fredr_set_key(Sys.getenv("FRED_KEY"))
econ_params = list(
    series = c("UNRATE", "A939RX0Q048SBEA", "CPIAUCSL", "AHETPI"),
    units = c("lin", "pc1", "pc1", "pc1")
)
econ = pmap_dfr(econ_params, ~ fredr(.x, observation_start=start_date-365,
                              frequency="q", units=.y)) %>%
    pivot_wider(names_from=series_id, values_from=value) %>%
    filter(date <= from_date) %>%
    transmute(unemp = as.numeric(UNRATE) / 100,
              gdp = as.numeric(A939RX0Q048SBEA) / 100,
              infl = as.numeric(CPIAUCSL) / 100,
              earn = as.numeric(AHETPI) / 100) %>%
    fill(everything()) %>%
    tail(1)
cli_alert_success("Economic data downloaded.")


natl_prior_d2020 = tibble(year=2020, pres=-1, midterm=0, unemp=econ$unemp,
                          gdp=econ$gdp, infl=econ$infl,
                          earn=econ$earn, appr=pres_appr,
                          lag_pop=tail(natl_d, 1)$logit_pop)


###

###
### Fitting prior model
###
cli_h1("Building prior")
suppressMessages(library(cmdstanr))
suppressMessages(library(rstanarm))
options(mc.cores=4)
source("model/get_models.R")

natl_model = get_natl_prior_m(natl_model_path, natl_d, opt$recompile)
natl_prior_pred = as.numeric(posterior_predict(natl_model, newdata=natl_prior_d2020))
cli_alert_success("National prior predictions made.")

race_model = get_race_prior_m(race_model_path, race_d, opt$recompile)
race_prior_pred = posterior_predict(race_model, newdata=race_prior_d2020)
race_prior_mean = colMeans(race_prior_pred)
race_x = apply(race_prior_pred, 2, function(x) x - mean(x))
race_prior_cov = (t(race_x) %*% race_x) / (nrow(race_x) - 1)
cli_alert_success("Race prior predictions made.")

###
### Fitting main model
###
cli_h1("Fitting main model")

model_d = compose_data(polls_d, .n_name = n_prefix("N"),
                       N_race = nrow(race_prior_d2020),
                       D = n_days,
                       W = n_weeks,
                       D_W = ceiling(n_days/n_weeks),
                       week_frac,
                       week_day = floor(wnum) + 1,
                       prior_natl_mean = mean(natl_prior_pred),
                       prior_natl_sd = sd(natl_prior_pred),
                       prior_race_mean = race_prior_mean,
                       prior_race_cov = race_prior_cov,
                       prior_rv_bias = 0.011,
                       prior_lv_bias = 0.0,
                       prior_a_bias = 0.02,
                       lv_rv_ratio = 8,
                       poll_errors,
)

polls_model = get_polls_m(polls_model_path, opt$recompile)
cli_alert_success("Model loaded.")

# TODO incorporate inv_metric stuff
fit_polls = polls_model$sample(data=model_d, num_chains=3, num_samples=opt$iter/3,
                               num_warmup=300, num_cores=4, adapt_delta=0.97,
                               stepsize=0.015)
cli_alert_success("Model successfully fit.")

raw_draws = posterior::as_draws_df(fit_polls$draws())

###
### Output predictions
###
cli_h1("Saving predictions")

natl_draws = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("natl_dem")) %>%
    pivot_longer(cols=starts_with("natl_dem"), names_to="day",
                 names_pattern="natl_dem\\[(.+)\\]", values_to="natl_dem") %>%
    mutate_if(is.character, as.numeric)

natl_final = filter(natl_draws, day==n_days) %>% pull(natl_dem)

race_draws = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("race_dem")) %>%
    pivot_longer(cols=starts_with("race_dem"), names_to=c("day", "race"),
                 names_pattern="race_dem\\[(.+),(.+)\\]", values_to="race_dem") %>%
    mutate(day = as.numeric(day),
           race_num = as.numeric(race),
           race = state_abbr$abbr[race_num])

race_summary = race_draws %>%
    filter(day == max(day)) %>%
    select(-.chain, -.iteration, -day, -race_num) %>%
    group_by(race) %>%
    summarize(prob = mean(race_dem > 0.5),
              dem_q05 = quantile(race_dem, 0.05),
              dem_q25 = quantile(race_dem, 0.25),
              dem_exp = median(race_dem),
              dem_q75 = quantile(race_dem, 0.75),
              dem_q95 = quantile(race_dem, 0.95)) %>%
    left_join(state_abbr, by=c("race"="abbr")) %>%
    left_join(select(race_prior_d2020, race=abbr, inc), by="race") %>%
    rename(race_name = state) %>%
    select(race, race_name, everything())

sims = race_draws %>%
    filter(day == max(day)) %>%
    group_by(.draw) %>%
    group_map(~ list(dem = as.integer(.$race_dem > 0.5)))

seats = race_draws %>%
    filter(day == max(day)) %>%
    left_join(select(race_prior_d2020, race=abbr, inc), by="race") %>%
    mutate(pickup_dem = inc != 1 & race_dem > 0.5,
           pickup_gop = inc == 1 & race_dem <= 0.5) %>%
    group_by(.draw) %>%
    summarize(dem_won = sum(race_dem > 0.5) + 2 - 1,
              dem_pickup = sum(pickup_dem),
              gop_pickup = sum(pickup_gop)) %>%
    mutate(dem_seats = fixed_dem + dem_won)

firm_ct = polls_d %>%
    group_by(firm) %>%
    summarize(n=n())
firm_ids = tibble(firm=polls_d$firm, id=model_d$firm) %>%
    distinct()
firms = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("house_effects")) %>%
    pivot_longer(cols=starts_with("house_effects"), names_to="id",
             names_pattern="house_effects\\[(.+)\\]", values_to="effect") %>%
    mutate_if(is.character, as.numeric) %>%
    group_by(id) %>%
    summarize(effect = median(effect)) %>%
    left_join(firm_ids, by="id") %>%
    left_join(firm_ct, by="firm") %>%
    select(firm, n, effect)


cli_alert_success("Outputs prepared.")

pr_presidency = read_json("../president/docs/estimate.json")$prob
entry = tibble(
    date = from_date,
    s_exp = median(seats$dem_seats),
    s_min = min(seats$dem_seats),
    s_q05 = quantile(seats$dem_seats, 0.05),
    s_q25 = quantile(seats$dem_seats, 0.25),
    s_q75 = quantile(seats$dem_seats, 0.75),
    s_q95 = quantile(seats$dem_seats, 0.95),
    s_max = max(seats$dem_seats),
    natl_exp = median(natl_final),
    natl_q05 = quantile(natl_final, 0.05),
    natl_q25 = quantile(natl_final, 0.25),
    natl_q75 = quantile(natl_final, 0.75),
    natl_q95 = quantile(natl_final, 0.95),
    pr_presidency = pr_presidency,
    pr_tie = mean(seats$dem_seats == 50), 
    prob = pr_tie*pr_presidency + mean(seats$dem_seats >= 51),
    dem_pickup_exp = median(seats$dem_pickup),
    dem_pickup_q05 = quantile(seats$dem_pickup, 0.05),
    dem_pickup_q25 = quantile(seats$dem_pickup, 0.25),
    dem_pickup_q75 = quantile(seats$dem_pickup, 0.75),
    dem_pickup_q95 = quantile(seats$dem_pickup, 0.95),
    gop_pickup_exp = median(seats$gop_pickup),
    gop_pickup_q05 = quantile(seats$gop_pickup, 0.05),
    gop_pickup_q25 = quantile(seats$gop_pickup, 0.25),
    gop_pickup_q75 = quantile(seats$gop_pickup, 0.75),
    gop_pickup_q95 = quantile(seats$gop_pickup, 0.95)
)
entry = mutate_if(entry, is.numeric, ~ round(., 4))

natl_intent_tbl = natl_draws %>%
    group_by(day) %>%
    summarize(natl_exp = median(natl_dem),
              natl_q05 = quantile(natl_dem, 0.05),
              natl_q25 = quantile(natl_dem, 0.25),
              natl_q75 = quantile(natl_dem, 0.75),
              natl_q95 = quantile(natl_dem, 0.95)) %>%
    mutate(date = to_date(day))

output = append(as.list(entry), list(
    time = Sys.time(),
    n_polls = nrow(polls_d),
    natl_intent = natl_intent_tbl,
    firm_effects = firms,
    races = race_summary,
    prior_natl_mean = mean(plogis(natl_prior_pred)),
    prior_natl_moe = 1.645*sd(plogis(natl_prior_pred)),
    hist = map_int(30:70, ~ sum(seats$dem_seats == .))
))

with(output, cat(glue("
    ===========================================
     2020 U.S. Senate Forecast
     {as.character(Sys.Date(), format='%B %d, %Y')}
    -------------------------------------------
     Forecast from: {as.character(from_date, format='%B %d, %Y')}
     {round((election_day - from_date)/7)} weeks until the election.
     {n_polls} polls.

     Median seat estimate:           {round(s_exp)}
     Estimated seat range:           {round(s_q05)} - {round(s_q95)}
     Median DEM pickup:              {dem_pickup_exp}
     Median GOP pickup:              {gop_pickup_exp}
     Probability of taking control:  {round(100*prob)}%
    ===========================================
    ")))
cat("\n\n")
system("osascript -e 'display notification \"Model run complete.\" with title \"Senate Model\"'")

if (opt$dry) quit("no")

if (file.exists(opt$history_file)) {
    history = bind_rows(
        suppressMessages(read_csv(opt$history_file)),
        entry
    ) %>%
        group_by(date) %>%
        slice(n()) %>%
        ungroup %>%
        arrange(date)
} else {
    history = entry
}
write_csv(history, opt$history_file, na="")

race_entry = race_summary %>%
    mutate(date = from_date) %>%
    mutate_if(is.numeric, ~ round(., 4)) %>%
    select(date, everything())

fname = str_c(dirname(opt$history_file), "/race_", basename(opt$history_file))
if (file.exists(fname)) {
    race_history = bind_rows(
        suppressMessages(read_csv(fname)),
        race_entry
    ) %>%
        group_by(race, date) %>%
        slice(n()) %>%
        ungroup %>%
        arrange(date)
} else {
    race_history = race_entry
}
write_csv(race_history, fname, na="")

# only save full output if current run
if (from_date != Sys.Date()) quit("no")

write_json(output, opt$output_file, auto_unbox=T, digits=7)
write_json(sims, opt$sims_file, auto_unbox=T, digis=4)

cli_alert_success("Model outputs saved.")
cat("\n")
