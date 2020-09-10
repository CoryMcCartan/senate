appr_url = "https://projects.fivethirtyeight.com/trump-approval-data/approval_polllist.csv"
elec_url = "https://projects.fivethirtyeight.com/polls-page/senate_polls.csv"
generic_url = "https://projects.fivethirtyeight.com/polls-page/generic_ballot_polls.csv"

get_approval = function() {
    pres_appr = suppressMessages(read_csv(appr_url)) %>%
        transmute(approval = approve/(approve + disapprove),
                  date = mdy(enddate)) %>%
        filter(date <= from_date) %>%
        arrange(date) %>%
        tail(20) %>%
        summarize(appr = mean(qlogis(approval))) %>%
        pull
}

get_elec_polls = function(write=F, min_date=ymd("2020-01-01")) {
    polls_raw = suppressWarnings(suppressMessages(read_csv(elec_url)))
    polls_d = polls_raw %>%
        filter(cycle == 2020, mdy(start_date) >= min_date) %>%
        mutate(state = if_else(race_id == 7780, "Georgia-S", state),
               candidate_party = str_to_lower(candidate_party)) %>%
        select(question_id, date1=start_date, date2=end_date, state,
               party=candidate_party, pct, sample_size, population, pollster) %>%
        mutate(party = if_else(state=="Maine" & party %in% c("ind", "gre"), "dem", party)) %>% # maine IRV
        filter(party %in% c("dem", "rep")) %>%
        group_by(question_id, date1, date2, state, party, sample_size,
                 population, pollster) %>%
        summarize(pct = sum(pct)) %>%
        #filter(n() == 1) %>%
        ungroup %>%
        pivot_wider(names_from=party, values_from=pct) %>%
        mutate(national = F,
               date = mdy(date1) + (mdy(date2) - mdy(date1))/2,
               day = n_days - ceiling(as.numeric(election_day - date) / 3),
               week = n_weeks - ceiling(as.numeric(election_day - date) / 21),
               total = dem + rep,
               dem = dem / total,
               sample = round(sample_size * total/100),
               logit_inflate = 1/dem + 1/(1 - dem),
               var_poll = logit_inflate^2 * dem * (1 - dem) / sample,
               firm = as.character(fct_lump(pollster, 50)),
               type_rv = population=="rv" | population == "v",
               type_lv = population=="lv",
               type_a = population=="a") %>%
        filter(date >= start_date, date <= from_date, !is.na(dem)) %>%
        select(race=state, national, date, day, week, firm, sample,
               type_rv, type_lv, type_a, dem, var_poll)

    if (write) write_csv(polls_d, "data/polls/race_polls.csv")

    generic_raw = suppressMessages(read_csv(generic_url))
    generic_d = generic_raw %>%
        ungroup %>%
        mutate(national = T, state=NA,
               date = mdy(start_date) + (mdy(end_date) - mdy(start_date))/2,
               day = n_days - ceiling(as.numeric(election_day - date) / 3),
               week = n_weeks - ceiling(as.numeric(election_day - date) / 21),
               twoparty = dem + rep,
               dem = dem / twoparty,
               sample = round(sample_size * twoparty/100),
               logit_inflate = 1/dem + 1/(1 - dem),
               var_poll = logit_inflate^2 * dem*(1 - dem) / sample,
               firm = as.character(fct_lump(pollster, 30)),
               type_rv = population=="rv" | population == "v",
               type_lv = population=="lv",
               type_a = population=="a") %>%
        select(race=state, national, date, day, week, firm, sample,
               type_rv, type_lv, type_a, dem, var_poll) %>%
        filter(date >= start_date, date <= from_date, !is.na(dem))

    if (write) write_csv(generic_d, "data/polls/generic_polls.csv")

    bind_rows(polls_d, generic_d) %>%
        arrange(date) %>%
        left_join(state_abbr, by=c("race"="state")) %>%
        group_by(national) %>%
        mutate(abbr = if_else(national, "MA", abbr),
               date = as.character(date)) %>%
        select(abbr, day, date, week, everything(), -race) %>%
        mutate(race = match(abbr, state_abbr$abbr)) %>%
        drop_na

}
