get_natl_prior_m = function(path, natl_d, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    if (file.exists(compiled_model)) {
        model = read_rds(compiled_model)
    } else {
        model = stan_glm(logit_pop ~ pres:(appr+earn+gdp+unemp) +
                                   midterm*pres + lag_pop,
                               data=natl_d, prior=normal(), chains=1, warmup=500)
        write_rds(model, compiled_model, compress="gz")
    }

    model
}

get_race_prior_m = function(path, race_d, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    if (file.exists(compiled_model)) {
        model = read_rds(compiled_model)
    } else {
        race_d = filter(race_d, year != 2020) %>%
            drop_na
        model = stan_lmer(level ~ ref_level + pres_avg + senate*inc +
                              pres*inc*midterm + (1|year),
                          data=race_d, prior=cauchy(), chains=1)

        write_rds(model, compiled_model, compress="gz")
    }

    model
}


get_polls_m = function(path, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    #if (file.exists(compiled_model)) {
    #    model = read_rds(compiled_model)
    #} else {
    #    model = stan_model("stan/polls.stan")
    #    write_rds(model, compiled_model, compress="gz")
    #}
    #model
    cmdstan_model("stan/polls.stan", cpp_options=list(CXXFLAGS="-Ofast"))
}

