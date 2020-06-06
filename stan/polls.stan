/************************************************
 * PRESIDENTIAL ELECTION POLLING MODEL          *
 * CORY McCARTAN                                *
 * (c) 2020                                     *
 ************************************************/

data {
    int N; // number of polls
    int D; // "days"
    int W; // "weeks"
    int D_W; // "days" per "week"
    int N_race;
    int N_firm;

    int<lower=1> day[N];
    int<lower=1> week[N];
    int<lower=-1> race[N];
    int<lower=1> firm[N];
    int<lower=1> week_day[D];
    real<lower=0, upper=1> week_frac[D];
    vector<lower=0, upper=1>[N] type_rv;
    vector<lower=0, upper=1>[N] type_lv;
    vector<lower=0, upper=1>[N] type_a;
    vector<lower=0, upper=1>[N] national;

    vector[N] dem;
    vector[N] var_poll;

    real prior_natl_mean;
    vector[N_race] prior_race_mean;
    cov_matrix[N_race] prior_race_cov; // prior model result correlation

    real<lower=0> prior_natl_sd;
    real<lower=0> prior_natl_poll_error;
    real<lower=0> prior_race_poll_error;
    real<lower=0> prior_generic_poll_error;
    real<lower=0> prior_generic_poll_bias;

    real prior_rv_bias;
    real prior_lv_bias;
    real prior_a_bias;
    real<lower=0> lv_rv_ratio;
}

transformed data {
    vector[N] logit_dem = logit(dem);
    cholesky_factor_cov[N_race] chol_race_cov
        = cholesky_decompose(prior_race_cov);
}

parameters {
    real<lower=0> sigma_natl; // day-to-day variance
    real<lower=0> sigma_race; // week-to-week variance

    real<lower=0> nonsamp_var;
    real generic_error;
    real natl_error;
    vector[N_race] race_error;
    real<lower=0> sd_firm; // hyperparameter for pollster errors
    real bias_rv; // registered voter poll bias
    real bias_lv; // likely voter poll bias
    real bias_a; // adult voter poll bias

    // reparametrization
    vector[D] delta_natl; // steps of random walk
    vector[N_race] delta_race[W]; // steps of randowm walk
    vector[N_firm] delta_firm; // polling firm standardized errors
}

transformed parameters {
    vector[D] mu_natl; // national voter intention, logit
    vector[N_race] mu_race[W]; // national voter intention, logit

    mu_natl[D] = prior_natl_mean + prior_natl_sd*delta_natl[D];
    mu_race[W] = prior_race_mean + chol_race_cov*delta_race[W];
    for (i in 1:(D-1)) {
        int d = D - i;
        mu_natl[d] = mu_natl[d+1] + sigma_natl*delta_natl[d];
    }
    for (i in 1:(W-1)) {
        int w = W - i;
        mu_race[w] = mu_race[w+1] + sqrt(D_W)*sigma_race*delta_race[w];
    }

}

model {
    // support for dem. in specific poll
    vector[N] val_poll = mu_natl[day] + 4*sd_firm*delta_firm[firm]
        + 4*bias_rv*type_rv + 4*bias_lv*type_lv + 4*bias_a*type_a +
        prior_generic_poll_error*generic_error + prior_generic_poll_bias;
    for (i in 1:N) {
        if (!national[i]) {
            val_poll[i] += prior_natl_poll_error*natl_error
                + prior_race_poll_error*race_error[race[i]];
            if (day[i] <= D_W)
                val_poll[i] += mu_race[week[i], race[i]];
            else
                val_poll[i] += week_frac[day[i]]*mu_race[week[i], race[i]]
                    + (1 - week_frac[day[i]])*mu_race[week[i] - 1, race[i]];
        }
    }

    logit_dem ~ normal(val_poll, sqrt(var_poll + nonsamp_var)); // modified binom. approx.

    // reparametrization
    delta_natl ~ std_normal();
    for (w in 1:W) {
        delta_race[w] ~ std_normal();
    }
    delta_firm ~ std_normal();

    // priors
    sigma_natl ~ gamma(3, 3/0.06);
    sigma_race ~ gamma(2, 2/0.02);
    generic_error ~ std_normal();
    natl_error ~ std_normal();
    race_error ~ std_normal();

    sd_firm ~ gamma(2, 2/0.05);
    bias_rv ~ normal(prior_rv_bias, 0.04);
    bias_lv ~ normal(prior_lv_bias, 0.04/lv_rv_ratio);
    bias_a ~ normal(prior_a_bias, 0.04);
}

generated quantities {
    real nonsamp_sd = sqrt(nonsamp_var);
    vector[D] natl_dem = inv_logit(mu_natl);
    vector[N_firm] house_effects = sd_firm * delta_firm;
    vector[N_race] race_dem[D];

    for (d in 1:D_W) {
        race_dem[d] = inv_logit(mu_natl[d] + mu_race[week_day[d]]);
    }
    for (d in (D_W+1):D) {
        race_dem[d] = inv_logit(mu_natl[d] + week_frac[d]*mu_race[week_day[d]]
            + (1 - week_frac[d])*mu_race[week_day[d] - 1]);
    }
}
