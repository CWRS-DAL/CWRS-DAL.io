
#------------------ tests ------------------

library("testthat")
library("brms")
library("stringr")
source("R/functions.R")

#------------------ test models ------------------

seed <- 1

# data:

data <- withr::with_seed(seed, {
  tibble(
    y = rnorm(10),
    ycens = "none"
  )
})

data_ar <- withr::with_seed(seed, {
  tibble(
    date = seq(as.Date("2020-01-01"), as.Date("2020-01-31"), by = "1 day"),
    A = arima.sim(list(ar = .7), length(date)),
    B = arima.sim(list(ar = .7), length(date))
  ) %>%
    mutate(
      date = as.numeric(date),
      date = date - min(date) + 1
    ) %>%
    pivot_longer(
      c(A, B),
      names_to = "series",
      values_to = "y"
    ) %>%
    arrange(series, date) %>%
    group_by(series) %>%
    mutate(
      d_x = date - lag(date),
      d_x = replace_na(d_x, 0),
      d_x = as.numeric(d_x)
    ) %>%
    ungroup()
})

phi_car1 <- .45
p_ret <- .6 # proportion retained

data_car1 <- withr::with_seed(seed, {
  tibble(
    x = 1:200,
    y = arima.sim(list(ar = phi_car1), length(x))
  ) %>%
    slice_sample(prop = p_ret) %>%
    arrange(x) %>%
    mutate(
      x_lag = lag(x),
      d_x = replace_na(x - x_lag, 0) # spacing of observations
    )
})

# fitted models:

fit <- fit_stan_model(
  "R/tests/test",
  seed,
  bf(y | cens(ycens) ~ 1),
  data,
  prior(normal(0, 1), class = Intercept),
  car1 = FALSE,
  save_warmup = FALSE,
  chains = 3
)

form_ar <- bf(y ~ ar(time = date, gr = series), sigma ~ series)
prior_ar <- prior(normal(0, 1), class = Intercept)

fit_ar <- fit_stan_model(
  "R/tests/test_ar",
  seed,
  form_ar,
  data_ar,
  prior_ar,
  save_warmup = FALSE,
  chains = 2
)

fit_ar2 <- brm(
  form_ar,
  data = data_ar,
  family = student(),
  seed = seed,
  chains = 2,
  file = "R/tests/test_ar_brm"
)

form_car1 <- bf(y ~ ar(time = x))

fit_car1 <- fit_stan_model(
  "R/tests/test_car1",
  seed,
  form_car1,
  data_car1,
  prior_ar,
  save_warmup = FALSE,
  chains = 2
)

# model and residuals from README.Rmd

rseed <- 2356

simdat <- readr::read_csv("data/data-simulated.csv")

bform <- bf(
  scaled_lead | cens(cens_lead) ~
  s(scaled_date_numeric, k = 10) +
    s(scaled_date_yday, bs = "cc", k = 10) +
    ar(time = scaled_date_numeric)
)

model <- fit_stan_model(
  "models/demo", rseed, bform, simdat,
  save_warmup = FALSE
)

#------------------ functions for figures and models ------------------

# calc_sigma():

test_that("calc_sigma yields the expected intercept for sigma", {
  extract_sigma <- function(x) {
    summary(x)$fixed %>%
      as_tibble(rownames = "param") %>%
      filter(str_detect(param, "^sigma_"))
  }

  extract_calc_sigma <- function(x) {
    x %>%
      calc_sigma() %>%
      filter(name == "b_sigma_Intercept") %>%
      select(-c(name, group)) %>%
      ggdist::point_interval(.point = median)
  }

  s1 <- extract_sigma(fit_ar)
  s2 <- extract_sigma(fit_ar2)

  s1a <- extract_calc_sigma(fit_ar)
  s2a <- extract_calc_sigma(fit_ar2)

  expect_equal(s1$Estimate[1], log(s1a$sigma), tolerance = 1e-2)
  expect_equal(s2$Estimate[1], log(s2a$sigma), tolerance = 1e-2)
})

test_that("calc_sigma() works for distributional and non-distributional models.", {
  non_dist <- calc_sigma(fit)
  dist <- calc_sigma(fit_ar)
  expect_equal(dim(dist), c(4e3, 6))
  expect_equal(dim(non_dist), c(3e3, 4))
})

# retrans():

test_that("retrans() works as expected.", {
  x <- rnorm(10, 10)
  x_trans <- scale(x)[, 1]
  x_ltrans <- scale(log(x))[, 1]
  x_retrans <- retrans(x_trans, x, log = FALSE)
  x_reltrans <- retrans(x_ltrans, x)
  # in a dataframe:
  data <- tibble(
    orig = x,
    trans = x_trans,
    ltrans = x_ltrans
  ) %>%
    mutate(
      re_trans = retrans(trans, orig, log = FALSE),
      re_ltrans = retrans(ltrans, orig),
    )
  expect_equal(x, x_retrans)
  expect_equal(x, x_reltrans)
  expect_equal(data$orig, data$re_trans)
  expect_equal(data$orig, data$re_ltrans)
})

# extract_resp():

test_that("extract_resp() returns error when input is wrong class.", {
  x <- tibble(formula = 1)
  expect_error(extract_resp(x))
})

test_that("extract_resp() extracts the correct model names.", {
  x <- extract_resp(model) %>%
    unlist()
  fit_resp <- extract_resp(fit) %>%
    unlist()
  fit_ar_resp <- extract_resp(fit_ar) %>%
    unlist()
  expect_equal(
    x,
    c(
      resp = "scaled_lead", cens = "cens_lead", gr_sigma = NA,
      gr_ar = NA, time_ar = "scaled_date_numeric"
    )
  )
  expect_equal(fit_resp, c(resp = "y", cens = "ycens", gr_sigma = NA))
  expect_equal(
    fit_ar_resp,
    c(resp = "y", cens = NA, gr_sigma = "series", gr_ar = "series", time_ar = "date")
  )
})

# ppc_km_nada():

test_that("ppc_km_nada() throws an error when the AR(1) term is missing", {
  expect_error(ppc_km_nada(data, fit, seed = seed))
})

test_that("ppc_km_nada() throws an error when the censoring term is missing", {
  expect_error(ppc_km_nada(data_ar, fit_ar, seed = seed))
})

test_that("ppc_km_nada() yields the same ecdf as NADA::cenfit()", {
  x <- ppc_km_nada(data, fit, draw_ids = 34, seed = seed, car1 = FALSE) %>%
    filter(type == "Observations") %>%
    select(-c(type, .draw))
  y <- cenfit(obs = data$y, censored = data$ycens == "left") %>%
    summary() %>%
    as_tibble()
  expect_equal(x, y)
})

# modify_stancode():

test_that("modify_stancode() adds expected number of characters to stan code.", {
  n1 <- str_count(fit$model)
  n1_mod <- modify_stancode(fit$model) %>%
    str_count()
  d1 <- n1_mod - n1

  n2 <- str_count(fit_ar$model)
  n2_mod <- modify_stancode(fit_ar$model) %>%
    str_count()
  d2 <- n2_mod - n2
  expect_equal(d1, 35) # 35 chars to data block
  expect_equal(d2, 73) # 73 chars overall
})

# add_car1():

extract_acf <- function(x) {
  x %>%
    acf(plot = FALSE) %>%
    with(acf) %>%
    as.numeric()
}

add_car1_input <- data_ar %>%
  mutate(.index = 1, `ar[1]` = .7, .epred = 0)

test_that(
  "add_car1() models an AR(1).",
  {
    car1 <- add_car1(add_car1_input, "y", gr_vars = c(".index", "series")) %>%
      mutate(r = y - .epred)
    autocorr <- extract_acf(car1$r)
    autocorr2 <- arima(add_car1_input$y, order = c(1, 0, 0)) %>%
      residuals() %>%
      extract_acf()
    expect_equal(autocorr, autocorr2, tolerance = .025)
  }
)

test_that(
  "add_car1() models an irregularly-spaced AR(1).",
  {
    add_car1_sub <- withr::with_seed(3526, {
      slice_sample(add_car1_input, prop = .7) %>%
        arrange(series, date) %>%
        group_by(series) %>%
        mutate(
          d_x = date - lag(date),
          d_x = replace_na(d_x, 0),
          d_x = as.numeric(d_x)
        ) %>%
        ungroup()
    })

    car1 <- add_car1(add_car1_sub, "y", gr_vars = c(".index", "series")) %>%
      mutate(r = y - .epred)

    car1_cor <- car1 %>%
      group_by(series) %>%
      mutate(
        r_lag = lag(r),
        y_lag = lag(y)
      ) %>%
      ungroup() %>%
      filter(d_x == 1) %>%
      summarize(
        cor_r = cor(r_lag, r),
        cor_y = cor(y_lag, y)
      )
    expect_equal(car1_cor$cor_r, 0, tolerance = .1)
  }
)

test_that(
  "add_car1() yields the same results as
  tidybayes::add_epred_draws() for an AR(1) model.",
  {
    tbl1 <- tidybayes::add_epred_draws(data_ar, fit_ar) %>%
      rename(.index = .draw) %>%
      ungroup()
    ar1 <- as_draws_df(fit_ar, "ar[1]") %>%
      as_tibble() %>%
      select(-c(.chain, .iteration))
    tbl2 <- tidybayes::add_epred_draws(data_ar, fit_ar, incl_autocor = FALSE) %>%
      left_join(ar1, by = ".draw") %>%
      rename(.index = .draw) %>%
      add_car1("y", gr_vars = c(".index", "series")) %>%
      select(-`ar[1]`)
    expect_equal(tbl1$.epred, tbl2$.epred)
  }
)

# filter_car1():

test_that("filter_car1() works the same as stats::filter() for regularly spaced data.", {
  x <- rnorm(100)
  d_x <- rep(1, length(x))
  f1 <- filter_car1(x, .5, s = d_x)
  f2 <- stats::filter(x, filter = .5, method = "recursive")
  expect_equal(as.numeric(f1), as.numeric(f2))
})

test_that("filter_car1() returns an error for NAs.", {
  x <- rnorm(100)
  x[50] <- NA
  d_x <- rep(1, length(x))
  expect_error(filter_car1(x, .5, s = d_x))
})

test_that("filter_car1() generates an irregularly-sampled AR(1) process.", {
  phi <- .8

  data <- withr::with_seed(356, {
    tibble(
      x = rnorm(200)
    ) %>%
      tibble::rowid_to_column()
  })

  sub <- withr::with_seed(341, {
    slice_sample(data, prop = .7)
  }) %>%
    arrange(rowid) %>%
    mutate(
      d_x = replace_na(rowid - lag(rowid), 0),
      car1 = filter_car1(x, phi, s = d_x)
    )

  fit <- data %>%
    left_join(sub) %>%
    with(arima(car1, order = c(1, 0, 0)))

  expect_equal(phi, as.numeric(fit$coef[1]), tolerance = .1)
})

# add_car1_err():

phi <- .7

car1_input <- crossing(.index = 1:2, location = letters[1:2], rep = 1:200) %>%
  mutate(
    `ar[1]` = phi,
    nu = 1e3,
    sigma = 1,
    .epred = 0
  )

test_that("add_car1_err() generates a regular AR(1) process.", {
  data_car1 <- withr::with_seed(32567, {
    car1_input %>%
      mutate(d_x = 1) %>%
      add_car1_err(gr_vars = c(".index", "location")) %>%
      filter(.index == 1, location == "a")
  })

  fit <- arima(data_car1$.prediction, order = c(1, 0, 0))

  expect_equal(as.numeric(fit$coef[1]), phi, tolerance = .1)
})

test_that("add_car1_err() generates an irregular AR(1) process.", {
  sub <- withr::with_seed(219, {
    car1_input %>%
      group_by(.index, location) %>%
      slice_sample(prop = .6) %>%
      ungroup() %>%
      arrange(.index, location, rep) %>%
      mutate(d_x = replace_na(rep - lag(rep), 0))
  })

  data_car1 <- withr::with_seed(32567, {
    sub %>%
      add_car1_err(gr_vars = c(".index", "location")) %>%
      filter(.index == 1, location == "a")
  })

  car1_input %>%
    left_join(data_car1)

  fit <- arima(data_car1$.prediction, order = c(1, 0, 0))

  expect_equal(as.numeric(fit$coef[1]), phi, tolerance = .1)
})

# extract_params():

test_that("extract_params() yields expected output.", {
  fit_ar_pars <- extract_params(fit_ar, draw_ids = 2000)
  reference <- tibble::tribble(
    ~`ar[1]`, ~.chain, ~.iteration, ~.draw, ~.index,
    0.748742,      2L,       1000L,  2000L,       1
  )
  expect_equal(fit_ar_pars, reference)
})


# add_pred_draws_car1():

test_that("add_pred_draws_car1() returns an error for incorrect 'type'", {
  expect_error(add_pred_draws_car1(data_ar, fit_ar, type = "wrong type"))
})

test_that(
  "add_pred_draws_car1() yields the same predictions as fitted.brmsfit()",
  {
    preds <- add_pred_draws_car1(
      input = simdat, object = model,
      car1 = FALSE, draw_ids = 4000
    )

    fitted_vals <- fitted(
      model, simdat,
      incl_autocor = FALSE,
      robust = TRUE, draw_ids = 4000
    ) %>%
      as_tibble()

    expect_equal(preds$.epred, fitted_vals$Estimate)
  }
)

test_that(
  "add_pred_draws_car1() adds the CAR(1) filter correctly.",
  {
    # function:
    preds <- add_pred_draws_car1(simdat, model, draw_ids = 4000)

    # by "hand":

    phi <- as_draws_df(model, variable = "ar[1]") %>%
      as_tibble() %>%
      filter(.draw == 4000) %>%
      pull(`ar[1]`)

    preds_noar <- add_pred_draws_car1(
      simdat, model,
      car1 = FALSE, draw_ids = 4000
    ) %>%
      # add the CAR(1) structure:
      group_by(.index) %>%
      mutate(
        r_lag = replace_na(lag(scaled_lead - .epred), 0),
        .epred = .epred + r_lag * phi^d_x
      ) %>%
      ungroup()

    expect_equal(preds$.epred, preds_noar$.epred)
  }
)

test_that(
  "add_pred_draws_car1() yields the same results as tidybayes::add_epred_draws()
  for a regular AR(1) fit.",
  {
    compare_preds <- function(x, ...) {
      x %>%
        select(date, series, y, d_x, .row, .epred, ...) %>%
        arrange(.draw) %>%
        ungroup()
    }

    preds1 <- tidybayes::add_epred_draws(data_ar, fit_ar) %>%
      compare_preds(.draw)

    preds2 <- add_pred_draws_car1(data_ar, fit_ar, draw_ids = 1:2000) %>%
      compare_preds(.draw = .index)

    expect_equal(preds1, preds2)
  }
)

test_that("add_pred_draws_car1() joins params correctly for README.Rmd model.", {
  these_ids <- c(452, 3298)
  preds <- add_pred_draws_car1(simdat, model, draw_ids = these_ids, type = "prediction")
  draws <- as_draws_df(model, c("sigma", "nu", "ar[1]")) %>%
    as_tibble() %>%
    filter(.draw %in% these_ids)
  sig1 <- draws %>%
    select(sigma)
  sig2 <- preds %>%
    ungroup() %>%
    distinct(sigma)
  nu1 <- draws %>%
    select(nu)
  nu2 <- preds %>%
    ungroup() %>%
    distinct(nu)
  ar1 <- draws %>%
    select(`ar[1]`)
  ar2 <- preds %>%
    ungroup() %>%
    distinct(`ar[1]`)
  expect_equal(sig1, sig2)
  expect_equal(nu1, nu2)
  expect_equal(ar1, ar2)
})

test_that("add_pred_draws_car1() joins params correctly for test model.", {
  these_ids <- c(452, 1298)
  preds <- add_pred_draws_car1(data, fit, draw_ids = these_ids, car1 = FALSE, type = "prediction")
  draws <- as_draws_df(fit, c("sigma", "nu")) %>%
    as_tibble() %>%
    filter(.draw %in% these_ids)
  sig1 <- draws %>%
    select(sigma)
  sig2 <- preds %>%
    ungroup() %>%
    distinct(sigma)
  nu1 <- draws %>%
    select(nu)
  nu2 <- preds %>%
    ungroup() %>%
    distinct(nu)
  expect_equal(sig1, sig2)
  expect_equal(nu1, nu2)
})

test_that("add_pred_draws_car1() joins params correctly for test distributional model.", {
  these_ids <- c(452, 1298)
  preds <- add_pred_draws_car1(data_ar, fit_ar, draw_ids = these_ids)
  draws <- as_draws_df(fit_ar, c("nu", "ar[1]")) %>%
    as_tibble() %>%
    filter(.draw %in% these_ids)
  sig1 <- calc_sigma(fit_ar) %>%
    filter(.draw %in% these_ids) %>%
    select(series = group, sigma)
  sig2 <- preds %>%
    ungroup() %>%
    distinct(series, sigma)
  nu1 <- draws %>%
    select(nu)
  nu2 <- preds %>%
    ungroup() %>%
    distinct(nu)
  ar1 <- draws %>%
    select(`ar[1]`)
  ar2 <- preds %>%
    ungroup() %>%
    distinct(`ar[1]`)
  expect_equal(sig1, sig2)
  expect_equal(nu1, nu2)
})

test_that(
  "add_pred_draws_car1() fits a CAR(1) model that accounts for the autocorrelation
  structure in an irregularly sampled AR(1).",
  {
    preds <- add_pred_draws_car1(data_car1, fit_car1, draw_ids = 1:2000) %>%
      ggdist::median_qi(.epred) %>%
      mutate(r = y - .epred)
    full <- tibble(
      x = seq_len(max(data_car1$x))
    ) %>%
      left_join(preds, by = "x")
    arima_car1_r <- arima(full$r, order = c(1, 0, 0))
    arima_car1 <- arima(full$y, order = c(1, 0, 0))
    expect_equal(coef(arima_car1_r)[1], c(ar1 = 0), tolerance = .1)
    expect_equal(coef(arima_car1)[1], c(ar1 = phi_car1), tolerance = .1)
  }
)

test_that(
  "add_pred_draws_car1() generates CAR(1) predictions.",
  {
    preds_car1 <- withr::with_seed(
      135,
      {
        add_pred_draws_car1(data_car1, fit_car1, draw_ids = 420, type = "prediction")
      }
    )
    full <- tibble(
      x = seq_len(max(data_car1$x))
    ) %>%
      left_join(preds_car1, by = "x")
    arima_car1 <- arima(full$.prediction, order = c(1, 0, 0))
    expect_equal(coef(arima_car1)[1], c(ar1 = unique(preds_car1$`ar[1]`)), tolerance = .1)
  }
)

# calc_ll():

ll_brm <- brms::log_lik(fit_ar)
ll_myfn <- add_pred_draws_car1(data_ar, fit_ar, draw_ids = 1:2000) %>%
  ungroup() %>%
  calc_ll("y", cens = FALSE) %>%
  pivot_wider(c(.draw, .chain, .iteration), names_from = .row, values_from = log_lik) %>%
  select(matches("^\\d")) %>%
  as.matrix()

colnames(ll_myfn) <- NULL

test_that("calc_ll() returns the same likelihoods as brms::log_lik()", {
  expect_equal(ll_brm, ll_myfn)
})

ll_brm <- brms::log_lik(model, incl_autocor = FALSE)
ll_myfn <- add_pred_draws_car1(simdat, model, draw_ids = 1:8000, car1 = FALSE) %>%
  ungroup() %>%
  calc_ll("scaled_lead", "cens_lead") %>%
  pivot_wider(.index, names_from = .row, values_from = log_lik) %>%
  select(matches("^\\d")) %>%
  as.matrix()

colnames(ll_myfn) <- NULL

test_that("calc_ll() returns the same likelihoods as brms::log_lik()", {
  expect_equal(ll_brm, ll_myfn)
})

# loo_cv():

loo1 <- loo(model, incl_autocor = FALSE)
loo2 <- loo_cv(simdat, model, car1 = FALSE, draw_ids = 1:8000)

loo1a <- loo(fit_ar)
loo2a <- loo_cv(data_ar, fit_ar, censoring = FALSE, draw_ids = 1:2000)

test_that("brms::loo() yields the same results as loo_cv()", {
  expect_equal(loo1$estimate, loo2$estimate)
  expect_equal(loo1a$estimate, loo2a$estimate)
})

# fit_stan_model():

test_that("function loads the correct model", {
  expect_equal(class(fit), "brmsfit")
  expect_equal(as.character(fit$formula)[1], "y | cens(ycens) ~ 1")
})

test_that("fit_stan_model() yields the same results as brms for AR(1) model.", {
  ar_mod1 <- summary(fit_ar)$cor_pars[, 1]
  ar_mod2 <- summary(fit_ar2)$cor_pars[, 1]
  expect_equal(ar_mod1, ar_mod2, tolerance = 1e-2)
})

test_that(
  "fit_stan_model() fits a CAR(1) model that recovers the parameters used to generate the data.",
  {
    expect_equal(summary(fit_car1)$cor_pars$Estimate, phi_car1, tolerance = .1)
  }
)

# impute_censored():

test_that("impute_censored() replaces censored values with predictions.", {
  x <- tibble(
    x = 1:10,
    y = rnorm(10),
    .prediction = rnorm(10),
    cens = "left"
  )
  input <- select(x, x, y, cens)
  imp <- impute_censored(x, input, "y", "cens")
  expect_equal(imp$y, x$.prediction)
})

# simulate_residuals():

sim <- simulate_residuals(
  data, fit, 2000, seed,
  "R/tests/test_resid",
  car1 = FALSE, model = "lm",
  save_warmup = FALSE
)

resids <- tidybayes::add_residual_draws(
  newdata = data, object = fit, method = "posterior_epred"
)

test_that("simulate_residuals() returns the expected output.", {
  expect_equal(sim$residuals, resids)
})

test_that("simulate_residuals() informs user when tidybayes::add_residual_draws() is used.", {
  expect_message(
    suppressWarnings(
      simulate_residuals(
        data, fit, 2000, seed,
        "R/tests/test_resid",
        car1 = FALSE, model = "lm",
        save_warmup = FALSE
      )
    ),
    regexp = "Generating residual draws using tidybayes::add_residual_draws()"
  )
})

test_that("posterior_smooths() works for CAR(1) model.", {
  sm1 <- posterior_smooths(model, smooth = "s(scaled_date_numeric,k=10)")
  sm2 <- posterior_smooths(model, smooth = "s(scaled_date_yday,bs=\"cc\",k=10)")
  intercept <- as_draws_df(model, "Intercept") %>%
    as_tibble()
  preds1 <- sm1 + sm2 + intercept$Intercept
  preds2 <- posterior_epred(model, incl_autocor = FALSE)
  preds_tidy1 <- preds1 %>%
    as_tibble(.name_repair = "unique") %>%
    pivot_longer(cols = everything())
  preds_tidy2 <- add_pred_draws_car1(simdat, model, car1 = FALSE, draw_ids = 1:8000) %>%
    ungroup() %>%
    arrange(.draw, .row)
  expect_equal(preds1, preds2)
  expect_equal(preds_tidy1$value, preds_tidy2$.epred)
})
