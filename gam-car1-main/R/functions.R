
library("dplyr")
library("purrr")
library("tidyr")
library("stringr")
library("data.table")
library("NADA")

#------------------ functions for figures and models ------------------

calc_sigma <- function(x) { # this function is currently used for testing only
  varnames <- extract_resp(x)
  distributional <- as.character(x$formula)[2] %>%
    str_detect("sigma ~ ")
  if (distributional) {
    groups <- unique(x$data[, varnames$gr_sigma])
  }
  draws <- x %>%
    as_draws_df("sigma", regex = TRUE) %>%
    as_tibble()
  if (distributional) {
    draws <- draws %>%
      pivot_longer(c(starts_with("b_sigma_"), -matches("^b_sigma_Intercept$"))) %>%
      mutate(value = b_sigma_Intercept + value) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      pivot_longer(starts_with("b_sigma_"), values_to = "sigma") %>%
      select(c(name, sigma, starts_with("."))) %>%
      mutate(
        sigma = exp(sigma),
        group = str_remove(name, varnames$gr_sigma) %>%
          str_extract("(?<=b_sigma_).+$") %>%
          str_replace("Intercept", groups[1])
      )
  }
  draws
}

# retransform predictions:

retrans <- function(x, scale_var, log = TRUE, recensor = FALSE, lcl = 2e-4) {
  rt <- if (log) {
    exp(
      x * sd(log(scale_var), na.rm = TRUE) +
        mean(log(scale_var), na.rm = TRUE)
    )
  } else {
    x * sd(scale_var, na.rm = TRUE) +
      mean(scale_var, na.rm = TRUE)
  }

  if (recensor) {
    rt <- if_else(rt < lcl, lcl, rt)
  }

  rt
}

# model a CAR(1) process:

add_car1 <- function(x, response, gr_vars = c(".index", "series")) {
  x <- x %>%
    ungroup() %>%
    as.data.table()

  x[,
    r_lag := replace_na(lag(get(response) - .epred), 0),
    by = c(gr_vars)
  ][,
    .epred := .epred + r_lag * `ar[1]`^d_x,
    by = c(gr_vars)
  ][
    ,
    r_lag := NULL
  ]

  as_tibble(x)
}

# CAR(1) filter:

filter_car1 <- function(x, phi, s) {
  if (sum(is.na(x)) > 0) stop("input series contains NAs")
  n <- length(x)
  xf <- rep(NA, n)
  xf[1] <- x[1]
  for (i in 2:n) {
    xf[i] <- phi^s[i] * xf[i - 1] + x[i]
  }
  xf
}

# add simulated CAR(1) error to model:

add_car1_err <- function(x, car1 = TRUE, gr_vars = c(".index", "series")) {
  x <- x %>%
    tibble::rowid_to_column() %>%
    as.data.table()

  x[,
    .err := rstudent_t(n = 1, df = nu, mu = 0, sigma = sigma),
    by = rowid
  ][
    ,
    rowid := NULL
  ]

  if (car1) {
    x[,
      .err := filter_car1(.err, phi = unique(`ar[1]`), s = d_x),
      by = c(gr_vars)
    ]
  }

  x %>%
    as_tibble() %>%
    mutate(.prediction = .epred + .err)
}

# extract parameters from model:

extract_params <- function(object, car1 = TRUE, draw_ids) {
  this_var <- if (car1) {
    "ar[1]"
  } else {
    "Intercept"
  }

  out <- as_draws_df(object, variable = this_var) %>%
    as_tibble() %>%
    filter(.draw %in% draw_ids) %>%
    mutate(.index = as.numeric(factor(.draw)))

  if (!car1) {
    # this adds .chain and .iteration to draws object, which would otherwise be missing
    # (but is needed for loo::loo())
    select(out, -Intercept)
  } else {
    out
  }
}

# calculate epred draws for car1 model:

add_pred_draws_car1 <- function(input,
                                object,
                                type = "epred",
                                car1 = TRUE,
                                draw_ids = 1:4000,
                                ...) {
  if (!type %in% c("epred", "prediction")) stop("'type' must be either 'prediction' or 'epred'")

  data_vars <- glue::glue("^{names(input)}$") %>% # group output by columns in "data"
    paste(collapse = "|")

  # extract variables from model object:

  varnames <- extract_resp(object) # responses, etc ...

  params <- extract_params(object, car1, draw_ids)

  gr_vars <- c(".index", varnames$gr_ar) %>% # for CAR1 error
    na.omit()

  order_vars <- c(".index", varnames$gr_ar, varnames$time_ar) %>%
    na.omit() # for CAR1 error

  # generate predictions without AR term:
  preds <- tidybayes::add_epred_draws(
    input, object,
    incl_autocor = FALSE,
    draw_ids = draw_ids, dpar = TRUE,
    ...
  ) %>%
    ungroup() %>%
    rename(.index = .draw) %>%
    select(-c(.chain, .iteration)) %>%
    left_join(params, by = ".index")

  # add CAR(1) process to mean:
  if (type == "epred" & car1) {
    setorderv(preds, order_vars)
    preds <- add_car1(preds, varnames$resp, gr_vars)
  }

  # add CAR(1) residual error:
  if (type == "prediction" & car1) {
    setorderv(preds, order_vars)
    preds <- add_car1_err(preds, car1, gr_vars)
  }

  # add grouping vars:
  preds %>%
    group_by(across(matches(data_vars)))
}

add_resid_draws_car1 <- function(input, object, yvar, ...) {
  add_pred_draws_car1(input, object, ...) %>%
    ungroup() %>%
    mutate(.residual = {{ yvar }} - .epred)
}

# extract response vars from brmsfit:

extract_resp <- function(x) {
  if (class(x)[1] != "brmsfit") stop("'x' must be a brmsfit object")

  prep <- prepare_predictions(x)

  as_char <- as.character(x$formula)
  response <- str_extract(as_char[1], "\\w+")
  cenvar <- str_extract(as_char[1], "(?<=cens\\()\\w+")
  gr_sig <- str_extract(as_char[2], "(?<=sigma ~ ).+(?=\\)$)")
  gr_ar <- prep$dpars$mu$ac$acef$gr
  time_ar <- prep$dpars$mu$ac$acef$time
  list(
    "resp" = response,
    "cens" = cenvar,
    "gr_sigma" = gr_sig,
    "gr_ar" = if_else(gr_ar == "NA", NA_character_, gr_ar),
    "time_ar" = time_ar
  )
}

# ecdf:

ppc_km_nada <- function(input,
                        object,
                        draw_ids = withr::with_seed(1256, {
                          sample(1:4000, 200)
                        }),
                        car1 = TRUE,
                        seed,
                        ...) {
  varnames <- extract_resp(object)

  if (is.na(varnames$cens)) stop("Model formula does not include censoring")

  cenfit_in <- tibble(
    y = pull(input, .data[[varnames$resp]]),
    ci = pull(input, .data[[varnames$cens]]) == "left"
  ) %>%
    filter(!is.na(y))

  cenfit_brms <- NADA::cenfit(
    obs = cenfit_in$y,
    censored = cenfit_in$ci
  )

  preds <- if (car1) {
    withr::with_seed(
      seed,
      {
        add_pred_draws_car1(input, object, draw_ids = draw_ids, type = "prediction", ...)
      }
    )
  } else {
    tidybayes::add_predicted_draws(input, object, draw_ids = draw_ids, seed = seed)
  }

  pp <- preds %>%
    ungroup() %>%
    group_by(.draw) %>%
    nest() %>%
    ungroup() %>%
    mutate(
      cenfit = map(
        data,
        ~ with(.x, NADA::cenfit(.prediction, censored = rep(FALSE, length(.prediction))))
      ),
      cenfit_summ = map(cenfit, summary)
    ) %>%
    unnest(cenfit_summ) %>%
    select(where(~ !is.list(.x)))

  bind_rows(
    "Posterior draws" = pp,
    "Observations" = summary(cenfit_brms),
    .id = "type"
  )
}

# modify stan code to fit CAR(1):

modify_stancode <- function(scode_raw) {
  scode <- scode_raw %>%
    # add time difference variable s:
    str_replace(
      "response variable\\\n",
      "response variable\n  vector[N] s;  // CAR(1) exponent\n"
    ) %>%
    # set lower bound of zero on ar param:
    str_replace(
      "vector\\[Kar\\] ar;",
      "vector<lower=0, upper=1>[Kar] ar;"
    ) %>%
    # convert AR process to CAR1:
    str_replace(
      "mu\\[n\\] \\+= Err\\[n, 1:Kar\\] \\* ar;",
      "mu[n] += Err[n, 1] * pow(ar[1], s[n]); // CAR(1)"
    ) %>%
    # modify prior to reflect 0 lower bound:
    # (n.b., this works for a normal prior on the ar param)
    str_replace(
      "(target.+\\(ar.+\\\n.+normal_lcdf\\()(-1)",
      "\\10"
    )

  class(scode) <- "brmsmodel"

  return(scode)
}

# fit stan model using rstan (checks for CSVs first; if present, load them instead):

fit_stan_model <- function(file,
                           seed,
                           bform,
                           bdata,
                           bpriors = NULL,
                           car1 = TRUE,
                           ...) {
  path <- str_remove(file, "\\/[^\\/]+$")
  csvfiles <- list.files(path = path, pattern = ".+\\.csv", full.names = TRUE)

  regex <- str_extract(file, "[^\\/]+$") %>%
    paste0("_\\d\\.csv")

  # generate stan data:

  data <- brms::make_standata(
    bform,
    data = bdata,
    prior = bpriors,
    family = student()
  )

  if (car1) {
    data$s <- bdata$d_x
  }

  # generate stan code:

  code <- brms::make_stancode(
    bform,
    data = bdata,
    prior = bpriors,
    family = student()
  )

  if (car1) {
    code <- modify_stancode(code)
  }

  # fit model:

  stanmod <- if (sum(str_detect(csvfiles, regex)) > 0) {
    csvfiles[str_detect(csvfiles, regex)] %>%
      rstan::read_stan_csv()
  } else {
    rstan::stan(
      model_code = code,
      data = data,
      sample_file = file, # output in csv format
      seed = seed,
      ...
    )
  }

  # feed back into brms:

  brmsmod <- brm(
    bform,
    data = bdata,
    prior = bpriors,
    family = student(),
    empty = TRUE
  )
  brmsmod$fit <- stanmod
  brmsmod <- rename_pars(brmsmod)

  return(brmsmod)
}

# simulate residuals

impute_censored <- function(x, input, yvar, ycens) {
  x %>%
    ungroup() %>%
    mutate(
      # replace left-censored values with predictions:
      "{yvar}" := if_else(.data[[ycens]] == "left", .prediction, .data[[yvar]]),
      "{ycens}" := "none"
    ) %>%
    select(matches(paste(names(input), collapse = "|")))
}

simulate_residuals <- function(input,
                               object,
                               draw_ids,
                               seed,
                               file,
                               car1 = TRUE,
                               car1_resid = TRUE,
                               model = "gam",
                               ...) {
  if (model != "gam") message("Generating residual draws using tidybayes::add_residual_draws()")

  # extract varnames from fit:

  varnames <- extract_resp(object)

  # generate a draw from the posterior predictive:

  preds_raw <- if (car1) {
    withr::with_seed(
      seed,
      {
        add_pred_draws_car1(
          input = input,
          object = object,
          draw_ids = draw_ids,
          type = "prediction"
        )
      }
    )
  } else {
    tidybayes::add_predicted_draws(input, object, draw_ids = draw_ids, seed = seed)
  }

  model_in_aug <- preds_raw %>%
    impute_censored(input, varnames$resp, varnames$cens)

  brmsmod <- fit_stan_model(
    file, seed, object$formula,
    model_in_aug, object$prior,
    car1 = car1,
    ...
  )

  resids <- if (model == "gam") {
    add_pred_draws_car1(input = input, object = brmsmod, car1 = car1_resid) %>%
      ungroup() %>%
      mutate(.residual = .data[[varnames$resp]] - .epred)
  } else {
    tidybayes::add_residual_draws(
      newdata = input,
      object = object,
      method = "posterior_epred"
    )
  }

  list(
    "residuals" = resids,
    "model" = brmsmod
  )
}

# loo method:

calc_ll <- function(x, response, censored = NULL, cens = TRUE) {
  x <- as.data.table(x)
  if (cens) {
    x[, log_lik := if_else(
      get(censored) == "left",
      pstudent_t(get(response), nu, .epred, sigma, log.p = TRUE),
      dstudent_t(get(response), nu, .epred, sigma, log = TRUE)
    )]
  } else {
    x[, log_lik := dstudent_t(
      get(response), nu, .epred, sigma,
      log = TRUE
    )]
  }
}

loo_cv <- function(input, object, censoring = TRUE, ...) {
  varnames <- extract_resp(object) # extract responses from model formula

  n <- nrow(input) %>%
    seq_len() %>%
    as.character()

  ll <- add_pred_draws_car1(input, object, ...) %>%
    calc_ll(varnames$resp, varnames$cens, cens = censoring) %>%
    dcast(.draw + .chain + .iteration ~ .row, value.var = "log_lik")

  ll_mat <- ll[, c(mget(n))] %>%
    as.matrix()

  rel_eff <- loo::relative_eff(x = exp(ll_mat), chain_id = ll$.chain)

  loo(ll_mat, r_eff = rel_eff)
}
