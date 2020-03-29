#' Simulate results from a single arm Bayesian Phase 2 Trial
#'
#' Trial is a Bayesian phase 2, sequential trial as outlined in
#' Thall, Peter F., and Richard Simon. "Practical Bayesian guidelines for phase
#' IIB clinical trials." Biometrics (1994): 337-349.
#'
#' @param sim_n Number of simulated trials to create
#' @param omega_s_mu Mean of prior distribution on omega_s
#' @param omega_s_width Width of the probability interval (set by `w_conf_level=`)
#' running from the `w_conf_level / 2` to the `1 - w_conf_level / 2` percentiles
#' @param omega_e_c Concentration parameter for prior distribution on omega_e
#' @param delta_0 Targeted improvement of treatment for treatment E over S
#' @param n_min Minimum number of patients that may be enrolled
#' @param n_max Maximum number of patients that may be enrolled
#' @param pr_low Lower probability limit for concluding treatment not promising
#' @param pr_high Upper probablity limit for concluding treatment promising
#' @param omega_s_w_conf_level Confidence level of width of probabilty specified
#' @param mu_e True rate of success with experimental treatment
#' in `omega_s_width`
#' @param verbose When TRUE, additional information is returned as an attribute
#' @param quiet Run with no notes, progress bars, etc.
#' @export

ph2_single_bayes_seq_sim <- function(omega_s_mu, omega_s_width,
                                     omega_s_w_conf_level = 0.90,
                                     omega_e_c, delta_0, n_min, n_max,
                                     sim_n = 1, pr_low = 0.05, pr_high = 0.95,
                                     mu_e = omega_s_mu + delta_0,
                                     verbose = FALSE, quiet = FALSE) {
  # input checks ---------------------------------------------------------------

  # prior distributions --------------------------------------------------------
  # getting parameters for prior distribution of standard treatment
  prior_omega_s <-
    cnvrt_w_to_beta_param(mu = omega_s_mu, w = omega_s_width,
                          w_conf_level = omega_s_w_conf_level)

  # getting parameters for prior distribution of experimental treatment
  prior_omega_e <- list(
    shape1 = omega_e_c * (omega_s_mu + delta_0 / 2),
    shape2 = omega_e_c * (1 - (omega_s_mu + delta_0 / 2) )
  )

  # simulating trial results ---------------------------------------------------
  if (quiet == FALSE) pb <- dplyr::progress_estimated(sim_n) # this is the progress bar
  result <-
    seq_len(sim_n) %>%
    map_dfr(
      function(.x) {
        if (quiet == FALSE) pb$tick()$print()
        smry_trial_result(n_min = n_min, n_max = n_max, mu_e = mu_e,
                        delta_0 = delta_0,
                        prior_omega_s = prior_omega_s,
                        prior_omega_e = prior_omega_e,
                        pr_low = pr_low, pr_high = pr_high) %>%
        as_tibble() %>%
        mutate(sim_id = .x) %>%
        select(.data$sim_id, everything())
      }
    )

  # including additional results if requested
  if (verbose == TRUE) {
    attr(result, "prior_omega_s") <- prior_omega_s
    attr(result, "prior_omega_e") <- prior_omega_e
  }

  result
}






# helper function to summarize trial results
smry_trial_result <- function(n_min, n_max, mu_e,
                              delta_0,
                              prior_omega_s,
                              prior_omega_e,
                              pr_low, pr_high) {


  df_trial_outcome <-
    tibble(
      id = seq_len(n_max),
      outcome = as.numeric(stats::runif(n_max) < mu_e),
      cum_sucess = cumsum(.data$outcome)
    ) %>%
    # decisions only made after n_min enrolled
    filter(.data$id >= n_min) %>%
    mutate(
      lambda_prob = map2_dbl(
        .data$id, .data$cum_sucess,
        ~lambda_fun(n = .x, x = .y, delta_0, prior_omega_s, prior_omega_e)
      ),
      stop_futility = .data$lambda_prob <= .env$pr_low,
      stop_success = .data$lambda_prob >= .env$pr_high
    )


  # returning trial results ----------------------------------------------------
  stop_futility <-
    which(df_trial_outcome$stop_futility) %>%
    {suppressWarnings(min(.))} %>%
    discard(is.infinite)
  stop_success <- which(df_trial_outcome$stop_success) %>%
    {suppressWarnings(min(.))} %>%
    discard(is.infinite)

  # inconclusive
  if (length(stop_futility) == 0 && length(stop_success) == 0) {
    return(list(
      result = "inconclusive",
      n_enrolled = n_max
    ))
  }
  # success
  if (
    (length(stop_futility) == 0 && length(stop_success) > 0) ||
    (length(stop_futility) > 0 && length(stop_success) > 0 && stop_success < stop_futility)
  ) {
    return(list(
      result = "success",
      n_enrolled = df_trial_outcome$id[stop_success]
    ))
  }
  # futility
  if (
    (length(stop_success) == 0 && length(stop_futility) > 0) ||
    (length(stop_success) > 0 && length(stop_futility) > 0 && stop_futility < stop_success)
  ) {
    return(list(
      result = "futile",
      n_enrolled = df_trial_outcome$id[stop_futility]
    ))
  }

  # there should be nothing here
  stop("There was an error")
}

# this is the lambda function from the paper that defines probs of success, futility
lambda_fun <- function(n, x, # number of obs = n, x = successes
                       delta_0, prior_omega_s, prior_omega_e) {
  inside_fun <- function(p) {
    cdf_portion <- stats::pbeta(p + delta_0,
                                shape1 = prior_omega_e$shape1 + x,
                                shape2 = prior_omega_e$shape2 + n - x)

    pdf_portion <- stats::dbeta(p,
                                shape1 = prior_omega_s$shape1,
                                shape2 = prior_omega_s$shape2)

    (1 - cdf_portion) * pdf_portion
  }

  stats::integrate(inside_fun, lower = 0, upper = 1 - delta_0)$value
}

# helper function that converts mu and w into a and b Beta parameterizations
cnvrt_w_to_beta_param <- function(w, w_conf_level, mu) {
  # this function takes an x in (0, 1), converts it to a in (0, Inf),
  # calculates the w, and returns teh difference between the calculated w
  # and our target w
  find_beta_params <- function(x, w, w_conf_level) {
    a = -log(x)
    b = a * (1 - mu) / mu

    w_est <-
      stats::qbeta((1 + w_conf_level) / 2, shape1 = a, shape2 = b) -
      stats::qbeta((1 - w_conf_level) / 2, shape1 = a, shape2 = b)

    abs(w - w_est)
  }

  # this finds the x that results in the w closest our target
  x <- golden_search(
    function(x) find_beta_params(x, w = w, w_conf_level = w_conf_level),
    a = 0, b = 1, tol = 1e-10
  )

  # return the a and b parameters (or the shape1 and shape2)
  list(
    shape1 = -log(x),
    shape2 = -log(x) * (1 - mu) / mu
  )
}


