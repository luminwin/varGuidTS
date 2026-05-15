## ============================================================================
## simulate_scenario() : data generators for the four canonical scenarios
## used in the simulation studies of Wang & Lu (2026), Section III.
##
## Ported from `simulator_scenarios.r` in the demo repository, with
## ggplot2/dplyr dependencies removed (plotting helpers live in
## inst/scripts/ instead).
## ============================================================================

# Standardized Student-t (unit variance)
.rstdt <- function(n, nu) {
    stopifnot(nu > 2)
    stats::rt(n, df = nu) / sqrt(nu / (nu - 2))
}
.pstdt <- function(q, nu) {
    stopifnot(nu > 2)
    stats::pt(q * sqrt(nu / (nu - 2)), df = nu)
}

# Piecewise regime path around a target run length L.
.gen_regime_piecewise <- function(Ttot, L,
                                  min_run = max(1L, floor(0.5 * L)),
                                  jitter  = floor(0.3 * L)) {
    reg <- integer(Ttot); t <- 1L; state <- 0L
    while (t <= Ttot) {
        run <- max(min_run,
                   L + sample.int(2L * jitter + 1L, 1L) - (jitter + 1L))
        t_to <- min(Ttot, t + run - 1L)
        reg[t:t_to] <- state
        state <- 1L - state
        t <- t_to + 1L
    }
    reg
}

# Covariate generators
.gen_iid <- function(S, Ttot, d) {
    array(stats::rnorm(S * Ttot * d), c(S, Ttot, d))
}
.gen_ar1_cov <- function(S, Ttot, d, phi = 0.6, sigma_eta = 1) {
    phi <- rep(phi, length.out = d)
    sigma_eta <- rep(sigma_eta, length.out = d)
    arr <- array(0, c(S, Ttot, d))
    for (k in seq_len(d)) {
        sd0 <- sqrt(sigma_eta[k]^2 / (1 - phi[k]^2))
        arr[, 1, k] <- stats::rnorm(S, 0, sd0)
        for (t in 2:Ttot)
            arr[, t, k] <- phi[k] * arr[, t - 1, k] +
                stats::rnorm(S, 0, sigma_eta[k])
    }
    arr
}
.gen_seasonal <- function(S, Ttot, d, period) {
    stopifnot(d >= 2)
    tt <- seq_len(Ttot)
    s1 <- sin(2 * pi * tt / period); c1 <- cos(2 * pi * tt / period)
    arr <- array(0, c(S, Ttot, d))
    arr[, , 1] <- matrix(rep(s1, each = S), S, Ttot)
    arr[, , 2] <- matrix(rep(c1, each = S), S, Ttot)
    if (d > 2) arr[, , 3:d] <- stats::rnorm(S * Ttot * (d - 2))
    arr
}

.auto_period <- function(Ttot) max(12, round(Ttot / 3))
.slice_X_time <- function(X, t) {
    S <- dim(X)[1]; d <- dim(X)[3]
    Xt <- X[, t, , drop = FALSE]; dim(Xt) <- c(S, d); Xt
}
.slice_X_cov_keep <- function(X, keep, j) {
    S <- dim(X)[1]; K <- length(keep)
    Xjk <- X[, keep, j, drop = FALSE]
    dim(Xjk) <- c(S, K)
    Xjk
}

# Internal configuration for the four canonical simulation scenarios.
.scenario_config <- function(scen) {
    stopifnot(scen %in% c(1L, 2L, 3L, 4L))
    switch(as.character(scen),
           "1" = list(
               name = "Baseline",
               theta = 0.4, a = 0.05, b = 0.50,
               innov = list(dist = "norm", nu = NA),
               alpha = list(kind = "normal", mean = 0, sd = 0.5),
               omega = list(kind = "const",  val  = 0.10),
               cov_kind = list(kind = "iid"),
               burn = 200, use_bin = FALSE, bin_def = NULL,
               beta_active  = c(X1 = 0.40),
               gamma_active = c(X1 = 0.20)
           ),
           "2" = list(
               name = "Mean-driven",
               theta = 0.5, a = 0.06, b = 0.50,
               innov = list(dist = "norm", nu = NA),
               alpha = list(kind = "normal", mean = 0, sd = 0.5),
               omega = list(kind = "const",  val  = 0.10),
               cov_kind = list(kind = "seasonal"),
               burn = 220, use_bin = FALSE, bin_def = NULL,
               beta_active  = c(X1 = 0.80, X2 = 0.60),
               gamma_active = c(X1 = 0.10)
           ),
           "3" = list(
               name = "Regime Contrast (Xbin via subject median)",
               theta = 0.40, a = 0.06, b = 0.20,
               innov = list(dist = "norm", nu = NA),
               alpha = list(kind = "normal", mean = 0, sd = 0.5),
               omega = list(kind = "const",  val  = 0.01),
               cov_kind = list(kind = "ar1"),
               burn = 180, use_bin = TRUE, bin_def = "median",
               beta_active  = c(X1 = 0.20, X2 = 0.10),
               gamma_active = c(Xbin = 0.40)
           ),
           "4" = list(
               name = "Seasonal + Piecewise Regime (Xbin via L)",
               theta = 0.40, a = 0.06, b = 0.20,
               innov = list(dist = "stdt", nu = 8),
               alpha = list(kind = "normal", mean = 0, sd = 0.5),
               omega = list(kind = "const",  val  = 0.01),
               cov_kind = list(kind = "seasonal"),
               burn = 180, use_bin = TRUE, bin_def = "piecewise",
               bin_params = list(min_run_frac = 0.5, jitter_frac = 0.3),
               beta_active  = c(X1 = 0.80, X2 = 0.60),
               gamma_active = c(Xbin = 0.50)
           )
    )
}

.draw_alpha <- function(S, alpha) {
    if (alpha$kind == "normal")
        stats::rnorm(S, alpha$mean, alpha$sd)
    else stop("alpha$kind must be 'normal'")
}
.draw_omega <- function(S, omega) {
    if (omega$kind == "const") rep(omega$val, S)
    else stop("omega$kind must be 'const'")
}

#' Simulate a panel ARX-GARCHX dataset under one of four scenarios
#'
#' @param scen Scenario id in 1:4.
#' @param S Number of subjects.
#' @param T Number of post-burn-in time points (per subject).
#' @param L Mean run length for the piecewise regime (Scenario 4).
#' @param d_noise Number of nuisance covariates with zero true effect.
#' @param noise_kind Either \code{"iid"} or \code{"ar1"}.
#' @param threshold_c Threshold \eqn{c} at which the oracle exceedance
#'   probability \code{pi_true} is computed and stored.
#' @param seed Random seed.
#' @param floor_sigma2 Lower clamp on the simulated variance.
#' @return A long-format data frame with columns \code{scenario},
#'   \code{s}, \code{t}, \code{y}, oracle \code{mu} and \code{sigma},
#'   the covariate columns \code{X1..X3}, optional \code{Noise1..Noiseq},
#'   optional \code{Xbin}, and the oracle exceedance \code{pi_true}.
#'   Generating parameters are stored in \code{attr(., "params")}.
#'
#' @examples
#' set.seed(7)
#' sim <- simulate_scenario(scen = 2, S = 5, T = 80,
#'                          d_noise = 5, noise_kind = "iid", seed = 7)
#' head(sim)
#' @export
simulate_scenario <- function(
    scen        = 1L,
    S           = 8L,
    T           = 240L,
    L           = 20L,
    d_noise     = 20L,
    noise_kind  = c("iid", "ar1"),
    threshold_c = 1.0,
    seed        = NULL,
    floor_sigma2 = 1e-10
) {
    noise_kind <- match.arg(noise_kind)
    if (!is.null(seed)) set.seed(seed)
    cfg  <- .scenario_config(scen)
    Ttot <- T + cfg$burn

    X_base <- switch(cfg$cov_kind$kind,
        "iid"      = .gen_iid(S, Ttot, 3),
        "ar1"      = .gen_ar1_cov(S, Ttot, 3, phi = 0.6, sigma_eta = 1),
        "seasonal" = .gen_seasonal(S, Ttot, 3,
                                   period = .auto_period(Ttot)),
        stop("Unknown cov_kind")
    )

    X_noise <- if (d_noise > 0)
        if (noise_kind == "iid") .gen_iid(S, Ttot, d_noise)
        else .gen_ar1_cov(S, Ttot, d_noise, phi = 0.5, sigma_eta = 1)
    else array(0, c(S, Ttot, 0))

    feat_names <- c(paste0("X", 1:3),
                    if (d_noise > 0) paste0("Noise", seq_len(d_noise)))
    if (isTRUE(cfg$use_bin)) feat_names <- c(feat_names, "Xbin")
    d_total <- length(feat_names)

    X <- array(0, c(S, Ttot, d_total))
    X[, , 1:3] <- X_base
    if (d_noise > 0) X[, , 4:(3 + d_noise)] <- X_noise

    if (isTRUE(cfg$use_bin)) {
        if ((cfg$bin_def %||% "median") == "median") {
            keep_idx <- (cfg$burn + 1):Ttot
            X1_keep  <- X_base[, keep_idx, 1, drop = FALSE]
            med_s    <- apply(X1_keep, 1, stats::median, na.rm = TRUE)
            med_mat  <- matrix(med_s, nrow = S, ncol = Ttot, byrow = FALSE)
            Xbin     <- 1L * (X_base[, , 1] > med_mat)
            X[, , d_total] <- Xbin
        } else if (cfg$bin_def == "piecewise") {
            min_run <- round((cfg$bin_params$min_run_frac %||% 0.5) * L)
            jitter  <- round((cfg$bin_params$jitter_frac  %||% 0.3) * L)
            Xbin    <- .gen_regime_piecewise(Ttot, L = L,
                                             min_run = min_run,
                                             jitter  = jitter)
            X[, , d_total] <- matrix(rep(Xbin, each = S), S, Ttot)
        } else stop("Unknown bin_def: ", cfg$bin_def)
    }

    beta  <- stats::setNames(rep(0, d_total), feat_names)
    gamma <- stats::setNames(rep(0, d_total), feat_names)
    set_actives <- function(vec, act_map) {
        for (nm in names(act_map))
            if (nm %in% names(vec)) vec[nm] <- act_map[[nm]]
        vec
    }
    beta  <- set_actives(beta,  cfg$beta_active)
    gamma <- set_actives(gamma, cfg$gamma_active)

    alpha_s <- .draw_alpha(S, cfg$alpha)
    omega_s <- .draw_omega(S, cfg$omega)

    y <- mu <- sig2 <- e <- matrix(0, S, Ttot)
    for (s in seq_len(S)) {
        sig2_init <- if (cfg$a + cfg$b < 1)
            omega_s[s] / (1 - cfg$a - cfg$b)
        else omega_s[s]
        sig2[s, 1] <- sig2_init
        y[s, 1]    <- alpha_s[s]
        mu[s, 1]   <- alpha_s[s]
        e[s, 1]    <- 0
    }
    for (t in 2:Ttot) {
        Xt   <- .slice_X_time(X, t)
        mu_x <- as.numeric(Xt %*% beta)
        vx   <- as.numeric(Xt %*% gamma)
        mu[, t]   <- alpha_s + cfg$theta * y[, t - 1] + mu_x
        sig2[, t] <- omega_s + cfg$a * (e[, t - 1]^2) +
                     cfg$b * sig2[, t - 1] + vx
        sig2[, t] <- pmax(sig2[, t], floor_sigma2)
        eps_t <- if (cfg$innov$dist == "norm")
            stats::rnorm(S) else .rstdt(S, cfg$innov$nu)
        y[, t] <- mu[, t] + sqrt(sig2[, t]) * eps_t
        e[, t] <- y[, t] - mu[, t]
    }

    keep <- (cfg$burn + 1):Ttot
    df <- data.frame(
        scenario = cfg$name,
        s = rep(seq_len(S), each = length(keep)),
        t = rep(seq_along(keep), times = S),
        y = as.vector(t(y[, keep])),
        mu = as.vector(t(mu[, keep])),
        sigma = sqrt(as.vector(t(sig2[, keep])))
    )
    for (j in seq_along(feat_names))
        df[[feat_names[j]]] <- as.vector(t(.slice_X_cov_keep(X, keep, j)))

    z <- (threshold_c - df$mu) / df$sigma
    df$pi_true <- if (cfg$innov$dist == "norm")
        (1 - stats::pnorm(z))
    else
        (1 - .pstdt(z, cfg$innov$nu))

    attr(df, "params") <- list(
        cfg = cfg, beta = beta, gamma = gamma,
        threshold = threshold_c, d_noise = d_noise,
        noise_kind = noise_kind, feat_names = feat_names, L = L
    )
    df
}
