## ============================================================================
## lmvt() : penalized panel ARX-GARCHX fit
##
## Reference implementation of Algorithm 1 in Wang & Lu (2026), "Variance-Aware
## Penalized Panel Models for Temporal Risk Detection from Wearable Sensor Data".
##
## This implementation is a direct port of `temporalVarGuid.r` from the
## companion demo repository (variance-guided-risk-demo). It performs:
##   (A)  weighted/penalized closed-form mean-block update for
##        {alpha_s, theta, beta} via glmnet (when lambda_beta > 0) or WLS;
##   (B1) joint L-BFGS-B update of (omega_s, a, b) with gamma fixed,
##        evaluated on the strict variance mask (rows with all required
##        lags available);
##   (B2) proximal-gradient update of gamma with backtracking and
##        feasibility repair to keep the variance linear predictor
##        strictly positive;
##   (C)  forward sigma^2 recursion within each subject.
## Per-iteration timing for the mean and variance blocks is recorded.
## ============================================================================

#' Fit a penalized panel ARX-GARCHX model
#'
#' Joint maximum-quasi-likelihood fit of the conditional mean and variance
#' across a balanced panel of subjects, with subject-specific intercepts,
#' shared population coefficients, and L1 penalties on the mean and
#' variance covariate vectors.
#'
#' @param data A data frame containing columns \code{s} (subject ID),
#'   \code{t} (integer time), \code{y} (outcome), and any number of
#'   covariate columns. Other columns (e.g. \code{mu}, \code{sigma},
#'   \code{pi_true}, \code{scenario}) from \code{simulate_scenario()} are
#'   ignored.
#' @param p Autoregressive order in the mean (lags of \eqn{y}).
#' @param q Distributed-lag order for covariates \eqn{x} (used in the
#'   mean equation, and in the variance equation when
#'   \code{use_x_in_variance = TRUE}).
#' @param r ARCH order in the variance (lags of \eqn{e^2}).
#' @param s_ord GARCH order in the variance (lags of \eqn{\sigma^2}).
#' @param lambda_beta L1 penalty for the mean covariate coefficients.
#' @param lambda_gamma L1 penalty for the variance covariate coefficients.
#' @param maxit Maximum number of outer block-coordinate iterations.
#' @param tol Absolute tolerance on the QML objective for convergence.
#' @param standardize_X Logical; if \code{TRUE}, covariates are centered
#'   and scaled (per column) before fitting; means and SDs are stored on
#'   the returned object so \code{predict()} can apply the same
#'   transformation to new data.
#' @param use_x_in_variance Logical; if \code{TRUE}, the same lag block of
#'   \eqn{x} is used in both the mean and the variance equations.
#' @param phi_cap Upper bound on the GARCH stationarity sum
#'   \eqn{\sum a + \sum b}.
#' @param omega_min Lower bound on subject-specific baseline variance.
#' @param omega_cap_mult Multiplier on \code{median(eps^2)} used as the
#'   upper bound for the subject-specific baseline variance \eqn{\omega_s}
#'   in the L-BFGS-B update.
#' @param gamma_steps Number of inner proximal-gradient steps in the
#'   variance covariate (B2) update per outer iteration.
#' @param gamma_step Step size for the proximal-gradient update of
#'   \eqn{\gamma}.
#' @param verbose Logical; print per-iteration diagnostic output.
#'
#' @return An object of class \code{"lmvt"}: a list with elements
#'   \code{alpha} (subject intercepts), \code{theta} (AR coefficients),
#'   \code{beta} (mean covariate coefficients on the lagged-X block),
#'   \code{omega}, \code{a}, \code{b} (GARCH parameters), \code{gamma}
#'   (variance covariate coefficients), \code{ids}, \code{xcols},
#'   \code{x_center}, \code{x_scale}, \code{orders}, \code{penalties},
#'   \code{call}, plus diagnostics \code{converged}, \code{iters},
#'   \code{mean_times}, \code{var_times}, \code{iter_times}.
#'
#' @references Wang Z. & Lu M. (2026). Variance-Aware Penalized Panel Models for
#'   Temporal Risk Detection from Wearable Sensor Data.
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_scenario(scen = 1, S = 5, T = 80, d_noise = 5,
#'                          noise_kind = "iid", seed = 1)
#' df  <- sim[, c("s", "t", "y", grep("^X|^Noise|^Xbin",
#'                                    names(sim), value = TRUE))]
#' fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
#'             lambda_beta = 0.05, lambda_gamma = 0.05)
#' print(fit)
#' @export
lmvt <- function(data,
                 p = 1, q = 0,
                 r = 1, s_ord = 1,
                 lambda_beta  = 0,
                 lambda_gamma = 0,
                 maxit = 200, tol = 1e-4,
                 standardize_X     = TRUE,
                 use_x_in_variance = TRUE,
                 phi_cap     = 0.995,
                 omega_min   = 1e-6,
                 omega_cap_mult = 10,
                 gamma_steps = 8,
                 gamma_step  = 1e-3,
                 verbose     = FALSE) {

    cl <- match.call()
    stopifnot(all(c("s", "t", "y") %in% names(data)))
    data <- data[order(data$s, data$t), , drop = FALSE]
    ids  <- sort(unique(data$s))
    S    <- length(ids)
    n    <- nrow(data)

    # autodetect covariates: drop simulator-side oracle columns if present
    keep_cols <- c("s", "t", "y", "mu", "sigma", "pi_true", "scenario")
    xcols <- setdiff(names(data), keep_cols)
    if (length(xcols) == 0L)
        stop("No covariates: provide columns beyond s, t, y.")
    Xraw <- data.matrix(data[, xcols, drop = FALSE])

    if (standardize_X) {
        x_center <- colMeans(Xraw, na.rm = TRUE)
        x_scale  <- apply(Xraw, 2, stats::sd, na.rm = TRUE)
        x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
        X <- scale(Xraw, center = x_center, scale = x_scale)
        attr(X, "scaled:center") <- NULL
        attr(X, "scaled:scale")  <- NULL
    } else {
        x_center <- rep(0, ncol(Xraw))
        x_scale  <- rep(1, ncol(Xraw))
        X <- Xraw
    }
    P <- ncol(X)
    Y <- as.numeric(data$y)

    # subject-aware lag of a numeric vector v by L positions within each id
    lag_by_id <- function(v, L) {
        out <- rep(NA_real_, length(v))
        if (L <= 0L) return(v)
        for (id in ids) {
            idx <- which(data$s == id)
            if (length(idx) > L) {
                out[idx[(L + 1L):length(idx)]] <-
                    v[idx[1:(length(idx) - L)]]
            }
        }
        out
    }

    # mean-side design: AR lags of y plus distributed lag block of X (lag 0..q)
    Ylags <- if (p > 0L)
        do.call(cbind, lapply(seq_len(p), function(L) lag_by_id(Y, L)))
    else NULL
    Xlags_list_mean <- lapply(0:q, function(L)
        if (L == 0L) X else apply(X, 2, lag_by_id, L = L))
    Xmean <- do.call(cbind, Xlags_list_mean)

    # variance-side design follows training choice
    if (isTRUE(use_x_in_variance)) {
        Xlags_list_var <- Xlags_list_mean
        Xvar  <- Xmean
        gamma <- rep(0, (q + 1L) * P)
    } else {
        Xlags_list_var <- list(matrix(0, n, 0))
        Xvar  <- matrix(0, n, 0)
        gamma <- numeric(0)
        lambda_gamma <- 0
    }

    # mean-block validity mask (rows with all required mean lags present)
    mask_mean <- rep(TRUE, n)
    if (p > 0L) mask_mean <- mask_mean & (rowSums(is.na(Ylags)) == 0)
    for (LL in 0:q)
        mask_mean <- mask_mean & (rowSums(is.na(Xlags_list_mean[[LL + 1L]])) == 0)

    # initial parameter values
    vl <- stats::var(Y, na.rm = TRUE)
    vl <- if (is.finite(vl) && vl > 0) vl else 1
    alpha <- tapply(Y, data$s, function(z) mean(z, na.rm = TRUE))
    alpha[is.na(alpha)] <- 0
    alpha <- as.numeric(alpha); names(alpha) <- ids
    omega <- rep(max(omega_min, vl / 10), S); names(omega) <- ids
    theta <- if (p > 0L) rep(0, p) else numeric(0)
    beta  <- rep(0, (q + 1L) * P)
    a_vec <- if (r > 0L)     rep(0.05, r)                                  else numeric(0)
    b_vec <- if (s_ord > 0L) rep((phi_cap - 0.1) / max(1L, s_ord), s_ord) else numeric(0)

    # initial residuals and sigma^2
    mu     <- alpha[match(data$s, ids)] + as.numeric(Xmean %*% beta)
    eps    <- Y - mu
    eps[!is.finite(eps)] <- 0
    sigma2 <- rep(max(vl, 1e-2), n)

    make_var_lags <- function(eps_v, sigma2_v) {
        E2lags <- if (r > 0L)
            do.call(cbind, lapply(seq_len(r), function(L) lag_by_id(eps_v^2, L)))
        else NULL
        S2lags <- if (s_ord > 0L)
            do.call(cbind, lapply(seq_len(s_ord), function(L) lag_by_id(sigma2_v, L)))
        else NULL
        list(E2lags = E2lags, S2lags = S2lags)
    }

    proj_ab <- function(a, b, cap = phi_cap) {
        a <- pmax(a, 0); b <- pmax(b, 0)
        s <- sum(a) + sum(b)
        if (s > cap && s > 0) {
            sc <- cap / s
            a <- a * sc; b <- b * sc
        }
        list(a = a, b = b)
    }

    # safe objective for (omega, a, b) with gamma fixed; evaluated only on valid_idx
    var_obj_abw <- function(par, eps, ids, s_idx, E2lags, S2lags, Xvar,
                            gamma, phi_cap, omega_min, omega_max, valid_idx) {
        BIG  <- 1e50
        Sloc <- length(unique(ids))
        n_a  <- if (is.null(E2lags)) 0L else ncol(E2lags)
        n_b  <- if (is.null(S2lags)) 0L else ncol(S2lags)

        omega <- par[1:Sloc]
        a <- if (n_a > 0L) par[(Sloc + 1L):(Sloc + n_a)] else numeric(0)
        b <- if (n_b > 0L) par[(Sloc + n_a + 1L):(Sloc + n_a + n_b)] else numeric(0)

        if (any(!is.finite(omega))) return(BIG)
        omega <- pmin(pmax(omega, omega_min), omega_max)
        a[!is.finite(a)] <- 0; b[!is.finite(b)] <- 0
        if (any(a < 0) || any(b < 0)) return(BIG)
        ssum <- sum(a) + sum(b)
        if (ssum > phi_cap && ssum > 0) { sc <- phi_cap / ssum; a <- a * sc; b <- b * sc }

        delta <- omega[s_idx]
        if (n_a > 0L) {
            if (anyNA(E2lags[valid_idx, , drop = FALSE])) return(BIG)
            delta <- delta + as.numeric(E2lags %*% a)
        }
        if (n_b > 0L) {
            if (anyNA(S2lags[valid_idx, , drop = FALSE])) return(BIG)
            delta <- delta + as.numeric(S2lags %*% b)
        }
        u <- delta + if (length(gamma)) as.numeric(Xvar %*% gamma) else 0
        v <- valid_idx
        if (any(!is.finite(u[v]))) return(BIG)
        if (any(u[v] <= 0))        return(BIG)

        val <- sum(log(u[v]) + (eps[v]^2) / u[v])
        if (!is.finite(val)) return(BIG)
        val
    }

    # per-iteration timing
    mean_times <- numeric(maxit)
    var_times  <- numeric(maxit)
    iter_times <- numeric(maxit)

    obj_prev <- Inf
    for (k in seq_len(maxit)) {

        iter_start <- proc.time()[["elapsed"]]

        ## ===== (A) Mean block ===================================================
        mean_start <- proc.time()[["elapsed"]]

        w        <- 1 / pmax(sigma2, 1e-10)
        AR_part  <- if (p > 0L) rowSums(Ylags %*% matrix(theta, ncol = 1)) else 0
        X_part   <- as.numeric(Xmean %*% beta)
        numer    <- tapply(w * (Y - AR_part - X_part), data$s, sum, na.rm = TRUE)
        denom    <- tapply(w, data$s, sum, na.rm = TRUE)
        alpha    <- as.numeric(numer / pmax(denom, 1e-12))
        names(alpha) <- ids
        alpha[!is.finite(alpha)] <- 0

        X_A   <- if (p > 0L) cbind(Ylags, Xmean) else Xmean
        y_A   <- Y - alpha[match(data$s, ids)]
        keepA <- mask_mean & is.finite(y_A) & is.finite(w) &
            is.finite(rowSums(if (is.null(dim(X_A))) cbind(X_A) else X_A))
        X_Ak  <- X_A[keepA, , drop = FALSE]
        y_Ak  <- y_A[keepA]
        wk    <- w[keepA]
        Xw    <- X_Ak * sqrt(wk)
        yw    <- y_Ak * sqrt(wk)

        if (lambda_beta > 0) {
            if (!requireNamespace("glmnet", quietly = TRUE))
                stop("Install 'glmnet' for penalized mean block (lambda_beta > 0).")
            pen_vec <- c(rep(0, p), rep(1, ncol(Xmean)))
            gfit <- glmnet::glmnet(Xw, yw, family = "gaussian",
                                   alpha = 1, lambda = lambda_beta,
                                   penalty.factor = pen_vec,
                                   standardize = FALSE, intercept = FALSE)
            coef_all <- as.numeric(as.matrix(glmnet::coef.glmnet(gfit)))[-1]
        } else {
            XtX <- crossprod(Xw); Xty <- crossprod(Xw, yw)
            coef_all <- tryCatch(
                as.numeric(solve(XtX, Xty)),
                error = function(e)
                    as.numeric(solve(XtX + 1e-8 * diag(ncol(XtX)), Xty))
            )
        }
        if (p > 0L) {
            theta <- coef_all[1:p]
            beta  <- coef_all[(p + 1L):length(coef_all)]
        } else {
            beta <- coef_all
        }

        AR_part <- if (p > 0L) rowSums(Ylags %*% matrix(theta, ncol = 1)) else 0
        X_part  <- as.numeric(Xmean %*% beta)
        mu      <- alpha[match(data$s, ids)] + AR_part + X_part
        eps     <- Y - mu
        eps[!is.finite(eps)] <- 0

        mean_times[k] <- proc.time()[["elapsed"]] - mean_start

        ## ===== Build strict variance mask =======================================
        var_start <- proc.time()[["elapsed"]]

        VL <- make_var_lags(eps, sigma2)
        E2lags <- VL$E2lags; S2lags <- VL$S2lags
        mask_var <- rep(TRUE, n)
        if (r > 0L)     mask_var <- mask_var & (rowSums(is.na(E2lags)) == 0)
        if (s_ord > 0L) mask_var <- mask_var & (rowSums(is.na(S2lags)) == 0)
        if (use_x_in_variance) {
            for (LL in 0:q)
                mask_var <- mask_var & (rowSums(is.na(Xlags_list_var[[LL + 1L]])) == 0)
        }
        valid_idx <- which(mask_var & is.finite(eps))
        if (!length(valid_idx))
            stop("No valid rows for variance update; check lags/orders.")

        ## ===== (B1) Optimize (omega, a, b) with gamma fixed =====================
        s_idx <- match(data$s, ids)
        omega_cap <- max(omega_cap_mult * stats::median(eps^2, na.rm = TRUE),
                         omega_min * 10)
        par0 <- c(omega, if (r > 0L) a_vec, if (s_ord > 0L) b_vec)
        lower <- c(rep(omega_min, S),
                   if (r > 0L)     rep(0, r),
                   if (s_ord > 0L) rep(0, s_ord))
        upper <- c(rep(omega_cap, S),
                   if (r > 0L)     rep(phi_cap, r),
                   if (s_ord > 0L) rep(phi_cap, s_ord))

        opt <- stats::optim(par0, var_obj_abw, method = "L-BFGS-B",
                            lower = lower, upper = upper,
                            control = list(maxit = 200),
                            eps = eps, ids = data$s, s_idx = s_idx,
                            E2lags = E2lags, S2lags = S2lags, Xvar = Xvar,
                            gamma = gamma, phi_cap = phi_cap,
                            omega_min = omega_min, omega_max = omega_cap,
                            valid_idx = valid_idx)

        par_star <- opt$par
        omega    <- par_star[1:S]
        if (r > 0L)     a_vec <- par_star[(S + 1L):(S + r)]
        if (s_ord > 0L) {
            off   <- S + if (r > 0L) r else 0L
            b_vec <- par_star[(off + 1L):(off + s_ord)]
        }
        abp <- proj_ab(a_vec, b_vec, phi_cap)
        a_vec <- abp$a; b_vec <- abp$b

        ## ===== (B2) Update gamma with delta fixed (prox-grad on valid_idx) =====
        delta <- omega[s_idx] +
            if (r > 0L)     as.numeric(E2lags %*% a_vec) else 0 +
            if (s_ord > 0L) as.numeric(S2lags %*% b_vec) else 0

        if (use_x_in_variance && length(gamma)) {
            v <- valid_idx
            u_floor <- 1e-8

            for (gi in seq_len(max(1L, gamma_steps))) {
                u  <- delta + as.numeric(Xvar %*% gamma)
                uv <- u[v]
                step <- gamma_step

                for (bt in 0:8) {
                    if (any(!is.finite(uv)) || any(uv <= u_floor)) {
                        shrink <- 1.0
                        repeat {
                            gamma_test <- gamma * shrink
                            uv_test <- delta[v] + as.numeric(Xvar[v, ] %*% gamma_test)
                            if (all(is.finite(uv_test)) && all(uv_test > u_floor)) {
                                gamma <- gamma_test
                                uv    <- uv_test
                                break
                            }
                            shrink <- shrink * 0.5
                            if (shrink < 1e-6) {
                                gamma <- 0 * gamma
                                uv    <- delta[v]
                                break
                            }
                        }
                        step <- step * 0.5
                    }

                    g_u      <- (1 / uv - (eps[v]^2) / (uv^2))
                    grad_gam <- as.numeric(crossprod(Xvar[v, , drop = FALSE], g_u))
                    gamma_try <- gamma - step * grad_gam
                    if (lambda_gamma > 0)
                        gamma_try <- sign(gamma_try) *
                            pmax(abs(gamma_try) - step * lambda_gamma, 0)

                    u_try <- delta[v] + as.numeric(Xvar[v, ] %*% gamma_try)
                    if (any(!is.finite(u_try)) || any(u_try <= u_floor)) {
                        step <- step * 0.5; next
                    }

                    obj_now <- sum(log(uv) + (eps[v]^2) / uv)
                    obj_try <- sum(log(u_try) + (eps[v]^2) / u_try)
                    if (is.finite(obj_try) && obj_try <= obj_now) {
                        gamma <- gamma_try; break
                    } else {
                        step <- step * 0.5
                    }
                }
            }

            uv_final <- delta[v] + as.numeric(Xvar[v, ] %*% gamma)
            if (any(!is.finite(uv_final)) || any(uv_final <= u_floor)) {
                shrink <- 1.0
                repeat {
                    gamma_test <- gamma * shrink
                    uv_test <- delta[v] + as.numeric(Xvar[v, ] %*% gamma_test)
                    if (all(is.finite(uv_test)) && all(uv_test > u_floor)) {
                        gamma <- gamma_test; break
                    }
                    shrink <- shrink * 0.5
                    if (shrink < 1e-6) { gamma <- 0 * gamma; break }
                }
            }
        } else {
            gamma <- numeric(0)
        }

        ## ===== (C) Forward sigma^2 recursion ===================================
        sigma2_new <- rep(NA_real_, n)
        for (id in ids) {
            idx <- which(data$s == id)
            for (ii in seq_along(idx)) {
                tt <- idx[ii]
                arch_part <- 0; garch_part <- 0
                if (r > 0L) {
                    for (L in seq_len(r)) {
                        tprev <- which(data$t[idx] == (data$t[tt] - L))
                        if (length(tprev) == 1L)
                            arch_part <- arch_part + a_vec[L] * (eps[idx[tprev]]^2)
                    }
                }
                if (s_ord > 0L) {
                    for (L in seq_len(s_ord)) {
                        tprev <- which(data$t[idx] == (data$t[tt] - L))
                        if (length(tprev) == 1L)
                            garch_part <- garch_part + b_vec[L] * sigma2_new[idx[tprev]]
                    }
                }
                xrow_var <- if (ncol(Xvar)) Xvar[tt, ] else numeric(0)
                vlin <- omega[match(id, ids)] + arch_part + garch_part +
                    if (length(gamma)) sum(gamma * xrow_var) else 0
                sigma2_new[tt] <- max(vlin, 1e-8)
            }
        }
        sigma2 <- sigma2_new

        var_times[k] <- proc.time()[["elapsed"]] - var_start

        # monitor (valid_idx only)
        u_now <- delta[valid_idx] +
            if (length(gamma)) as.numeric(Xvar[valid_idx, ] %*% gamma) else 0
        obj_now <- sum(log(u_now) + (eps[valid_idx]^2) / u_now)

        iter_times[k] <- proc.time()[["elapsed"]] - iter_start

        if (verbose)
            cat(sprintf(
                "iter %d: obj=%s  sum(a)+sum(b)=%.4f  mean(omega)=%.4f  t_mean=%.4fs  t_var=%.4fs\n",
                k, format(obj_now, digits = 6), sum(a_vec) + sum(b_vec),
                mean(omega), mean_times[k], var_times[k]))

        if (!is.finite(obj_now)) break
        if (abs(obj_prev - obj_now) < tol) break
        obj_prev <- obj_now
    }

    iters_used <- k
    mean_times <- mean_times[seq_len(iters_used)]
    var_times  <- var_times[seq_len(iters_used)]
    iter_times <- iter_times[seq_len(iters_used)]

    structure(list(
        call          = cl,
        orders        = list(p = p, q = q, r = r, s = s_ord),
        penalties     = list(lambda_beta = lambda_beta,
                             lambda_gamma = lambda_gamma),
        ids           = ids,
        alpha         = alpha,
        theta         = theta,
        beta          = beta,
        omega         = omega,
        a             = a_vec,
        b             = b_vec,
        gamma         = gamma,
        xcols         = xcols,
        standardize_X = standardize_X,
        use_x_in_variance = use_x_in_variance,
        x_center      = x_center,
        x_scale       = x_scale,
        converged     = (iters_used < maxit),
        iters         = iters_used,
        mean_times    = mean_times,
        var_times     = var_times,
        iter_times    = iter_times
    ), class = "lmvt")
}

# ---- S3 helpers --------------------------------------------------------------

#' @export
print.lmvt <- function(x, ...) {
    cat("Penalized panel ARX-GARCHX fit (varGuidTS::lmvt)\n")
    cat(sprintf("  Subjects:           %d\n", length(x$ids)))
    cat(sprintf("  Mean lags  (p, q):  (%d, %d)\n", x$orders$p, x$orders$q))
    cat(sprintf("  Var  lags  (r, s):  (%d, %d)\n", x$orders$r, x$orders$s))
    cat(sprintf("  Penalties (lambda_beta, lambda_gamma): (%.3g, %.3g)\n",
                x$penalties$lambda_beta, x$penalties$lambda_gamma))
    cat(sprintf("  Iterations / converged: %d / %s\n",
                x$iters, isTRUE(x$converged)))
    cat(sprintf("  Nonzero beta:  %d / %d\n",
                sum(abs(x$beta) > 1e-8), length(x$beta)))
    cat(sprintf("  Nonzero gamma: %d / %d\n",
                sum(abs(x$gamma) > 1e-8), length(x$gamma)))
    if (length(x$iter_times))
        cat(sprintf("  Mean time/iter (mean / var / total): %.4f / %.4f / %.4fs\n",
                    mean(x$mean_times), mean(x$var_times),
                    mean(x$iter_times)))
    invisible(x)
}

#' @export
coef.lmvt <- function(object, ...) {
    list(alpha = object$alpha, theta = object$theta,
         beta  = object$beta,  omega = object$omega,
         a     = object$a,     b     = object$b,
         gamma = object$gamma)
}
