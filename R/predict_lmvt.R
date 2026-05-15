## ============================================================================
## predict.lmvt() : conditional mean, variance, and exceedance probability
##
## Reference implementation ported from `temporalVarGuid.r` in the
## variance-guided-risk-demo repository, aligned with the package's lmvt()
## API. The forward variance recursion runs subject-by-subject and
## includes a feasibility-repair step that mirrors the training-time
## guards in the (B2) gamma update, so that the variance linear predictor
## is strictly positive even on out-of-sample rows.
## ============================================================================

#' Predict from a fitted \code{lmvt} object
#'
#' Computes the conditional mean \eqn{\hat\mu_{s,t}}, conditional standard
#' deviation \eqn{\hat\sigma_{s,t}}, and the exceedance-based risk
#' score \eqn{\hat\pi_{s,t}(c) = P(Y_{s,t} > c \,|\, F_{t-1})} on a new
#' (or held-out) panel data frame.
#'
#' @param object A fitted \code{lmvt} object.
#' @param newdata A data frame containing \code{s}, \code{t}, \code{y}
#'   columns and the same covariate columns used during fitting.
#' @param threshold Threshold \eqn{c} at which to evaluate the exceedance
#'   probability. May be (a) a single numeric, (b) a vector of length
#'   \code{nrow(newdata)}, or (c) a per-subject vector (named or aligned
#'   with the unique subjects in \code{newdata}). Defaults to the pooled
#'   90th percentile of \code{newdata$y}.
#' @param innov_g Logical; return Gaussian-innovation exceedance
#'   probability \code{risk_g}.
#' @param innov_t Logical; also return Student-t-innovation exceedance
#'   probability \code{risk_t}.
#' @param df_t Degrees of freedom for the standardized Student-t
#'   innovation distribution. Must be greater than 2 when
#'   \code{innov_t = TRUE}.
#' @param ... Unused, kept for S3 compatibility.
#'
#' @return A data frame with columns \code{s}, \code{t},
#'   \code{yhat = muhat}, \code{muhat}, \code{sigmahat},
#'   \code{threshold_used}, and one or both of \code{risk_g}, \code{risk_t}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- simulate_scenario(scen = 2, S = 4, T = 70, d_noise = 4,
#'                          noise_kind = "iid", seed = 1)
#' df  <- sim[, c("s", "t", "y", grep("^X|^Noise|^Xbin",
#'                                    names(sim), value = TRUE))]
#' fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
#'             lambda_beta = 0.05, lambda_gamma = 0.05, maxit = 8)
#'
#' # 1) Pooled/global threshold: one cutoff is applied to every row.
#' pr_global <- predict(fit, df, threshold = stats::quantile(df$y, 0.90))
#' head(pr_global)
#'
#' # 2) Subject-specific thresholds: use a named vector indexed by subject ID.
#' subject_cutoffs <- tapply(df$y, df$s, stats::quantile, probs = 0.75,
#'                           na.rm = TRUE)
#' pr_subject <- predict(fit, df, threshold = subject_cutoffs)
#' head(pr_subject)
#'
#' # 3) Return both Gaussian and Student-t exceedance risk scores.
#' pr_t <- predict(fit, df, threshold = subject_cutoffs,
#'                 innov_g = TRUE, innov_t = TRUE, df_t = 8)
#' head(pr_t[, c("s", "t", "muhat", "sigmahat", "risk_g", "risk_t")])
#' }
#' @export
predict.lmvt <- function(object, newdata, threshold = NULL,
                         innov_g = TRUE, innov_t = FALSE, df_t = 6, ...) {

    stopifnot(inherits(object, "lmvt"))
    stopifnot(all(c("s", "t", "y") %in% names(newdata)))
    newdata <- newdata[order(newdata$s, newdata$t), , drop = FALSE]
    ids <- object$ids

    if (is.null(threshold))
        threshold <- stats::quantile(newdata$y, 0.9, na.rm = TRUE)

    # guard: unseen subjects
    unseen <- setdiff(unique(newdata$s), ids)
    if (length(unseen))
        stop(sprintf("predict.lmvt: newdata contains unseen subject IDs: %s",
                     paste(unseen, collapse = ", ")))

    thr_per_row <- .resolve_threshold(threshold, newdata$s, ids)

    # align covariates exactly like training
    Xraw <- matrix(0, nrow = nrow(newdata), ncol = length(object$xcols))
    colnames(Xraw) <- object$xcols
    present <- intersect(object$xcols, names(newdata))
    if (length(present))
        Xraw[, match(present, object$xcols)] <-
            data.matrix(newdata[, present, drop = FALSE])
    X <- if (isTRUE(object$standardize_X))
        scale(Xraw, center = object$x_center, scale = object$x_scale)
    else
        Xraw
    attr(X, "scaled:center") <- NULL
    attr(X, "scaled:scale")  <- NULL

    Y <- as.numeric(newdata$y)

    lag_by_id <- function(v, L) {
        out <- rep(NA_real_, length(v))
        if (L <= 0L) return(v)
        for (id in ids) {
            idx <- which(newdata$s == id)
            if (length(idx) > L) {
                out[idx[(L + 1L):length(idx)]] <-
                    v[idx[1:(length(idx) - L)]]
            }
        }
        out
    }

    p <- object$orders$p; q <- object$orders$q
    r <- object$orders$r; s_ord <- object$orders$s

    # mean-side lags
    Ylags <- if (p > 0L)
        do.call(cbind, lapply(seq_len(p), function(L) lag_by_id(Y, L)))
    else NULL
    Xlags_list_mean <- lapply(0:q, function(L)
        if (L == 0L) X else apply(X, 2, lag_by_id, L = L))
    Xmean <- do.call(cbind, Xlags_list_mean)

    # variance-side design follows training choice
    Xvar <- if (isTRUE(object$use_x_in_variance))
        Xmean else matrix(0, nrow(newdata), 0)

    # conditional mean & residuals
    AR_part <- if (p > 0L)
        rowSums(Ylags %*% matrix(object$theta, ncol = 1), na.rm = TRUE) else 0
    X_part  <- as.numeric(Xmean %*% object$beta)
    muhat   <- object$alpha[match(newdata$s, ids)] + AR_part + X_part
    eps_hat <- Y - muhat
    eps_hat[!is.finite(eps_hat)] <- 0

    # forward variance recursion with feasibility repair
    u_floor <- 1e-8
    sigma2 <- rep(NA_real_, nrow(newdata))
    for (id in ids) {
        idx <- which(newdata$s == id)
        for (ii in seq_along(idx)) {
            tt <- idx[ii]

            # ARCH/GARCH parts using available lags within subject
            arch_part <- 0; garch_part <- 0
            if (r > 0L) {
                for (L in seq_len(r)) {
                    tprev <- which(newdata$t[idx] == (newdata$t[tt] - L))
                    if (length(tprev) == 1L)
                        arch_part <- arch_part +
                            object$a[L] * (eps_hat[idx[tprev]]^2)
                }
            }
            if (s_ord > 0L) {
                for (L in seq_len(s_ord)) {
                    tprev <- which(newdata$t[idx] == (newdata$t[tt] - L))
                    if (length(tprev) == 1L)
                        garch_part <- garch_part +
                            object$b[L] * sigma2[idx[tprev]]
                }
            }

            # variance linear predictor
            u0 <- object$omega[match(id, ids)] + arch_part + garch_part
            xrow_var <- if (ncol(Xvar)) Xvar[tt, ] else numeric(0)
            xg <- if (length(object$gamma))
                sum(object$gamma * xrow_var) else 0
            vlin <- u0 + xg

            # Feasibility repair: prediction-time analogue of training B2 guards
            if (!is.finite(vlin) || vlin <= u_floor) {
                if (is.finite(xg) && xg < 0) {
                    s <- (u_floor - u0) / xg  # xg < 0
                    s <- min(1, max(0, s))
                    vlin_try <- u0 + s * xg
                    if (is.finite(vlin_try) && vlin_try > u_floor) {
                        vlin <- vlin_try
                    } else {
                        vlin <- max(u0, u_floor)
                    }
                } else {
                    vlin <- max(u0, u_floor)
                }
            }

            sigma2[tt] <- vlin
        }
    }

    sigmahat <- sqrt(pmax(sigma2, u_floor))
    z <- (thr_per_row - muhat) / sigmahat

    out <- data.frame(s = newdata$s, t = newdata$t,
                      yhat = muhat, muhat = muhat,
                      sigmahat = sigmahat,
                      threshold_used = thr_per_row)
    if (innov_g) out$risk_g <- 1 - stats::pnorm(z)
    if (innov_t) {
        if (!is.numeric(df_t) || length(df_t) != 1L || !is.finite(df_t) || df_t <= 2)
            stop("df_t must be a single number greater than 2 for standardized Student-t innovations.")
        out$risk_t <- 1 - stats::pt(z * sqrt(df_t / (df_t - 2)), df = df_t)
    }
    out
}

# ---- threshold resolution helper --------------------------------------------

.resolve_threshold <- function(threshold, new_s, train_ids) {
    n <- length(new_s)
    # 1) single fixed threshold
    if (length(threshold) == 1L) {
        thr <- rep(as.numeric(threshold), n)
        if (!all(is.finite(thr))) stop("Non-finite threshold.")
        return(thr)
    }
    # 2) per-row vector (length == nrow(newdata))
    if (length(threshold) == n) {
        thr <- as.numeric(threshold)
        if (any(!is.finite(thr)))
            stop("Non-finite values in threshold vector.")
        return(thr)
    }
    # 3) per-subject vector (preferred)
    if (!is.null(names(threshold))) {
        thr_map <- stats::setNames(as.numeric(threshold),
                                   as.character(names(threshold)))
        thr <- thr_map[as.character(new_s)]
        if (any(is.na(thr))) {
            missing_ids <- unique(new_s[is.na(thr)])
            stop(sprintf("Thresholds missing for subjects: %s",
                         paste(missing_ids, collapse = ", ")))
        }
        if (any(!is.finite(thr)))
            stop("Non-finite values in named threshold vector.")
        return(thr)
    }
    # 4) no names: try matching by subject set size
    new_ids_unique <- unique(new_s)
    if (length(threshold) == length(new_ids_unique)) {
        thr_map <- stats::setNames(as.numeric(threshold),
                                   as.character(new_ids_unique))
        thr <- thr_map[as.character(new_s)]
        if (any(!is.finite(thr)))
            stop("Non-finite values in per-subject threshold vector.")
        return(thr)
    }
    if (length(threshold) == length(train_ids)) {
        thr_map <- stats::setNames(as.numeric(threshold),
                                   as.character(train_ids))
        thr <- thr_map[as.character(new_s)]
        if (any(is.na(thr))) {
            missing_ids <- unique(new_s[is.na(thr)])
            stop(sprintf("Thresholds (aligned to training IDs) missing for subjects: %s",
                         paste(missing_ids, collapse = ", ")))
        }
        if (any(!is.finite(thr)))
            stop("Non-finite values in per-subject threshold vector.")
        return(thr)
    }
    stop("`threshold` must be: length 1; length nrow(newdata); ",
         "or a per-subject vector whose length equals the number of ",
         "subjects (in newdata or in training). Provide names(threshold) = ",
         "subject IDs to be explicit.")
}
