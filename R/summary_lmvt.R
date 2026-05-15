## ============================================================================
## summary.lmvt() : model and coefficient summary tables
## ============================================================================

# Build lagged-X term labels in the same order used by lmvt():
# cbind(X lag 0, X lag 1, ..., X lag q).
.lmvt_lagged_x_terms <- function(xcols, q, ncoef = NULL) {
    if (length(xcols) == 0L) {
        terms <- character(0)
        covars <- character(0)
        lags <- integer(0)
    } else {
        lags <- rep(0:q, each = length(xcols))
        covars <- rep(xcols, times = q + 1L)
        terms <- ifelse(lags == 0L, covars, paste0(covars, "_lag", lags))
    }

    if (!is.null(ncoef) && length(terms) != ncoef) {
        terms <- paste0("x", seq_len(ncoef))
        covars <- terms
        lags <- rep(NA_integer_, ncoef)
    }
    list(term = terms, covariate = covars, lag = lags)
}

.lmvt_coef_table <- function(object, tol = 1e-8) {
    stopifnot(inherits(object, "lmvt"))

    make_block <- function(equation, block, parameter, term, estimate,
                           subject = NA_character_, covariate = NA_character_,
                           lag = NA_integer_) {
        n <- length(estimate)
        if (n == 0L) return(NULL)
        data.frame(
            equation = rep(equation, n),
            block = rep(block, n),
            parameter = rep(parameter, n),
            term = as.character(term),
            subject = rep(subject, length.out = n),
            covariate = rep(covariate, length.out = n),
            lag = rep(lag, length.out = n),
            estimate = as.numeric(estimate),
            stringsAsFactors = FALSE
        )
    }

    ids_chr <- as.character(object$ids)
    q <- object$orders$q
    p <- object$orders$p
    r <- object$orders$r
    s_ord <- object$orders$s

    parts <- list()
    parts[[length(parts) + 1L]] <- make_block(
        equation = "mean", block = "subject_intercept", parameter = "alpha",
        term = paste0("alpha[", ids_chr, "]"), estimate = object$alpha,
        subject = ids_chr
    )

    if (length(object$theta)) {
        lag_theta <- seq_along(object$theta)
        parts[[length(parts) + 1L]] <- make_block(
            equation = "mean", block = "autoregressive", parameter = "theta",
            term = paste0("theta_lag", lag_theta), estimate = object$theta,
            lag = lag_theta
        )
    } else if (p > 0L) {
        parts[[length(parts) + 1L]] <- NULL
    }

    if (length(object$beta)) {
        xt <- .lmvt_lagged_x_terms(object$xcols, q, ncoef = length(object$beta))
        parts[[length(parts) + 1L]] <- make_block(
            equation = "mean", block = "covariate", parameter = "beta",
            term = paste0("beta[", xt$term, "]"), estimate = object$beta,
            covariate = xt$covariate, lag = xt$lag
        )
    }

    parts[[length(parts) + 1L]] <- make_block(
        equation = "variance", block = "subject_baseline", parameter = "omega",
        term = paste0("omega[", ids_chr, "]"), estimate = object$omega,
        subject = ids_chr
    )

    if (length(object$a)) {
        lag_a <- seq_along(object$a)
        parts[[length(parts) + 1L]] <- make_block(
            equation = "variance", block = "ARCH", parameter = "a",
            term = paste0("a_lag", lag_a), estimate = object$a,
            lag = lag_a
        )
    } else if (r > 0L) {
        parts[[length(parts) + 1L]] <- NULL
    }

    if (length(object$b)) {
        lag_b <- seq_along(object$b)
        parts[[length(parts) + 1L]] <- make_block(
            equation = "variance", block = "GARCH", parameter = "b",
            term = paste0("b_lag", lag_b), estimate = object$b,
            lag = lag_b
        )
    } else if (s_ord > 0L) {
        parts[[length(parts) + 1L]] <- NULL
    }

    if (length(object$gamma)) {
        xt <- .lmvt_lagged_x_terms(object$xcols, q, ncoef = length(object$gamma))
        parts[[length(parts) + 1L]] <- make_block(
            equation = "variance", block = "covariate", parameter = "gamma",
            term = paste0("gamma[", xt$term, "]"), estimate = object$gamma,
            covariate = xt$covariate, lag = xt$lag
        )
    }

    out <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
    if (is.null(out) || !nrow(out)) {
        out <- data.frame(equation = character(0), block = character(0),
                          parameter = character(0), term = character(0),
                          subject = character(0), covariate = character(0),
                          lag = integer(0), estimate = numeric(0),
                          stringsAsFactors = FALSE)
    }
    out$abs_estimate <- abs(out$estimate)
    out$selected <- out$abs_estimate > tol
    rownames(out) <- NULL
    out
}

#' Summarize a fitted \code{lmvt} object
#'
#' Produces a compact model summary and a coefficient table for all fitted
#' parameter blocks in a penalized panel ARX--GARCHX model. The coefficient
#' table is returned as \code{summary(fit)$table} and includes subject-specific
#' intercepts/baseline variances, autoregressive and GARCH parameters, and
#' mean/variance covariate coefficients. Because \code{lmvt()} fits a penalized
#' quasi-likelihood model, the table reports estimates and variable-selection
#' indicators rather than classical standard errors or p-values.
#'
#' @param object A fitted \code{lmvt} object.
#' @param tol Absolute value used to flag a coefficient as selected/nonzero.
#' @param include_zero Logical; if \code{FALSE}, only coefficients with
#'   \code{abs(estimate) > tol} are returned in the printed/returned table.
#' @param sort_by Sorting rule for the coefficient table. Use \code{"block"}
#'   for model order, \code{"abs_estimate"} for decreasing absolute estimate,
#'   or \code{"none"} to keep construction order.
#' @param top Optional integer limiting the returned coefficient table to the
#'   first \code{top} rows after filtering and sorting.
#' @param ... Unused, kept for S3 compatibility.
#'
#' @return An object of class \code{"summary.lmvt"}, a list containing
#'   \code{model} (model-level diagnostics), \code{nonzero} (nonzero counts by
#'   parameter block), and \code{table} (the coefficient summary table).
#'
#' @examples
#' \donttest{
#' set.seed(12)
#' sim <- simulate_scenario(scen = 1, S = 4, T = 60, d_noise = 4,
#'                          noise_kind = "iid", seed = 12)
#' df  <- sim[, c("s", "t", "y", grep("^X|^Noise|^Xbin",
#'                                    names(sim), value = TRUE))]
#' fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
#'             lambda_beta = 0.05, lambda_gamma = 0.05, maxit = 8)
#' sm <- summary(fit, include_zero = FALSE)
#' sm$table
#' }
#' @export
summary.lmvt <- function(object, tol = 1e-8, include_zero = TRUE,
                         sort_by = c("block", "abs_estimate", "none"),
                         top = NULL, ...) {
    stopifnot(inherits(object, "lmvt"))
    sort_by <- match.arg(sort_by)

    full_table <- .lmvt_coef_table(object, tol = tol)
    coef_table <- full_table
    if (!isTRUE(include_zero))
        coef_table <- coef_table[coef_table$selected, , drop = FALSE]

    if (sort_by == "abs_estimate") {
        coef_table <- coef_table[order(-coef_table$abs_estimate,
                                       coef_table$equation,
                                       coef_table$parameter,
                                       coef_table$term), , drop = FALSE]
    }
    # sort_by == "block" keeps the natural model-block construction order.

    if (!is.null(top)) {
        if (length(top) != 1L)
            stop("`top` must be a positive integer when supplied.")
        top <- as.integer(top)
        if (is.na(top) || top < 1L)
            stop("`top` must be a positive integer when supplied.")
        coef_table <- coef_table[seq_len(min(top, nrow(coef_table))), , drop = FALSE]
    }
    rownames(coef_table) <- NULL

    # Counts are based on the full coefficient table, not the filtered table.
    params <- c("alpha", "theta", "beta", "omega", "a", "b", "gamma")
    nonzero <- data.frame(
        parameter = params,
        n_total = vapply(params, function(z)
            as.integer(sum(full_table$parameter == z)), integer(1)),
        n_selected = vapply(params, function(z)
            as.integer(sum(full_table$parameter == z & full_table$selected)),
            integer(1)),
        stringsAsFactors = FALSE
    )
    nonzero <- nonzero[nonzero$n_total > 0L, , drop = FALSE]
    rownames(nonzero) <- NULL

    model <- data.frame(
        metric = c("subjects", "covariates", "mean_orders_p_q",
                   "variance_orders_r_s", "lambda_beta", "lambda_gamma",
                   "x_in_variance", "iterations", "converged",
                   "mean_time_per_iter", "variance_time_per_iter",
                   "total_time_per_iter"),
        value = c(
            as.character(length(object$ids)),
            as.character(length(object$xcols)),
            sprintf("(%d, %d)", object$orders$p, object$orders$q),
            sprintf("(%d, %d)", object$orders$r, object$orders$s),
            format(object$penalties$lambda_beta, digits = 4),
            format(object$penalties$lambda_gamma, digits = 4),
            as.character(isTRUE(object$use_x_in_variance)),
            as.character(object$iters),
            as.character(isTRUE(object$converged)),
            if (length(object$mean_times))
                format(mean(object$mean_times), digits = 4) else NA_character_,
            if (length(object$var_times))
                format(mean(object$var_times), digits = 4) else NA_character_,
            if (length(object$iter_times))
                format(mean(object$iter_times), digits = 4) else NA_character_
        ),
        stringsAsFactors = FALSE
    )

    out <- list(call = object$call, model = model, nonzero = nonzero,
                table = coef_table, full_table = full_table, tol = tol,
                include_zero = include_zero, sort_by = sort_by)
    class(out) <- "summary.lmvt"
    out
}

# Internal S3 print method; registered in NAMESPACE.
print.summary.lmvt <- function(x, digits = max(3L, getOption("digits") - 3L),
                               n = 30L, ...) {
    stopifnot(inherits(x, "summary.lmvt"))
    cat("Summary of varGuidTS::lmvt fit\n")
    if (!is.null(x$call)) {
        cat("Call: ")
        print(x$call)
    }

    cat("\nModel diagnostics:\n")
    print(x$model, row.names = FALSE)

    cat("\nNonzero coefficient counts (|estimate| > ", x$tol, "):\n", sep = "")
    print(x$nonzero, row.names = FALSE)

    if (length(n) != 1L) n <- 30L
    n <- as.integer(n)
    if (is.na(n) || n < 1L) n <- 30L
    cat("\nCoefficient summary table")
    if (nrow(x$table) > n)
        cat(sprintf(" (showing first %d of %d rows)", n, nrow(x$table)))
    cat(":\n")

    tab <- x$table[seq_len(min(n, nrow(x$table))), , drop = FALSE]
    if (nrow(tab)) {
        tab$estimate <- signif(tab$estimate, digits)
        tab$abs_estimate <- signif(tab$abs_estimate, digits)
    }
    print(tab, row.names = FALSE)
    invisible(x)
}
