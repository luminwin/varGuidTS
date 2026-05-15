# Internal utilities for varGuidTS.
# Not exported. Lightweight helpers used by lmvt() and predict.lmvt().

`%||%` <- function(x, y) if (is.null(x)) y else x

# Lag a numeric vector v within each subject, where ids and a sorted
# (s, t) data frame `df` are passed in. Returns NA for the first L
# observations of each subject.
.lag_by_id <- function(v, L, s_vec, ids) {
    out <- rep(NA_real_, length(v))
    if (L <= 0L) return(v)
    for (id in ids) {
        idx <- which(s_vec == id)
        if (length(idx) > L) {
            out[idx[(L + 1L):length(idx)]] <- v[idx[1:(length(idx) - L)]]
        }
    }
    out
}

# Soft-thresholding operator used in proximal updates.
.soft_threshold <- function(z, lambda) sign(z) * pmax(abs(z) - lambda, 0)

# Build the lagged-X design block: cbind(X, lag(X, 1), ..., lag(X, q)).
.lagged_X <- function(X, q, s_vec, ids) {
    if (q == 0L) return(X)
    blocks <- vector("list", q + 1L)
    blocks[[1L]] <- X
    for (L in seq_len(q)) {
        blocks[[L + 1L]] <- apply(X, 2L, .lag_by_id, L = L,
                                  s_vec = s_vec, ids = ids)
    }
    do.call(cbind, blocks)
}
