test_that("summary.lmvt returns model diagnostics and coefficient table", {
    set.seed(4)
    sim <- simulate_scenario(scen = 1, S = 4, T = 60,
                             d_noise = 3, noise_kind = "iid", seed = 4)
    df  <- sim[, c("s", "t", "y",
                   grep("^X|^Noise|^Xbin", names(sim), value = TRUE))]
    fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
                lambda_beta = 0.05, lambda_gamma = 0.05,
                maxit = 8L)

    sm <- summary(fit)
    expect_s3_class(sm, "summary.lmvt")
    expect_s3_class(sm$table, "data.frame")
    expect_s3_class(sm$model, "data.frame")
    expect_s3_class(sm$nonzero, "data.frame")
    expect_true(all(c("equation", "block", "parameter", "term", "estimate",
                      "abs_estimate", "selected") %in% names(sm$table)))
    expect_true(any(sm$table$parameter == "beta"))
    expect_true(any(sm$table$parameter == "gamma"))

    sm_nz <- summary(fit, include_zero = FALSE, sort_by = "abs_estimate")
    expect_true(all(sm_nz$table$selected))

    out <- capture.output(print(sm, n = 5))
    expect_true(any(grepl("Summary of varGuidTS::lmvt fit", out)))
    expect_true(any(grepl("Coefficient summary table", out)))
})
