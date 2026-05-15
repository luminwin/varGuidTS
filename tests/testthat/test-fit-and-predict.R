test_that("simulate_scenario returns expected columns and shapes", {
    set.seed(1)
    sim <- simulate_scenario(scen = 1, S = 4, T = 60,
                             d_noise = 3, noise_kind = "iid", seed = 1)
    expect_s3_class(sim, "data.frame")
    expect_true(all(c("s","t","y","mu","sigma","pi_true",
                      "X1","X2","X3") %in% names(sim)))
    expect_equal(nrow(sim), 4 * 60)
    expect_true(all(is.finite(sim$y)))
})

test_that("lmvt returns a fitted object that print() can format", {
    set.seed(2)
    sim <- simulate_scenario(scen = 1, S = 4, T = 60,
                             d_noise = 3, noise_kind = "iid", seed = 2)
    df  <- sim[, c("s","t","y",
                   grep("^X|^Noise|^Xbin", names(sim), value = TRUE))]
    fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
                lambda_beta = 0.05, lambda_gamma = 0.05,
                maxit = 8L)
    expect_s3_class(fit, "lmvt")
    expect_named(fit, c("call","orders","penalties","ids",
                        "alpha","theta","beta","omega","a","b","gamma",
                        "xcols","standardize_X","use_x_in_variance",
                        "x_center","x_scale",
                        "converged","iters",
                        "mean_times","var_times","iter_times"),
                 ignore.order = TRUE)
    expect_length(fit$alpha, length(unique(df$s)))
    out <- capture.output(print(fit))
    expect_true(any(grepl("Penalized panel ARX-GARCHX fit", out)))
})

test_that("predict.lmvt returns finite muhat, sigmahat, and risk_g", {
    set.seed(3)
    sim <- simulate_scenario(scen = 1, S = 4, T = 60,
                             d_noise = 3, noise_kind = "iid", seed = 3)
    df  <- sim[, c("s","t","y",
                   grep("^X|^Noise|^Xbin", names(sim), value = TRUE))]
    fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
                lambda_beta = 0.05, lambda_gamma = 0.05,
                maxit = 8L)
    pr <- predict(fit, df, threshold = stats::quantile(df$y, 0.9))
    expect_true(all(c("muhat","sigmahat","risk_g","threshold_used")
                    %in% names(pr)))
    expect_true(all(is.finite(pr$muhat)))
    expect_true(all(pr$sigmahat > 0))
    expect_true(all(pr$risk_g >= 0 & pr$risk_g <= 1))
})

