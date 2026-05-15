# Extracted R code from intro.Rmd. Code chunks are marked eval=FALSE in the vignette.

library(varGuidTS)

set.seed(1)
sim <- simulate_scenario(scen = 2, S = 4, T = 60,
                         d_noise = 4, noise_kind = "ar1", seed = 1)

df <- sim[, c("s", "t", "y",
              grep("^X|^Noise|^Xbin", names(sim), value = TRUE))]

fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
            lambda_beta = 0.05, lambda_gamma = 0.05,
            use_x_in_variance = TRUE, maxit = 8)

fit

sm <- summary(fit, include_zero = FALSE, sort_by = "abs_estimate")
sm$model
sm$nonzero
head(sm$table)

# Pooled threshold
pred_global <- predict(fit, df, threshold = quantile(df$y, 0.90))
head(pred_global)

# Subject-specific thresholds
subject_cutoffs <- tapply(df$y, df$s, quantile, probs = 0.75, na.rm = TRUE)

pred_subject <- predict(fit, df,
                        threshold = subject_cutoffs,
                        innov_g = TRUE,
                        innov_t = TRUE,
                        df_t = 8)
head(pred_subject)

s1 <- simulate_scenario(scen = 1, S = 3, T = 50, d_noise = 3, seed = 10)
s2 <- simulate_scenario(scen = 2, S = 3, T = 50, d_noise = 3, seed = 20)
s3 <- simulate_scenario(scen = 3, S = 3, T = 50, d_noise = 3, seed = 30)
s4 <- simulate_scenario(scen = 4, S = 3, T = 50, d_noise = 3, seed = 40)
