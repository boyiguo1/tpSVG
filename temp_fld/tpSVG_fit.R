tpSVG_fit <- function(
    x, y,
    weights = rep.int(1, nobs),    # Reserved for weighted analysis
    offset = rep.int(0, nobs), family = gaussian(),
    # Parallezation parameters
    n_threads = 1, BPPARAM = NULL,
    verbose = FALSE
){
  out_tp <- bplapply(ix, function(i) {
    # fit model (intercept-only model if x is NULL)
    fit.df$tmp_y <- y[i, ]
    runtime <- system.time({
      tp_mdl <- gam(tmp_y~s(coord_x, coord_y, bs="tp"), data = fit.df)
    })

    res_i <- c(
      # TODO: these implementation wont work for situation where covariates are allowed
      # F_stat = anova(out_i)$s.table[,"F"],
      # loglik = out_i$log_likelihood,
      F_stat = anova(tp_mdl)$s.table[,"F"],
      raw_p = anova(tp_mdl)$s.table[,"p-value"],
      GE_mean = tp_mdl$coefficients[1] |> unname(),
      tp_edf = tp_mdl$edf |> sum() - 1,
      tp_ref.df = tp_mdl$edf1 |> sum() - 1,
      residual_var = tp_mdl$residuals |> var(),
      runtime = runtime[["elapsed"]]
    )
    res_i
  }, BPPARAM = BPPARAM)

  browser()
  # collapse output list into matrix
  mat_tp <- do.call("rbind", out_tp)
}
