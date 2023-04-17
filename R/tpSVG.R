#' Title
#'
#' @param input
#' @param spatial_coords
#' @param X
#' @param assay_name
#' @param n_threads
#' @param BPPARAM
#' @param verbose
#'
#' @return
#'
#' @importFrom mgcv gam anova.gam
#' @importFrom stats var anova
#' @importFrom methods is
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment counts logcounts
#' @importFrom SummarizedExperiment assayNames assays rowData 'rowData<-'
#' @importFrom BiocParallel bplapply MulticoreParam
#'
#' @export
#'
#'
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
#' library(scran)
#' library(nnSVG)
#'
#' # load example dataset from STexampleData package
#' spe <- Visium_humanDLPFC()
#'
#' # preprocessing steps
#'
#' # keep only spots over tissue
#' spe <- spe[, colData(spe)$in_tissue == 1]
#'
#' # skip spot-level quality control, since this has been performed previously
#' # on this dataset
#'
#' # filter low-expressed and mitochondrial genes
#' spe <- filter_genes(spe)
#'
#' # calculate logcounts (log-transformed normalized counts) using scran package
#' # using library size factors
#' spe <- computeLibraryFactors(spe)
#' spe <- logNormCounts(spe)
#'
#' # select small number of genes for faster runtime in this example
#' set.seed(123)
#' ix <- sample(seq_len(nrow(spe)), 4)
#' spe <- spe[ix, ]
#'
#' # run nnSVG
#' set.seed(123)
#' spe <- tpSVG(spe)
#'
#' # show results
#' # for more details see extended example in vignette
#' rowData(spe)
tpSVG <- function(input, spatial_coords = NULL, X = NULL,
                  assay_name = "logcounts",
                  n_threads = 1, BPPARAM = NULL,
                  verbose = FALSE) {


  # NOTE: Some code blocks are borrowed from nnSVG by Lukas M Weber.
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
  }

  if (!is.null(X)) {
    stop("Not implemented yet")
    stopifnot(nrow(X) == ncol(input))
  }

  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }


  if (is(input, "SpatialExperiment")) {
    y <- assays(spe)[[assay_name]]
    coords <- spatialCoords(spe)
  } else {
    y <- input
    coords <- spatial_coords
    row_names <- rownames(input)
  }

  # scale coordinates proportionally
  # range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  # coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)




  stopifnot(ncol(coords)==2)


  colnames(coords) <- c("coord_x", "coord_y")

  fit.df <- data.frame(coords)

  ix <- seq_len(nrow(y))
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

  # --------------------
  # calculate statistics
  # --------------------

  # if (is(input, "SpatialExperiment") & ("logcounts" %in% assayNames(spe))) {
  #   lc <- logcounts(spe)
  #   # mean logcounts
  #   mat_tp <- cbind(
  #     mat_tp,
  #     mean = rowMeans(lc)
  #   )
  #   # variance of logcounts
  #   mat_tp <- cbind(
  #     mat_tp,
  #     var = rowVars(as.matrix(lc))
  #   )
  #   # spatial coefficient of variation
  #   mat_tp <- cbind(
  #     mat_tp,
  #     spcov = sqrt(mat_tp[, "sigma.sq"]) / mat_tp[, "mean"]
  #   )
  # } else {
  #   # return NAs if logcounts not provided
  #   mat_tp <- cbind(
  #     mat_tp, mean = NA, var = NA, spcov = NA
  #   )
  # }

  # proportion of spatial variance out of total variance
  # mat_tp <- cbind(
  #   mat_tp,
  #   prop_sv = mat_tp[, "sigma.sq"] / (mat_tp[, "sigma.sq"] + mat_tp[, "tau.sq"])
  # )

  # ------------------------------------------
  # likelihood ratio (LR) statistics and tests
  # ------------------------------------------

  # if (is(input, "SpatialExperiment")) {
  #   nrows <- nrow(spe)
  #   ncols <- ncol(spe)
  # } else {
  #   nrows <- nrow(input)
  #   ncols <- ncol(input)
  # }

  # calculate log likelihoods for nonspatial models

  # loglik_lm <- vapply(seq_len(nrows), function(i) {
  #   y_i <- y[i, ]
  #   if (is.null(X)) {
  #     X <- rep(1, ncols)
  #   }
  #   as.numeric(logLik(lm(y_i ~ X)))
  # }, numeric(1))
  #
  # mat_tp <- cbind(
  #   mat_tp,
  #   loglik_lm = loglik_lm
  # )

  # calculate LR statistics and tests (Wilks' theorem, asymptotic chi-square
  # with 2 degrees of freedom)

  # LR_stat <- -2 * (mat_tp[, "loglik_lm"] - mat_tp[, "loglik"])

  # pval <- 1 - pchisq(LR_stat, df = 2)
  # padj <- p.adjust(pval, method = "BH")

  # rank SVGs according to LR statistics
  # LR_rank <- rank(-1 * LR_stat)

  # mat_tp <- cbind(
  #   mat_tp,
  #   LR_stat = LR_stat,
  #   rank = LR_rank,
  #   pval = pval,
  #   padj = padj
  # )

  # --------------
  # return outputs
  # --------------

  if (is(input, "SpatialExperiment")) {
    # return in rowData of spe object
    stopifnot(nrow(spe) == nrow(mat_tp))
    rowData(spe) <- cbind(rowData(spe), mat_tp)
    spe
  } else {
    # return as numeric matrix
    stopifnot(nrow(input) == nrow(mat_tp))
    rownames(mat_tp) <- row_names
    mat_tp
  }
}
