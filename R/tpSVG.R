#' Title
#'
#' @param input \code{SpatialExperiment} or \code{numeric} matrix: Input data,
#'   which can either be a \code{SpatialExperiment} object or a \code{numeric}
#'   matrix of values. If it is a \code{SpatialExperiment} object, it is assumed
#'   to have an \code{assay} slot containing either logcounts (e.g. from the
#'   \code{scran} package) or deviance residuals (e.g. from the \code{scry}
#'   package), and a \code{spatialCoords} slot containing spatial coordinates of
#'   the measurements. If it is a \code{numeric} matrix, the values are assumed
#'   to already be normalized and transformed (e.g. logcounts), formatted as
#'   \code{rows = genes} and \code{columns = spots}, and a separate
#'   \code{numeric} matrix of spatial coordinates must also be provided with the
#'   \code{spatial_coords} argument.
#' @param spatial_coords \code{numeric} matrix: Matrix containing columns of
#'   spatial coordinates, formatted as \code{rows = spots}. This must be
#'   provided if \code{input} is provied as a \code{numeric} matrix of values,
#'   and is ignored if \code{input} is provided as a \code{SpatialExperiment}
#'   object. Default = NULL.
#' @param X \code{numeric} matrix: Optional design matrix containing columns of
#'   covariates per spatial location, e.g. known spatial domains. Number of rows
#'   must match the number of spatial locations. Default = NULL, which fits an
#'   intercept-only model.
#' @param assay_name \code{character}: If \code{input} is provided as a
#'   \code{SpatialExperiment} object, this argument selects the name of the
#'   \code{assay} slot in the input object containing the preprocessed gene
#'   expression values. For example, \code{logcounts} for log-transformed
#'   normalized counts from the \code{scran} package, or
#'   \code{binomial_deviance_residuals} for deviance residuals from the
#'   \code{scry} package. Default = \code{"logcounts"}, or ignored if
#'   \code{input} is provided as a \code{numeric} matrix of values.
#' @param family TODO: fixme
#' @param n_threads \code{integer}: Number of threads for parallelization.
#'   Default = 1. We recommend setting this equal to the number of cores
#'   available (if working on a laptop or desktop) or around 10 or more (if
#'   working on a compute cluster).
#' @param BPPARAM \code{BiocParallelParam}: Optional additional argument for
#'   parallelization. This argument is provided for advanced users of
#'   \code{BiocParallel} for further flexibility for parallelization on some
#'   operating systems. If provided, this should be an instance of
#'   \code{BiocParallelParam}. For most users, the recommended option is to use
#'   the \code{n_threads} argument instead. Default = NULL, in which case
#'   \code{n_threads} will be used instead.
#' @param verbose \code{logical}: Whether to display verbose output for model
#'   fitting and parameter estimation from \code{BRISC}. Default = FALSE.
#' @param ... Reserved for future arguments.
#'
#' @return If the input was provided as a \code{SpatialExperiment} object, the
#'   output values are returned as additional columns in the \code{rowData} slot
#'   of the input object. If the input was provided as a \code{numeric} matrix
#'   of values, the output is returned as a \code{numeric} matrix. The output
#'   values include p-values without any adjustment and statistics reporting
#'   reporting the thinplate spline model.
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
#' Add library size
#' spe <- addPerCellQCMetrics(spe)
#'
#' # filter low-expressed and mitochondrial genes
#' spe <- filter_genes(spe)
#'
#' # calculate logcounts (log-transformed normalized counts) using scran package
#' # using library size factors
#' spe <- computeLibraryFactors(spe)
#' spe <- logNormCounts(spe)
#'
#'
#'
#' # select small number of genes for faster runtime in this example
#' set.seed(123)
#' ix <- sample(seq_len(nrow(spe)), 4)
#' spe <- spe[ix, ]
#'
#' # run nnSVG
#' set.seed(123)
#' spe_gaus <- tpSVG(spe)
#' spe_poisson  <- tpSVG(spe, family = poisson,
#'  assay_name = "counts", offset = log(spe$total))
#'
#' # show results
#' # for more details see extended example in vignette
#' rowData(spe_gaus)
tpSVG <- function(input, spatial_coords = NULL, X = NULL,
                  family = gaussian(),
                  offset = NULL,
                  assay_name = "logcounts",
                  n_threads = 1, BPPARAM = NULL,
                  verbose = FALSE, ...) {


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


  # TODO: create isolated internal function
  # Check family
  # NOTE: the code is copied from "stats::glm"
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  # Copy end
  if (!(family$family %in% c("gaussian", "negative binomial", "poisson"))){
    print(family)
    stop("'family' has to be one of the following distributions: gaussian,",
         "negative binomial (negbin), poisson")
  }
  # New function ends



  flag_count_mdl <- family$family %in% c("negative binomial", "poisson")

  if (is.null(offset))
    offset <- rep(0, ncol(y))

  browser()
  if(flag_count_mdl){
    if (is.null(offset)){
      warning("Using count-based model without supplying offset.",
              "Library size is calculated with column sums of the count matrix")
      offset <- log(colSums2(y))
    }
    # TODO: Check if y matrix is count matrix
  }



  stopifnot(ncol(coords)==2)
  colnames(coords) <- c("coord_x", "coord_y")

  fit.df <- data.frame(coords)

  ix <- seq_len(nrow(y))
  out_tp <- bplapply(ix, function(i) {
    # fit model (intercept-only model if x is NULL)
    fit.df$tmp_y <- y[i, ]
    runtime <- system.time({
      tp_mdl <- gam(
        tmp_y~s(coord_x, coord_y, bs="tp"),
        family = family, data = fit.df,
        offset = offset, weights = weights)
    })

    browser()

    if(flag_count_mdl) {
      res_i <- c(
        # TODO: these implementation wont work for situation where covariates are allowed
        # F_stat = anova(out_i)$s.table[,"F"],
        # loglik = out_i$log_likelihood,
        Chi.sq_stat = anova(tp_mdl)$s.table[,"Chi.sq"],
        raw.p = anova(tp_mdl)$s.table[,"p-value"],
        # F_stat = anova(tp_mdl)$s.table[,"F"],
        # raw_p = anova(tp_mdl)$s.table[,"p-value"],
        # GE_mean = tp_mdl$coefficients[1] |> unname(),
        # tp_edf = tp_mdl$edf |> sum() - 1,
        # tp_ref.df = tp_mdl$edf1 |> sum() - 1,
        # residual_var = tp_mdl$residuals |> var(),
        runtime = runtime[["elapsed"]]
      )
    } else{
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
    }



    res_i
  }, BPPARAM = BPPARAM)

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
