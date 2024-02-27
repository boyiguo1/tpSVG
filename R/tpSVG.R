#' Thin Plate Spline Model to Detect Spatially Variable Genes
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
#' @param family a description of the error distribution and link function
#'   to be used in the model. Currently support two distributions \code{poisson}
#'   and \code{gaussian}
#' @param offset This can be used to account for technician variation when
#'   \code{family = poisson} model is used to model raw counts. \code{offset}
#'   should take in the log-transformed scale factor, e.g.
#'   \code{offset = log(spe$sizeFactor)}, library size, or other normalization
#'   factor.
#' @param weights Reserved for future development, e.g. correcting mean-var
#'   relationship for Gaussian models. Please use with caution.
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
#'   reporting the thinplate spline model. The \code{test_stat} entry of the
#'   returned object is the test statistic for the corresponding model,
#'    that is F statistics for the gaussian model and the Chi-squared statistics
#'    for generalized models.
#'
#' @importFrom mgcv gam anova.gam  negbin
#' @importFrom stats var anova gaussian poisson as.formula family
#' @importFrom methods is
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment counts logcounts
#' @importFrom SummarizedExperiment assayNames assays rowData 'rowData<-'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom MatrixGenerics colSums2
#'
#' @export
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
#' # Add library size
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
#' # run tpSVG
#' set.seed(123)
#'
#' # Gaussian Model
#' spe_gaus <- tpSVG(
#'  spe,
#'  family = gaussian(),
#'  assay_name = "logcounts"
#'  )
#'
#' # Poisson Model
#' spe_poisson  <- tpSVG(
#'  spe,
#'  family = poisson,
#'  assay_name = "counts",
#'  offset = log(spe$sizeFactor)   # Natural log library size
#'  )

tpSVG <- function(
    input,
    spatial_coords = NULL,
    X = NULL,
    family = poisson(),
    offset = log(input$sizeFactor),
    weights = NULL,
    assay_name = "counts",
    n_threads = 1,
    BPPARAM = NULL,
    verbose = FALSE,
    ...
) {


  # NOTE: Some code blocks are borrowed from nnSVG by Lukas M Weber.
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot("Can't find assay in spe" = assay_name %in% assayNames(spe))
  }

  if (!is.null(X)) {
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


  # Check family
  # NOTE: the code is copied from "stats::glm"
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    stop("'family' not recognized")
  }
  # Copy end
  if (!(family$family %in% c("gaussian", "negative binomial", "poisson"))){
    stop(
      "'family' has to be one of the following distributions: gaussian,",
      "negative binomial (negbin), poisson"
    )
  }
  # New function ends



  flag_count_mdl <- family$family %in% c("negative binomial", "poisson")

  if (is.null(offset))
    offset <- rep(0, ncol(y))
  if (is.null(weights))
    weights <- rep(1, ncol(y))

  if(flag_count_mdl){
    if (is.null(offset)){
      warning("Using count-based model without providing offset.",
              "Library size is calculated with column sums of the count matrix")
      offset <- log(colSums2(y))
    }
  }



  stopifnot(ncol(coords)==2)
  colnames(coords) <- c("coord_x", "coord_y")


  if(!is.null(ncol(X)))
    stop("not implemented for complex desgin matrix.")
  if(is.null(ncol(X)))
    fit.df <- data.frame( coords)
  else {
    fit.df <- data.frame(
      coords,
      X
    )
  }


  ix <- seq_len(nrow(y))
  out_tp <- bplapply(ix, function(i) {
    # fit model (intercept-only model if x is NULL)
    fit.df$tmp_y <- y[i, ]
    runtime <- system.time({
      tp_formula <- "tmp_y~s(coord_x, coord_y, bs='tp')"
      if(!is.null(X)){
        tp_formula <- paste0(tp_formula, " + ",
                             paste0( "X", collapse = "+")
        )
      }

      tp_mdl <- gam(
        formula = tp_formula|> as.formula() ,
        family = family, data = fit.df,
        offset = offset, weights = weights)
    })

    if(flag_count_mdl) {
      res_i <- c(
        test_stat = anova(tp_mdl)$s.table[,"Chi.sq"],
        raw_p = anova(tp_mdl)$s.table[,"p-value"],
        runtime = runtime[["elapsed"]]
      )
    } else{
      res_i <- c(
        test_stat = anova(tp_mdl)$s.table[,"F"],
        raw_p = anova(tp_mdl)$s.table[,"p-value"],
        runtime = runtime[["elapsed"]]
      )
    }



    res_i
  }, BPPARAM = BPPARAM)

  # collapse output list into matrix
  mat_tp <- do.call("rbind", out_tp)

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
