test_that("Argument error prevention works", {

  # Create some mock data for testing
  input <- matrix(runif(100), nrow = 10)
  row_names <- paste0("Row", 1:10)
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = input))
  empty_spe <- SpatialExperiment::SpatialExperiment()


  expect_error(
    tpSVG(empty_spe),
    "Can't find assay in spe"
  )

  expect_error(tpSVG(data.frame(sizeFactor = c(1))))

  expect_error(tpSVG(data.frame(sizeFactor = c(1,2,3,4)), X = c(1)))
})

test_that("Farmily error prevention works", {
  expect_error(
    tpSVG(data.frame(sizeFactor = c(1,2,3,4)),
          family = "poisson"))

  expect_error(tpSVG(data.frame(), family = "tmp"))
})
