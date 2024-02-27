test_that("argument error prevention works", {
  expect_error(
    tpSVG(SpatialExperiment()),
    "Can't find assay in spe"
  )
  expect_error(
    tpSVG(SpatialExperiment()))
  expect_error(tpSVG(data.frame(sizeFactor = c(1))),
               "Can't find assay in spe")
  expect_error(tpSVG(data.frame(sizeFactor = c(1,2,3,4)), X = c(1)))
})

test_that("farmily error prevention works", {
  expect_error(
    tpSVG(data.frame(sizeFactor = c(1,2,3,4)),
          family = "poisson"))
  expect_error(tpSVG(data.frame(), family = "tmp"))
  expect_error(tpSVG(data.frame(), family = "negbin"))
})
