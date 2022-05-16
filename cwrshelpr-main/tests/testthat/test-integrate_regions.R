test_that("integrate_regions() returns the expected values", {
  file <- list.files(
    path = system.file("extdata", package = "cwrshelpr"),
    full.names = TRUE,
    pattern = ".+\\.csv"
  )
  feem <- read_feem(file[1], truncate = TRUE)
  expect_equal(
    integrate_regions(feem)$integrated,
    c(46.64113, 246.98857, 2420.94329, 744.74940, 13312.19190, 21597.09472)
  )
})
