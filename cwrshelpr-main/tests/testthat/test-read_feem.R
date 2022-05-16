test_that("read_feem() returns the expected number of rows", {
  file <- list.files(
    path = system.file("extdata", package = "cwrshelpr"),
    full.names = TRUE,
    pattern = ".+\\.csv"
  )
  expect_equal(
    dim(read_feem(file[1]))[1],
    15125
  )
})
