test_that("read_icp() yields a tibble with the expected number of rows", {
  file <- list.files(
    path = system.file("extdata", package = "cwrshelpr"),
    full.names = TRUE
  )
  expect_equal(nrow(read_icp(file[1])), 1440)
})
