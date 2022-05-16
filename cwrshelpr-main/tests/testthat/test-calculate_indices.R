test_that("calculate_indices() returns the expected values", {
  file <- list.files(
    path = system.file("extdata", package = "cwrshelpr"),
    full.names = TRUE,
    pattern = ".+\\.csv"
  )
  data <- read_feem(file[2])
  ix <- calculate_indices(data)
  expect_equal(
    ix$value,
    c(0.5081594, 13.6475984)
  )
})
