test_that("package loads and has no exports yet", {
  expect_true(requireNamespace("seqbench", quietly = TRUE))
  expect_length(getNamespaceExports("seqbench"), 0L)
})
