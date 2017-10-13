context("diffExpressedVariants_count")
test_that("diffExpressedVariants function works as expected on count file", {
  # test on count file
  fpath1 <- system.file("extdata", "table_counts_alt_splicing.txt", package="kissDE")
  tableCounts <- read.table(fpath1, head = TRUE)
  conditions <- c("C1", "C1", "C2", "C2")
  diff <- diffExpressedVariants(tableCounts, conditions)
  expect_equal(names(diff), c("finalTable", "correctedPVal", "uncorrectedPVal", "resultFitNBglmModel", "f/psiTable", "k2rgFile"))
  expect_equal(dim(diff$finalTable)[2], 13)
  expect_equal(dim(diff$finalTable)[1], dim(tableCounts)[1]/2)
  expect_equal(dim(diff$finalTable)[1], length(diff$correctedPVal))
  expect_equal(dim(diff$finalTable)[1], dim(diff$`f/psiTable`)[1])
  expect_equal(names(diff$correctedPVal), names(diff$uncorrectedPVal))
  expect_equal(names(diff$correctedPVal), rownames(diff$resultFitNBglmModel))
  expect_null(diff$k2rgFile)
  # TODO : tests the final results -> number of significants events...
})
  