context("diffExpressedVariants_KS")
test_that("diffExpressedVariants function works as expected on kissplice fasta file", {
  # test on kissplice fasta file
  fpath1 <- system.file("extdata", "output_kissplice_SNV.fa", package = "kissDE")
  mySNVcounts <- kissplice2counts(fpath1, pairedEnd = TRUE)
  mySNVconditions <- c("C1", "C1", "C2", "C2")
  diffSNV <- diffExpressedVariants(mySNVcounts, mySNVconditions)
  expect_equal(names(diffSNV), c("finalTable", "correctedPVal", "uncorrectedPVal", "resultFitNBglmModel", "f/psiTable", "k2rgFile"))
  expect_equal(dim(diffSNV$finalTable)[2], 13)
  expect_equal(dim(diffSNV$finalTable)[1], length(diffSNV$correctedPVal))
  expect_equal(dim(diffSNV$finalTable)[1], dim(diffSNV$`f/psiTable`)[1])
  expect_equal(names(diffSNV$correctedPVal), names(diffSNV$uncorrectedPVal))
  expect_equal(names(diffSNV$correctedPVal), rownames(diffSNV$resultFitNBglmModel))
  expect_null(diffSNV$k2rgFile)
  # TODO : tests the final results -> number of significants events...
})