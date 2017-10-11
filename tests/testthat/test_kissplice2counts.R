context("kisslice2counts")
test_that("Loading data with kissplice2counts work as expected", {
  # test load kissplice fasta file
  fpath1 <- system.file("extdata", "output_kissplice_SNV.fa", package = "kissDE")
  mySNVcounts <- kissplice2counts(fpath1, pairedEnd = TRUE)
  expect_equal(names(mySNVcounts), c("countsEvents", "psiInfo", "exonicReadsInfo", "k2rgFile"))
  expect_null(mySNVcounts$k2rgFile)
  expect_true(mySNVcounts$exonicReadsInfo)
  expect_equal(mySNVcounts$countsEvents[, 1], mySNVcounts$psiInfo[, 1])
  expect_equal(dim(mySNVcounts$countsEvents)[2], 6)
  expect_equal(dim(mySNVcounts$countsEvents)[1], 126)
  expect_equal(dim(mySNVcounts$psiInfo)[2], 1)
  
  # test load kissplice2refgenome file
  fpath2 <- system.file("extdata", "output_k2rg_alt_splicing.txt", package = "kissDE")
  mySplicingcounts <- kissplice2counts(fpath2, pairedEnd = TRUE, k2rg = TRUE, counts = 2, exonicReads = FALSE)
  expect_equal(names(mySplicingcounts), c("countsEvents", "psiInfo", "exonicReadsInfo", "k2rgFile"))
  matches(mySplicingcounts$k2rgFile, "output_k2rg_alt_splicing.txt")
  expect_false(mySplicingcounts$exonicReadsInfo)
  expect_equal(mySplicingcounts$countsEvents[, 1], mySplicingcounts$psiInfo[, 1])
  expect_equal(dim(mySplicingcounts$countsEvents)[2], 6)
  expect_equal(dim(mySplicingcounts$countsEvents)[1], 212)
  expect_equal(dim(mySplicingcounts$psiInfo)[2], 5)
})
