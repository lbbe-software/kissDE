context("writeOutputKissDE")
test_that("writeOutputKissDE work as expected", {
  # input = kissplice fasta file
  fpath1 <- system.file("extdata", "output_kissplice_SNV.fa", package = "kissDE")
  mySNVcounts <- kissplice2counts(fpath1, counts = 0, pairedEnd = TRUE)
  mySNVconditions <- c("C1", "C1", "C2", "C2")
  diffSNV <- diffExpressedVariants(mySNVcounts, mySNVconditions)
  writeOutputKissDE(diffSNV, output = "kissDE_output_SVN.tab")
  expect_true(file.exists("kissDE_output_SVN.tab"))
  file.remove("kissDE_output_SVN.tab")
  
  # input = kissplice2refgenome file
  fpath2 <- system.file("extdata", "output_k2rg_alt_splicing.txt", package = "kissDE")
  mySplicingcounts <- kissplice2counts(fpath2, pairedEnd = TRUE, k2rg = TRUE, counts = 2, exonicReads = FALSE)
  mySplicingconditions <- c("C1", "C1", "C2", "C2")
  diffSplicing <- diffExpressedVariants(mySplicingcounts, mySplicingconditions)
  writeOutputKissDE(diffSplicing, output = "kissDE_output_Splicing.tab")
  expect_true(file.exists("kissDE_output_Splicing.tab"))
  file.remove("kissDE_output_Splicing.tab")
})