context("kisslice2counts")
test_that("loading data with kissplice2counts work as expected", {
  ## test load kissplice fasta file ##
  fpath1 <- system.file("extdata", "output_kissplice_SNV.fa", package = "kissDE")
  mySNVcounts <- kissplice2counts(fpath1, counts = 0, pairedEnd = TRUE)
  expect_equal(names(mySNVcounts), c("countsEvents", "psiInfo", "exonicReadsInfo", "k2rgFile"))
  expect_null(mySNVcounts$k2rgFile)
  expect_true(mySNVcounts$exonicReadsInfo)
  expect_equal(mySNVcounts$countsEvents[, 1], mySNVcounts$psiInfo[, 1])
  expect_equal(dim(mySNVcounts$countsEvents)[2], 6)
  expect_equal(dim(mySNVcounts$countsEvents)[1], 126)
  expect_equal(dim(mySNVcounts$psiInfo)[2], 1)
  ## test computed value ##
  realCountsSNV <- data.frame(counts1 = c(910, 26), counts2 = c(1687, 22), counts3 = c(5, 28), counts4 = c(70, 8569))
  row.names(realCountsSNV) <- NULL
  computedCountsSNV <- mySNVcounts$countsEvents[which(mySNVcounts$countsEvents[,1] == "bcc_44787|Cycle_421687"),c(3:6)]
  row.names(computedCountsSNV) <- NULL
  expect_equal(computedCountsSNV, realCountsSNV)
  
  ## test load kissplice2refgenome file ##
  fpath2 <- system.file("extdata", "output_k2rg_alt_splicing.txt", package = "kissDE")
  ## exonicReads = FALSE ##
  mySplicingcounts <- kissplice2counts(fpath2, pairedEnd = TRUE, k2rg = TRUE, counts = 2, exonicReads = FALSE)
  expect_equal(names(mySplicingcounts), c("countsEvents", "psiInfo", "exonicReadsInfo", "k2rgFile"))
  expect_match(mySplicingcounts$k2rgFile, "output_k2rg_alt_splicing.txt")
  expect_false(mySplicingcounts$exonicReadsInfo)
  expect_equal(mySplicingcounts$countsEvents[, 1], mySplicingcounts$psiInfo[, 1])
  expect_equal(dim(mySplicingcounts$countsEvents)[2], 6)
  expect_equal(dim(mySplicingcounts$countsEvents)[1], 212)
  expect_equal(dim(mySplicingcounts$psiInfo)[2], 5)
  ## test computed value (with ASSB) ##
  realCountsSplicing_withASSB <- data.frame(counts1 = c(2, 33), counts2 = c(1, 14), counts3 = c(23, 6), counts4 = c(8, 3))
  row.names(realCountsSplicing_withASSB) <- NULL
  computedCountsSplicing_withASSB <- mySplicingcounts$countsEvents[which(mySplicingcounts$countsEvents[,1] == "bcc_68965|Cycle_4"),c(3:6)]
  row.names(computedCountsSplicing_withASSB) <- NULL
  expect_equal(computedCountsSplicing_withASSB, realCountsSplicing_withASSB)
  ## test computed value ##
  realCountsSplicing <- data.frame(counts1 = c(54, 49), counts2 = c(21, 23), counts3 = c(8, 19), counts4 = c(7, 41), row.names = c(25,26))
  row.names(realCountsSplicing) <- NULL
  computedCountsSplicing <- mySplicingcounts$countsEvents[which(mySplicingcounts$countsEvents[,1] == "bcc_140028|Cycle_1"),c(3:6)]
  row.names(computedCountsSplicing) <- NULL
  expect_equal(computedCountsSplicing, realCountsSplicing)
  ## exonicReads = TRUE ##
  mySplicingcounts2 <- kissplice2counts(fpath2, pairedEnd = TRUE, k2rg = TRUE, counts = 2, exonicReads = TRUE)
  ## test computed value (with exonic reads) ##
  realCountsSplicing2 <- data.frame(counts1 = c(150, 49), counts2 = c(98, 23), counts3 = c(35, 19), counts4 = c(31, 41))
  row.names(realCountsSplicing2) <- NULL
  computedCountsSplicing2 <- mySplicingcounts2$countsEvents[which(mySplicingcounts2$countsEvents[,1] == "bcc_140028|Cycle_1"),c(3:6)]
  row.names(computedCountsSplicing2) <- NULL
  expect_equal(computedCountsSplicing2, realCountsSplicing2)
})
