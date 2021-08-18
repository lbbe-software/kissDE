context("diffExpressedVariants_count")
test_that("diffExpressedVariants function works as expected on count file", {
  # test on count file
  fpath1 <- system.file("extdata/table_counts_alt_splicing.txt", package="kissDE", mustWork=TRUE)
  tableCounts <- read.table(fpath1, head = TRUE)
  conditions <- c("C1", "C1", "C2", "C2")
  diff <- diffExpressedVariants(tableCounts, conditions)
  expect_equal(names(diff), c("finalTable", "correctedPVal", "uncorrectedPVal", "resultFitNBglmModel", "f/psiTable", "k2rgFile"))
  expect_equal(dim(diff$finalTable)[2], 13)
  expect_equal(dim(diff$finalTable)[1], length(diff$correctedPVal))
  expect_equal(dim(diff$finalTable)[1], dim(diff$`f/psiTable`)[1])
  expect_equal(names(diff$correctedPVal), names(diff$uncorrectedPVal))
  expect_equal(names(diff$correctedPVal), rownames(diff$resultFitNBglmModel))
  expect_null(diff$k2rgFile)
  expect_equal(dim(diff$finalTable[which(diff$finalTable$Adjusted_pvalue <= 0.05),])[1], 19)
  expect_equal(diff$finalTable[which(diff$finalTable$ID == "event1"), "Deltaf/DeltaPSI"], -0.7824)
})

test_that("diffExpressedVariants function works as expected on count file with 2 cores", {
  # don't test on bioconductor
  skip_on_bioc()
  # test nbCore parameter
  if (detectCores() - 1 >= 2){
    fpath1 <- system.file("extdata/table_counts_alt_splicing.txt", package="kissDE", mustWork=TRUE)
    tableCounts <- read.table(fpath1, head = TRUE)
    conditions <- c("C1", "C1", "C2", "C2")
    diff <- diffExpressedVariants(tableCounts, conditions, nbCore = 2)
    expect_equal(names(diff), c("finalTable", "correctedPVal", "uncorrectedPVal", "resultFitNBglmModel", "f/psiTable", "k2rgFile"))
    expect_equal(dim(diff$finalTable)[2], 13)
    expect_equal(dim(diff$finalTable)[1], length(diff$correctedPVal))
    expect_equal(dim(diff$finalTable)[1], dim(diff$`f/psiTable`)[1])
    expect_equal(names(diff$correctedPVal), names(diff$uncorrectedPVal))
    expect_equal(names(diff$correctedPVal), rownames(diff$resultFitNBglmModel))
    expect_null(diff$k2rgFile)
    expect_equal(dim(diff$finalTable[which(diff$finalTable$Adjusted_pvalue <= 0.05),])[1], 19)
    expect_equal(diff$finalTable[which(diff$finalTable$ID == "event1"), "Deltaf/DeltaPSI"], -0.7824)
  }
})

test_that("pvalue parameter of diffExpressedVariants function works as expected", {
  # test pvalue parameter
  fpath1 <- system.file("extdata/table_counts_alt_splicing.txt", package="kissDE", mustWork=TRUE)
  tableCounts <- read.table(fpath1, head = TRUE)
  conditions <- c("C1", "C1", "C2", "C2")
  diff <- diffExpressedVariants(tableCounts, conditions, pvalue = 0.05)
  expect_equal(names(diff), c("finalTable", "correctedPVal", "uncorrectedPVal", "resultFitNBglmModel", "f/psiTable", "k2rgFile"))
  expect_equal(dim(diff$finalTable)[2], 13)
  expect_equal(dim(diff$finalTable)[1], 19)
  expect_equal(names(diff$correctedPVal), names(diff$uncorrectedPVal))
  expect_equal(names(diff$correctedPVal), rownames(diff$resultFitNBglmModel))
  expect_null(diff$k2rgFile)
  expect_equal(diff$finalTable[which(diff$finalTable$ID == "event1"), "Deltaf/DeltaPSI"], -0.7824)
})
