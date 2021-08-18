context("qualityControl")
test_that("qualityControl work as expected", {
  fpath1 <- system.file("extdata", "output_kissplice_SNV.fa", package = "kissDE")
  mySNVcounts <- kissplice2counts(fpath1, counts = 0, pairedEnd = TRUE)
  mySNVconditions <- c("C1", "C1", "C2", "C2")
  # test qualityControl without storing figures
  qualityControl(mySNVcounts, mySNVconditions, storeFigs = FALSE)
  # test qualityControl storing figures in the default directory
  qualityControl(mySNVcounts, mySNVconditions, storeFigs = TRUE)
  dirname <- tempdir()
  expect_true(file.exists(paste0(dirname, "/kissDEFigures/heatmap.png")))
  expect_true(file.exists(paste0(dirname, "/kissDEFigures/pca.png")))
  unlink(paste0(dirname,"/kissDEFigures"), recursive = TRUE)
  # test qualityControl storing figures in a path choosen by the user
  qualityControl(mySNVcounts, mySNVconditions, storeFigs = "Rplots")
  expect_true(file.exists("Rplots/heatmap.png"))
  expect_true(file.exists("Rplots/pca.png"))
  unlink("Rplots",recursive = T)
})
