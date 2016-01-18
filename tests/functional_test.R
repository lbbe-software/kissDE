library("kissDE")

conditions <- c("C1", "C1", "C2", "C2")
file_fa <- "tests/resultsKissSknsh10M.fa"
file_df <- "tests/pretraitement.txt"

#### test for input = data.frame ####
counts_df <- read.table(file_df)
# qualityControl(counts_df, conditions)
signifVariantsKis_df <- diffExpressedVariants(counts_df, conditions)
# save(signifVariantsKis_df, file = "signifVariantsKis_df_save.rdata")
# load("signifVariantsKis_df_save.rdata")
# signifVariantsRef == signifVariantsKis

#### test for input = .fa file from KisSplice ####
counts_fa <- kissplice2counts(file_fa)
qualityControl(counts_fa, conditions)
signifVariantsKis_fa <- diffExpressedVariants(counts_fa, conditions)

# save(signifEventsRef,file="signifEventsRef.RData")
# load("signifEventsRef.RData")
# signifVariantsRef == signifVariantsKis
