library("kissDE")


## tests 'kissplice2counts'
fpath <- system.file("extdata", "output_kissplice_SNP.fa", package = "kissDE")
mySNPcounts <- kissplice2counts(fpath, pairedEnd = TRUE)

stopifnot(
  all(names(mySNPcounts) == c("countsEvents", "psiInfo", "discoInfo", "exonicReadsInfo", "dupBcc")),
  all(mySNPcounts$countsEvents[, 1] == mySNPcounts$psiInfo[, 1])
)

## tests 'diffExpressedVariants'
# diffSNP <- diffExpressedVariants(mySNPcounts, mySNPconditions)

res <- diffExpressedVariants(table_counts_alt_splicing, c(rep("condition1", 2), rep("condition2", 2)))

# 
# conditions <- c("C1", "C1", "C2", "C2")
# 
# #### test for input = data.frame ####
# file_df <- "tests/pretraitement.txt"
# counts_df <- read.table(file_df)
# # qualityControl(counts_df, conditions)
# # save(signifVariantsKis_df_save, file = "/home/aurelie/Bureau/save-kissDE/signifVariantsKis_df_save.rdata")
# # load("/home/aurelie/Bureau/save-kissDE/signifVariantsKis_df_save.rdata")
# # identical(signifVariantsKis_df, signifVariantsKis_df_save)
# # names(signifVariantsKis_df)
# 
# signifVariantsKis_df1 <- diffExpressedVariants(counts_df, conditions)
# pvaldf1 <- signifVariantsKis_df1$correctedPVal
# pvalsigndf1 <- which(pvaldf1 < 0.05)
# 
# signifVariantsKis_df2 <- diffExpressedVariants(counts_df, conditions)
# pvaldf2 <- signifVariantsKis_df2$correctedPVal
# pvalsigndf2 <- which(pvaldf2 < 0.05)
# 
# signifVariantsKis_df3 <- diffExpressedVariants(counts_df, conditions)
# pvaldf3 <- signifVariantsKis_df3$correctedPVal
# pvalsigndf3 <- which(pvaldf3 < 0.05)
# 
# all(pvaldf1 == pvaldf2)
# all(pvaldf1 == pvaldf3)
# all(pvaldf2 == pvaldf3)
# 
# all(pvalsigndf1 == pvalsigndf2)
# all(pvalsigndf1 == pvalsigndf3)
# all(pvalsigndf2 == pvalsigndf3)
# 
# all(pvaldf1[pvalsigndf1] == pvaldf2[pvalsigndf2])
# all(pvaldf1[pvalsigndf1] == pvaldf3[pvalsigndf3])
# all(pvaldf2[pvalsigndf2] == pvaldf3[pvalsigndf3])
# 
# 
# ## test with the same seed
# set.seed(40)
# signifVariantsKis_df4 <- diffExpressedVariants(counts_df, conditions)
# set.seed(40)
# signifVariantsKis_df5 <- diffExpressedVariants(counts_df, conditions)
# 
# all(signifVariantsKis_df4$correctedPVal == signifVariantsKis_df5$correctedPVal)
# all(signifVariantsKis_df4$uncorrectedPVal == signifVariantsKis_df5$uncorrectedPVal)
# all(as.vector(na.omit(as.vector(signifVariantsKis_df4$finalTable == signifVariantsKis_df5$finalTable))))
# all(signifVariantsKis_df4$resultFitNBglmModel == signifVariantsKis_df5$resultFitNBglmModel)
#   
# 
# 
# #### test for input = .fa file from KisSplice ####
# file_fa <- "tests/resultsKissSknsh10M.fa"
# counts_fa <- kissplice2counts(file_fa)
# # save(counts_fa, file = "/home/aurelie/Bureau/save-kissDE/counts_fa.RData")
# # load("/home/aurelie/Bureau/save-kissDE/counts_fa.RData")
# # qualityControl(counts_fa, conditions)
# # save(signifVariantsKis_fa, file = "/home/aurelie/Bureau/save-kissDE/signifVariantsKis_fa_save.rdata")
# # load("/home/aurelie/Bureau/save-kissDE/signifVariantsKis_fa_save.rdata")
# # signifVariantsRef == signifVariantsKis
# 
# signifVariantsKis_fa1 <- diffExpressedVariants(counts_fa, conditions)
# pvalfa1 <- signifVariantsKis_fa1$correctedPVal
# pvalsignfa1 <- which(pvalfa1 < 0.05)
# 
# signifVariantsKis_fa2 <- diffExpressedVariants(counts_fa, conditions)
# pvalfa2 <- signifVariantsKis_fa2$correctedPVal
# pvalsignfa2 <- which(pvalfa2 < 0.05)
# 
# signifVariantsKis_fa3 <- diffExpressedVariants(counts_fa, conditions)
# pvalfa3 <- signifVariantsKis_fa3$counts_fa
# pvalsignfa3 <- which(pvalfa3 < 0.05)
# 
# all(pvalsignfa1 == pvalsignfa2)
# all(pvalsignfa1 == pvalsignfa3)
# all(pvalsignfa2 == pvalsignfa3)
# 
# all(pval1[pvalsignfa1] == pval2[pvalsignfa2])
# all(pval1[pvalsignfa1] == pval3[pvalsignfa3])
# all(pval2[pvalsignfa2] == pval3[pvalsignfa3])
# 
# 
# 
# 
# 
# signifVariantsKis_df1 <- diffExpressedVariants(counts_df, conditions)
# signifVariantsKis_df2 <- diffExpressedVariants(counts_df, conditions)
# signifVariantsKis_df3 <- diffExpressedVariants(counts_df, conditions)
# names(signifVariantsKis_df1)
# # [1] "pALLGlobalPhi.glm.nb" "sing.events"          "dataPart3"            "ASSBinfo"             "allEventtables"       "lengths"              "phi"                  "dispData"            
# # # [9] "dispDataMeanCond"  
# # signifVariantsKis_df1$pALLGlobalPhi.glm.nb[2, c(3,  4 , 9, 10, 11, 12)] == signifVariantsKis_df2$pALLGlobalPhi.glm.nb[2,c(3,  4 , 9, 10, 11, 12)]
# # 
# # rbind(signifVariantsKis_df1$pALLGlobalPhi.glm.nb[2, c(3,  4 , 9, 10, 11, 12)],
# #       signifVariantsKis_df2$pALLGlobalPhi.glm.nb[2, c(3,  4 , 9, 10, 11, 12)],
# #       signifVariantsKis_df3$pALLGlobalPhi.glm.nb[2, c(3,  4 , 9, 10, 11, 12)])
# # 
# # # identical(signifVariantsKis_df1$dispData , signifVariantsKis_df2$dispData)
# # # identical(signifVariantsKis_df1$dispDataMeanCond , signifVariantsKis_df2$dispDataMeanCond)
# 
# 
# identical(signifVariantsKis_df1@normalizationFactor , signifVariantsKis_df2@normalizationFactor)
# identical(signifVariantsKis_df1@dispersion, signifVariantsKis_df2@dispersion)
# identical(signifVariantsKis_df1@experimentData, signifVariantsKis_df2@experimentData)
# identical(signifVariantsKis_df1@phenoData, signifVariantsKis_df2@phenoData)
# identical(signifVariantsKis_df1@featureData, signifVariantsKis_df2@featureData)
# identical(signifVariantsKis_df1@annotation, signifVariantsKis_df2@annotation)
# identical(signifVariantsKis_df1@protocolData, signifVariantsKis_df2@protocolData)
# identical(signifVariantsKis_df1@.__classVersion__, signifVariantsKis_df2@.__classVersion__)
# identical(signifVariantsKis_df1@assayData$exprs, signifVariantsKis_df2@assayData$exprs)
# 
# set.seed(40)
# disp1 <- DSS::estDispersion(signifVariantsKis_df1)
# disp2 <- DSS::estDispersion(signifVariantsKis_df2)
# disp3 <- DSS::estDispersion(signifVariantsKis_df3)
# 
# 
# 
# identical(disp1@normalizationFactor , disp2@normalizationFactor)
# identical(disp1@dispersion, disp2@dispersion)
# identical(disp1@experimentData, disp2@experimentData)
# identical(disp1@phenoData, disp2@phenoData)
# identical(disp1@featureData, disp2@featureData)
# identical(disp1@annotation, disp2@annotation)
# identical(disp1@protocolData, disp2@protocolData)
# identical(disp1@.__classVersion__, disp2@.__classVersion__)
# identical(disp1@assayData$exprs, disp2@assayData$exprs)
# 
