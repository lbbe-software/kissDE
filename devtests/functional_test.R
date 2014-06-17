install.packages("../kissDE_0.1.tar.gz")
library("kissDE")

conditions <- c("C1","C1","C2","C2")
file <- ("resultsKissSknsh10M.fa")
#### test for input = data.frame ####
# counts <- read.table("pretraitement.txt")

#### test for input = .fa file from KisSplice ####
# signifEventsRef <-diffExpressedEvents(2,c(2,2),data)
counts <- kissplice2counts(file)

qualityControl(counts,conditions)
signifVariantsKis <- diffExpressedVariants(counts,conditions)

# save(signifEventsRef,file="signifEventsRef.RData")
load("signifEventsRef.RData")
signifVariantsRef == signifVariantsKis
