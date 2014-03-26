install.packages("../../kissDE_1.0.tar.gz")
library("kissDE")

conditions <- c("C1","C1","C2","C2")

#### test for input = data.frame ####
# counts <- read.table("pretraitement.txt")
# colnames(counts) = c("ID","length",conditions)

#### test for input = .fa file from KisSplice ####
# signifEventsRef <-diffExpressedEvents(2,c(2,2),data)
counts <- kissplice2counts("resultsKissSknsh10M.fa",conditions)

qualityControl(counts)
signifEventsKis<-diffExpressedEvents(counts)

# save(signifEventsRef,file="signifEventsRef.RData")
load("signifEventsRef.RData")
signifEventsRef==signifEventsKis
