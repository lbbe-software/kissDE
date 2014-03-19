install.packages("kissDE_1.0.tar.gz")
library("kissDE")


#### test for input = data.frame ####
data<-read.table("pretraitement.txt")
# signifEventsRef <-diffExpressedEvents(2,c(2,2),data)
signifEventsDF<-diffExpressedEvents(2,c(2,2),data)

# save(signifEventsRef,file="signifEventsRef.RData")

load("signifEventsRef.RData")
signifEventsRef==signifEventsDF


#### test for input = .fa file from KisSplice ####
# signifEventsRef <-diffExpressedEvents(2,c(2,2),data)
signifEventsKis<-diffExpressedEvents(2,c(2,2),"resultsKissSknsh10M.fa")

# save(signifEventsRef,file="signifEventsRef.RData")

load("signifEventsRef.RData")
signifEventsRef==signifEventsKis
