install.packages("kissDE_1.0.tar.gz")
library("kissDE")
data<-("pretraitement.txt")
# signifEventsRef <-diffExpressedEvents(2,c(2,2),data)
signifEvents<-diffExpressedEvents(2,c(2,2),data)

# save(signifEventsRef,file="signifEventsRef.RData")

load("signifEventsRef.RData")
signifEventsRef==signifEvents