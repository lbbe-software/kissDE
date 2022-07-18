exploreResults <- function(rdsFile) {
  ## rdsFile is outputed by "writeOutputKissDE" with a ".rds" extension

  ## Input check
  if(is.na(rdsFile)) {
    stop("Input error: 'rdsFile' must be specified.")
  }
  
  if(!is.character(rdsFile)) {
    stop("Input error: 'rdsFile' must be a character.")
  }
  
  if(tail(strsplit(rdsFile,split = "\\.")[[1]],n=1)!="rds") {
    stop(paste("Input error: 'rdsFile' \"",rdsFile, "\" does not have the \".rds\" extension. Is that a kissDE result rds file?" , sep=""))
  }
  
  if(!file.exists(rdsFile)) {
    stop(paste("Input error: 'rdsFile' \"",rdsFile, "\" does not exists.", sep=""))
  }
  
  ## Read rds file
  res <- readRDS(rdsFile)
  
  ## Rds file check
  if(!all(names(res)%in%c("finalTable", "correctedPVal", "uncorrectedPVal", "resultFitNBglmModel", "f/psiTable", "k2rgFile", "k2rgRes","dfSS"))) {
    stop(paste("Input error: 'rdsFile' \"",rdsFile, "\" does not correspond to the expected format. Is that a kissDE result rds file?", sep=""))
  }
  
  resK2RG <- res$k2rgRes
  
  ## Format tables
  # All tables must have these first columns if k2rg is present:
  # bccID geneID geneName Position strand Event Repeat
  # Else, only bccID is mandatory
  
  ### PSIs
  PSItable <- res$`f/psiTable`
  PSItable[,-1] <- round(PSItable[,-1]*100,1)
  conditions <- unlist(lapply(strsplit(colnames(PSItable)[-1],split="_repl"),"[[",1))
  conditions <- factor(conditions, levels = unique(conditions))
  condRepl <- paste(levels(conditions),"_repl\\d+",sep="")
  ### Mean PSIs
  C1 <- levels(conditions)[1]
  C2 <- levels(conditions)[2]
  nC1 <- sum(conditions==C1)
  nC2 <- sum(conditions==C2)
  nC <- nC1+nC2
  minC1 <- round(nC1/2)
  minC2 <- round(nC2/2)
  okC1 <- rowSums(!is.na(PSItable[,grep(condRepl[1],colnames(PSItable))]))>=minC1
  okC2 <- rowSums(!is.na(PSItable[,grep(condRepl[2],colnames(PSItable))]))>=minC2
  meanPSI <- data.frame(ID=PSItable$ID)
  C1lab <- paste("meanPSI.",C1,sep="")
  C2lab <- paste("meanPSI.",C2,sep="")
  C1labVar <- paste("sd.",C1,sep="")
  C2labVar <- paste("sd.",C2,sep="")
  meanPSI[[C1lab]] <- ifelse(okC1,
                                                   rowMeans(PSItable[,grep(condRepl[1],colnames(PSItable))],na.rm = T),
                                                   NA)
  meanPSI[[C1labVar]] <- ifelse(okC1,
                             apply(PSItable[,grep(condRepl[1],colnames(PSItable))],1,sd,na.rm = T),
                             NA)
  meanPSI[[C2lab]] <- ifelse(okC2,
                                                   rowMeans(PSItable[,grep(condRepl[2],colnames(PSItable))],na.rm = T),
                                                   NA)
  meanPSI[[C2labVar]] <- ifelse(okC2,
                                apply(PSItable[,grep(condRepl[2],colnames(PSItable))],1,sd,na.rm = T),
                                NA)
  meanPSI[,-1] <- round(meanPSI[,-1],1)
  ### FDR/dPSIs
  diffTable <- res$finalTable[,c("ID","Adjusted_pvalue","Deltaf/DeltaPSI","lowcounts")]
  colnames(diffTable) <- c("ID","Adjusted_pvalue","DeltaPSI","lowcounts")
  rownames(diffTable) <- NULL
  diffTable[,2] <- as.numeric(formatC(diffTable[,2],digits = 1, format = "e"))
  diffTable[,3] <- round(diffTable[,3]*100,1)
  ## Add the meanPSIs
  o <- diffTable$ID # order for after the merge
  diffTable <- merge(meanPSI, diffTable, by=1, all.x=F, all.y=T)
  showCol <- colnames(diffTable)
  #hideCol <- c("lowcounts")
  
  filterPanelSize <- checkboxInput("plotSize","Adapt the size of the points to the mean event coverage (up to 100)?",T)
  filterPanelEvents <- ""
  filterPanelBiotypes <- ""
  keepPanelEvents <- ""
  keepPanelBiotypes <- ""
  filterPanelRepeats <- ""
  filterPanelGenomicWindow <- ""
  filterPanelChrom <- ""
  filterPanelStart <- ""
  filterPanelEnd <- ""
  filterPanelEventsDiff <- ""
  filterPanelBiotypesDiff <- ""
  keepPanelEventsDiff <- ""
  keepPanelBiotypesDiff <- ""
  filterPanelRepeatsDiff <- ""
  filterPanelGenomicWindowDiff <- ""
  filterPanelChromDiff <- ""
  filterPanelStartDiff <- ""
  filterPanelEndDiff <- ""
  filterPanelCoverageDiff <- numericInput("fCoverDiff","Minimum mean event coverage:",value = 0,min = 0)
  
  dfAddInfoCov <- res$finalTable[,c(1,grep("Variant",names(res$finalTable)))]
  normMeanC1V1 <- rowMeans(dfAddInfoCov[,grep(paste("Variant1_",C1,"_repl",sep=""),names(dfAddInfoCov))],na.rm = T)
  normMeanC2V1 <- rowMeans(dfAddInfoCov[,grep(paste("Variant1_",C2,"_repl",sep=""),names(dfAddInfoCov))],na.rm = T)
  normMeanC1V2 <- rowMeans(dfAddInfoCov[,grep(paste("Variant2_",C1,"_repl",sep=""),names(dfAddInfoCov))],na.rm = T)
  normMeanC2V2 <- rowMeans(dfAddInfoCov[,grep(paste("Variant2_",C2,"_repl",sep=""),names(dfAddInfoCov))],na.rm = T)
  normMeanC1 <- normMeanC1V1+normMeanC1V2
  normMeanC2 <- normMeanC2V1+normMeanC2V2
  C1eC <- paste("EventCoverageMean",C1,sep=".")
  C2eC <- paste("EventCoverageMean",C2,sep=".")
  dfAddInfoCov[[C1eC]] <- round(normMeanC1,1)
  dfAddInfoCov[[C2eC]] <- round(normMeanC2,1)
  dfAddInfoCov[["EventCoverageMean"]] <- round(rowMeans(dfAddInfoCov[,tail(names(dfAddInfoCov),2)]),1)
  dfAddInfoCov <- dfAddInfoCov[,c(1,grep("EventCoverageMean",names(dfAddInfoCov)))]

  ### ADDITIONAL INFORMATIONS given by k2rg
  if(!is.null(res$k2rgFile)) {
    ## Make the required columns
    dfInfo <- resK2RG[,c(16,1,2,3,4,5,9)]
    colnames(dfInfo) <- c("ID","GeneID","GeneName","EventPosition","Strand","EventType","Biotype")
    # Get splice sites info
    if("dfSS"%in%names(res)) {
      dfSS <- res$dfSS
    } else {
      dfSS <- resK2RG[resK2RG$X3.Chromosome_and_genomic_position!="multiple" & !resK2RG$X5.Event_type%in%c("insertion","deletion","indel"),c(16,3,12,18)]
      colnames(dfSS)[1]="ID"
      dfSS$chrom <- gsub("(.*):\\d+-\\d+","\\1",dfSS$X3.Chromosome_and_genomic_position)
      dfSS$mergeSS <- ifelse(dfSS$X12.genomic_position_of_each_splice_site_.upper_path..of_each_SNP=="-",
                             dfSS$X18.genomic_position_of_each_splice_site_.lower_path.,
                             paste(dfSS$X12.genomic_position_of_each_splice_site_.upper_path..of_each_SNP,dfSS$X18.genomic_position_of_each_splice_site_.lower_path.,sep=","))
      dfSS$mergeUniqSS <- unlist(lapply(lapply(lapply(lapply(strsplit(dfSS$mergeSS,split = ","),function(x){t <- table(x); return(names(t)[t==1])}),as.numeric),sort),paste,collapse=","))
      const1 <- paste(gsub("(.*:\\d+)-\\d+","\\1",dfSS$X3.Chromosome_and_genomic_position),
                      gsub("(\\d+),\\d+","\\1",dfSS$X18.genomic_position_of_each_splice_site_.lower_path.),
                      sep="-")
      const2 <- paste(dfSS$chrom,
                      ":",
                      gsub("\\d+,(\\d+)","\\1",dfSS$X18.genomic_position_of_each_splice_site_.lower_path.),
                      gsub(".*:\\d+(-\\d+)","\\1",dfSS$X3.Chromosome_and_genomic_position),
                      sep="")
      dfSS$constitutiveBlocs <- paste(const1,const2,sep="_")
      dfSS$alternativeBlocs <- NA
      for(i in c(1:nrow(dfSS))) {
        d <- dfSS[i,]
        chrom <- d$chrom
        cBloc <- strsplit(d$mergeUniqSS,split = ",")[[1]]
        s1=0
        e1=0
        if(length(grep(cBloc[1],d$constitutiveBlocs))!=0) {
          s1=1
        }
        if(length(grep(cBloc[2],d$constitutiveBlocs))!=0) {
          e1=1
        }
        cAltBloc <- paste(chrom,":",as.numeric(cBloc[1])+s1,"-",as.numeric(cBloc[2])-e1,sep="")
        j=3
        while(j<length(cBloc)) {
          s1=0
          e1=0
          if(length(grep(cBloc[j],d$constitutiveBlocs))!=0) {
            s1=1
          }
          if(length(grep(cBloc[j+1],d$constitutiveBlocs))!=0) {
            e1=1
          }
          cAltBloc <- paste(chrom,":",as.numeric(cBloc[j])+s1,"-",as.numeric(cBloc[j+1])-e1,sep="")
          j=j+2
        }
        dfSS$alternativeBlocs[i] <- paste(cAltBloc,collapse = "_")
      }
      dfSS <- dfSS[,c("ID","constitutiveBlocs","alternativeBlocs")]
      res$dfSS <- dfSS
      saveRDS(res,rdsFile)
    }
    
    ## Less important info
    dfAddInfo <- resK2RG[, c(16, 14, 6:8, 13, 12, 18, 10, 
                             22:length(resK2RG))]
    if (length(dfAddInfo)==12) {
      colnames(dfAddInfo) <- c("ID", "ComplexEvent", "VariablePartLength", 
                               "Frameshift", "inCDS", "Paralogs", "upperPathSS", 
                               "lowerPathSS", "unkownSS", "SS_IR", "Protein ID",
                               "Description (from annotaiton)")
    }
    else {
      colnames(dfAddInfo) <- c("ID", "ComplexEvent", "VariablePartLength", 
                             "Frameshift", "inCDS", "Paralogs", "upperPathSS", 
                             "lowerPathSS", "unkownSS", "SS_IR", "Protein ID",
                             "Description (from annotaiton)", "Busco Rank",
                             "Busco ID", "Busco Score", "Description (from BUSCO)")
    }
    asNumCompEvents <- as.numeric(dfAddInfo$ComplexEvent[dfAddInfo$ComplexEvent!="-"])
    dfAddInfo$ComplexEvent <- factor(dfAddInfo$ComplexEvent, levels = c("-",as.character(unique(sort(asNumCompEvents)))))
    ## Merge
    PSItable <- merge(dfInfo,PSItable,by=1,all.x=F,all.y=T)
    diffTable <- merge(dfAddInfoCov,diffTable,by=1,all.x=F,all.y=T)
    diffTable <- merge(dfAddInfo,diffTable,by=1,all.x=F,all.y=T)
    diffTable <- merge(dfSS,diffTable,by=1,all.x=F,all.y=T)
    diffTable <- merge(dfInfo,diffTable,by=1,all.x=F,all.y=T)
    #hideCol <- colnames(diffTable)[c(5,7:19,21,23,26)]
    #showCol <- colnames(diffTable)[c(1:4,6,20,22,24:25)]
    showCol <- c("ID","GeneID","GeneName","EventPosition","EventType","EventCoverageMean",C1lab,C2lab,"Adjusted_pvalue","DeltaPSI")
    events <- unique(PSItable$EventType)
    events <- events[-grep(",",events)]
    biotypes <- unique(PSItable$Biotype)
    biotypes <- biotypes[-grep(",",biotypes)]
    lPos <- strsplit(PSItable$EventPosition[PSItable$EventPosition!="multiple"],":|-")
    chromPos <- sort(unique(unlist(lapply(lPos,"[[",1))))
    posMin <- 0 # aletrnatively, min(as.integer(unlist(lapply(lPos,"[[",2))))
    posMax <- max(as.integer(unlist(lapply(lPos,"[[",3))))
    eventsDiff <- unique(diffTable$EventType)
    eventsDiff <- eventsDiff[-grep(",",eventsDiff)]
    biotypesDiff <- unique(diffTable$Biotype)
    biotypesDiff <- biotypesDiff[-grep(",",biotypesDiff)]
    lPosDiff <- strsplit(PSItable$EventPosition[diffTable$EventPosition!="multiple"],":|-")
    chromPosDiff <- sort(unique(unlist(lapply(lPosDiff,"[[",1))))
    posMinDiff <- 0 # aletrnatively, min(as.integer(unlist(lapply(lPos,"[[",2))))
    posMaxDiff <- max(as.integer(unlist(lapply(lPosDiff,"[[",3))))
    filterPanelEvents <- selectizeInput("fEvents","Filter this type of event:",choices = events,selected = NULL,multiple=T)
    filterPanelBiotypes <- selectizeInput("fBio","Filter this type of biotype:",choices = biotypes,selected = NULL,multiple=T)
    keepPanelEvents <- selectizeInput("fEventsK","Keep this type of event:",choices = events,selected = NULL,multiple=T)
    keepPanelBiotypes <- selectizeInput("fBioK","Keep this type of biotype:",choices = biotypes,selected = NULL,multiple=T)
    filterPanelRepeats <- checkboxInput("fRepeats","Filter repeat-linked event",F)
    filterPanelGenomicWindow <- checkboxInput("fGenomicWindow","Activate genomic window filter",F)
    filterPanelChrom <- selectizeInput("fChrom","Keep this chromosome:",choices = chromPos,selected = NULL,multiple=T)
    filterPanelStart <- numericInput("fStartPos","Start of the genomic window:",value = posMin,min = posMin,max = posMax)
    filterPanelEnd <- numericInput("fEndPos","End of the genomic window:",value = posMax,min = posMin,max = posMax)
    filterPanelEventsDiff <- selectizeInput("fEventsDiff","Filter this type of event:",choices = eventsDiff,selected = NULL,multiple=T)
    filterPanelBiotypesDiff <- selectizeInput("fBioDiff","Filter this type of biotype:",choices = biotypesDiff,selected = NULL,multiple=T)
    keepPanelEventsDiff <- selectizeInput("fEventsDiffK","Keep this type of event:",choices = eventsDiff,selected = NULL,multiple=T)
    keepPanelBiotypesDiff <- selectizeInput("fBioDiffK","Keep this type of biotype:",choices = biotypesDiff,selected = NULL,multiple=T)
    filterPanelRepeatsDiff <- checkboxInput("fRepeatsDiff","Filter repeat-linked event",F)
    filterPanelChromDiff <- selectizeInput("fChromDiff","Keep this chromosome:",choices = chromPosDiff,selected = NULL,multiple=T)
    filterPanelGenomicWindowDiff <- checkboxInput("fGenomicWindowDiff","Activate genomic window filter",F)
    filterPanelStartDiff <- numericInput("fStartPosDiff","Start of the genomic window:",value = posMinDiff,min = posMinDiff,max = posMaxDiff)
    filterPanelEndDiff <- numericInput("fEndPosDiff","End of the genomic window:",value = posMaxDiff,min = posMinDiff,max = posMaxDiff)
  } else {
    # NO K2RG FILE
    # We still have info about event coverage
    diffTable <- merge(dfAddInfoCov,diffTable,by=1,all.x=F,all.y=T)
    showCol <- c("ID","EventCoverageMean",C1lab,C2lab,"Adjusted_pvalue","DeltaPSI")
  }
  diffTable <- diffTable[match(o, diffTable$ID),]
  colDiff <- colnames(diffTable)[-grep("ID|Name|Position|Blocs|Counts|SS$",colnames(diffTable))]
  colDiffContinuous <- colDiff[grep("EventCoverage|meanPSI|pvalue|DeltaPSI|Length",colDiff)]
  colDiffDiscrete <- colDiff[-grep("EventCoverage|meanPSI|pvalue|DeltaPSI|Length",colDiff)]
  colMeanPSI <- colDiff[grep("meanPSI",colDiff)]
  ## Shiny interface
  # ui
  ui <- fluidPage(
    title="KissDE results",
    
    navbarPage(paste(C1," vs ",C2,sep=""),
               tabPanel("Differential analysis",
                        fluidRow(
                          column(width=12,
                                 h1(strong(paste("kissDE results (",C1," vs ",C2,")",sep="")),
                                 ),
                                 sidebarLayout(position="right",
                                               sidebarPanel(
                                                 h2("Filter Events"),
                                                            width = 3,
                                                            
                                                            numericInput("fFDR",
                                                                         label = "Maximum adjusted p-value (FDR):",
                                                                         min = 0,
                                                                         max = 1,
                                                                         value = 1,
                                                                         step = 0.001),
                                                            checkboxInput("fNADPSI","Filter events with an NA deltaPSI",T),
                                                            checkboxInput("fAbsDPSI","Use absolute value for dPSI filter",F),
                                                            numericInput("fDPSI",
                                                                         label = "Minimum deltaPSI value (maximum value if negative)",
                                                                         min = -100,
                                                                         max = 100,
                                                                         value = 0,
                                                                         step = 1),
                                                            sliderInput("fPSIC1",
                                                                        label = paste("PSI range for condition ",C1,sep=""),
                                                                        min = 0,
                                                                        max = 100,
                                                                        value = c(0,100),
                                                                        step = 1),
                                                            sliderInput("fPSIC2",
                                                                        label = paste("PSI range for condition ",C2,sep=""),
                                                                        min = 0,
                                                                        max = 100,
                                                                        value = c(0,100),
                                                                        step = 1),
                                                 filterPanelCoverageDiff,
                                                            keepPanelEventsDiff,
                                                            filterPanelEventsDiff,
                                                            keepPanelBiotypesDiff,
                                                            filterPanelBiotypesDiff,
                                                            filterPanelRepeatsDiff,
                                                 h2("Select Genomic Window"),
                                                 filterPanelGenomicWindowDiff,
                                                 filterPanelChromDiff,
                                                 filterPanelStartDiff,
                                                 filterPanelEndDiff
                                               ),
                                               mainPanel(width=9,
                                                         tagList(
                                                           DT::dataTableOutput('diffTable')
                                                         )
                                               )
                                 )
                          ),
                          column(width=12,
                                 column(width=6,
                                        wellPanel(
                                          checkboxGroupInput("fCol", "Show these columns:", choices = unique(colnames(diffTable)), selected = showCol,inline = T)
                                        )
                                 ),
                                 column(width=6,
                                        downloadButton("dlDiff","Download this Dataset")
                                 )
                          ),
                          column(width=12,
                                 h1("Plot"),
                                 sidebarLayout(position="right",
                                               sidebarPanel(h2("Plot parameters"),
                                                            p("By default, plot a Volcano plot seperated by event type."),
                                                            p("x-axis: discrete or continuous values are allowed. If the 'Plot density' checbox is selected, only continuous values are allowed and multiple choices are allowed."),
                                                            p("y-axis: only continuous values are allowed. If coupled with a discrete x-axis, the density will be printed (box-violin plot)"),
                                                            p("In plots with two continuous values, the bigger the points, the higher the mean coverage of the event (EventCoverageMean column, the max is reached at 100)."),
                                                            width=3,
                                                            checkboxInput("plotDensity","Plot box-violin plot of x-axis continuous values?",F),
                                                            conditionalPanel("input.plotDensity",
                                                                             selectInput("plotXdensity","Choose the values to print on the x-axis:",choices = colDiffContinuous,multiple = T,selected = colMeanPSI)
                                                            ),
                                                            conditionalPanel("!input.plotDensity",
                                                                             selectInput("plotX","Choose the values to print on the x-axis:",choices = colDiff,multiple = F,selected = "DeltaPSI"),
                                                                             selectInput("plotY","Choose the values to print on the y-axis:",choices = colDiffContinuous, multiple = F,selected = "Adjusted_pvalue")),
                                                            radioButtons("plotXtrans","Transformation to use for the x-axis:",choices = c("None","log10","-log10"),selected = "None"),
                                                            radioButtons("plotYtrans","Transformation to use for the y-axis:",choices = c("None","log10","-log10"),selected = "-log10"),
                                                            numericInput("plotFDR","Maximum adjusted p-value (FDR) to highlight a point:",min = 0,max = 1,value = 0.05,step = 0.001),
                                                            numericInput("plotDPSI","Minimum absolute deltaPSI value to highlight a point:",min = 0,max = 100,value = 10,step=1),
                                                            selectInput("plotXgroup","Wrap 1:",choices = c("None",colDiffDiscrete), selected = "EventType"),
                                                            selectInput("plotYgroup","Wrap 2:",choices = c("None",colDiffDiscrete)),
                                                            filterPanelSize
                                                            
                                               ),
                                               mainPanel(width=9,
                                                         actionButton("aPlot","Update plot with the data from the above table"),
                                                         withSpinner(plotOutput('Plot',
                                                                                click = "plotClick",
                                                                                brush = brushOpts(
                                                                                  id = "plotBrush"
                                                                                )
                                                         )
                                                         ),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         br(),
                                                         DT::dataTableOutput('selectedPoints')
                                               )
                                 )
                          )
                        )
                        
               ),
               
               
               tabPanel("PSI table",
                        fluidRow(
                          column(width=12,
                                 h1(strong("kissDE PSIs for each ASE")),
                                 sidebarLayout(position="right",
                                               sidebarPanel(h2("Filter Events"),
                                                            width = 3,
                                                            
                                                            checkboxInput("fCompleteCases", "Only show complete cases (filter ASE with one or more NA PSI value)", F),
                                                            keepPanelEvents,
                                                            filterPanelEvents,
                                                            keepPanelBiotypes,
                                                            filterPanelBiotypes,
                                                            filterPanelRepeats,
                                                            h2("Select Genomic Window"),
                                                            filterPanelGenomicWindow,
                                                            filterPanelChrom,
                                                            filterPanelStart,
                                                            filterPanelEnd
                                               ),
                                               mainPanel(width=9,
                                                         tagList(
                                                           DT::dataTableOutput('PSItable')
                                                         )
                                               )
                                 )
                                 
                          ),
                          column(width=6,
                                 downloadButton("dlPSI","Download the printed Dataset")
                          )
                        ),
                        column(width=12,
                               h1("PCA analysis on the selected lines"),
                               sidebarLayout(position="right",
                                             sidebarPanel(h2("PCA parameters"),
                                                          width=3,
                                                          numericInput("PCA1","Choose the PC for the x-axis:",value = 1,min = 1,max=nC-1,step = 1),
                                                          uiOutput("PCA2")
                                             ),
                                             mainPanel(width=9,
                                                       radioButtons("PCAtype","Choose the input data for the PCA:",choiceNames = c("All ASE","n most variables ASE"), choiceValues = c("all","n")),
                                                       conditionalPanel("input.PCAtype == 'n'",
                                                                        numericInput("PCAn",label = "Number of most variable ASE to use:",value = 500,min=2,step=1)
                                                       ),
                                                       actionButton("aPCA","Update PCA"),
                                                       withSpinner(plotOutput('PCAplot',
                                                                              click = "plotClickPCA",
                                                                              brush = brushOpts(
                                                                                id = "plotBrushPCA"
                                                                              )
                                                       )
                                                       ),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       DT::dataTableOutput('selectedPointsPCA')
                                             )
                               )
                        )
               )
    )
    
  )
  
  server <- function(input, output) {
    
    PSIcol <- grep("_repl\\d+",colnames(PSItable))
    
    cDataPSI <- reactiveValues(data = NULL)
    cDataDiff <- reactiveValues(data = NULL)
    
    # Keep only the selected PSI lines
    dataPSI <- reactive({
      data <- PSItable
      if(input$fCompleteCases) {
        data <- data[complete.cases(data[,PSIcol]),]
      }
      if(!is.null(res$k2rgFile)) {
        if(!is.null(input$fEventsK)) {
          grepEvents <- paste("\\b",paste(input$fEventsK,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[grep(grepEvents,data$EventType),]
        }
        if(!is.null(input$fEvents)) {
          grepEvents <- paste("\\b",paste(input$fEvents,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[-grep(grepEvents,data$EventType),]
        }
        if(!is.null(input$fBioK)) {
          grepEvents <- paste("\\b",paste(input$fBioK,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[grep(grepEvents,data$Biotype),]
        }
        if(!is.null(input$fBio)) {
          grepEvents <- paste("\\b",paste(input$fBio,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[-grep(grepEvents,data$Biotype),]
        }
        if(input$fRepeats) {
          data <- data[data$EventPosition!="multiple",]
        }
        if(input$fGenomicWindow) {
          if(!is.null(input$fChrom)) {
            data <- data[data$EventPosition!="multiple",]
            lChrom <- unlist(lapply(strsplit(data$EventPosition,":"),"[[",1))
            data <- data[lChrom%in%input$fChrom,]
          }
          if(input$fStartPos!=posMin | input$fEndPos!=posMax) {
            data <- data[data$EventPosition!="multiple",]
            lPos <- strsplit(data$EventPosition,":|-")
            data <- data[as.integer(unlist(lapply(lPos,"[[",2)))>=input$fStartPos & as.integer(unlist(lapply(lPos,"[[",3)))<=input$fEndPos,]
          }
        }
      }
      cDataPSI$data <- data
      data
    })
    
    # Keep only the selected diff lines
    dataDiff <- reactive({
      data <- diffTable[as.numeric(diffTable$Adjusted_pvalue)<=input$fFDR,]
      if(input$fNADPSI | input$fDPSI!=0) {
        data <- data[!is.na(data$DeltaPSI),]
      }
      if(input$fDPSI!=0) {
        if(input$fAbsDPSI) {
          data <- data[abs(data$DeltaPSI) >= abs(input$fDPSI),]
        } else {
          if(input$fDPSI<0) {
            data <- data[data$DeltaPSI <= input$fDPSI,]
          } else {
            data <- data[data$DeltaPSI >= input$fDPSI,]
          }
        }
      }
      if(input$fPSIC1[1]!=0) {
        data <- data[!is.na(data[[C1lab]]) & data[[C1lab]] >= input$fPSIC1[1] & data[[C1lab]] <= input$fPSIC1[2],]
      }
      if(input$fPSIC2[1]!=0) {
        data <- data[!is.na(data[[C2lab]]) & data[[C2lab]] >= input$fPSIC2[1] & data[[C2lab]] <= input$fPSIC2[2],]
      }
      data <- data[data$EventCoverageMean>=input$fCoverDiff,]
      if(!is.null(res$k2rgFile)) {
        if(!is.null(input$fEventsDiffK)) {
          grepEvents <- paste("\\b",paste(input$fEventsDiffK,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[grep(grepEvents,data$EventType),]
        }
        if(!is.null(input$fEventsDiff)) {
          grepEvents <- paste("\\b",paste(input$fEventsDiff,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[-grep(grepEvents,data$EventType),]
        }
        if(!is.null(input$fBioDiffK)) {
          grepEvents <- paste("\\b",paste(input$fBioDiffK,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[grep(grepEvents,data$Biotype),]
        }
        if(!is.null(input$fBioDiff)) {
          grepEvents <- paste("\\b",paste(input$fBioDiff,collapse = "\\b|\\b"),"\\b",sep="")
          data <- data[-grep(grepEvents,data$Biotype),]
        }
        if(input$fRepeatsDiff) {
          data <- data[data$EventPosition!="multiple",]
        }
        if(input$fGenomicWindowDiff) {
          if(!is.null(input$fChromDiff)) {
            data <- data[data$EventPosition!="multiple",]
            lChrom <- unlist(lapply(strsplit(data$EventPosition,":"),"[[",1))
            data <- data[lChrom%in%input$fChromDiff,]
          }
          if(input$fStartPosDiff!=posMin | input$fEndPosDiff!=posMax) {
            data <- data[data$EventPosition!="multiple",]
            lPos <- strsplit(data$EventPosition,":|-")
            data <- data[as.integer(unlist(lapply(lPos,"[[",2)))>=input$fStartPosDiff & as.integer(unlist(lapply(lPos,"[[",3)))<=input$fEndPosDiff,]
          }
        }
      }
      cDataDiff$data <- data
      data <- data[,which(colnames(data)%in%input$fCol)]
      data
    })
    
    
    # display the PSI table
    output$PSItable <- DT::renderDataTable(DT::datatable(
      dataPSI(),
      extensions = c('Buttons','ColReorder'),
      rownames = FALSE,
      filter = "top",
      options = list(
        dom = 'Bfrtip',
        pageLength = 20,
        buttons = c('copy','csv','excel'),
        scrollX = TRUE,
        colReorder = TRUE
      )
    ))
    
    # display the diff table
    output$diffTable <- DT::renderDataTable(DT::datatable(
      dataDiff(),
      extensions = c('Buttons','ColReorder'),
      rownames = FALSE,
      filter = "top",
      options = list(
        dom = 'Bfrtip',
        buttons = c('copy','csv','excel'),
        pageLength = 20,
        scrollX = TRUE,
        colReorder = TRUE
      )
    ))
    
    # Downloads
    output$dlPSI <- downloadHandler(
      filename = paste("kissDE.PSI.",C1,"_vs_",C2,".tab",sep=""),
      content = function(file) {
        write.table(dataPSI(),file,sep="\t",quote=F,row.names = F)
      }
    )
    
    output$dlDiff <- downloadHandler(
      filename = paste("kissDE.results.",C1,"_vs_",C2,".tab",sep=""),
      content = function(file) {
        write.table(dataDiff(),file,sep="\t",quote=F,row.names = F)
      }
    )
    
    
    ## PCA
    output$PCA2 <- renderUI({
      n1 <- input$PCA1
      n2 <- input$PCA2
      if(is.null(n2)) {
        n2 <- 2
      }
      numericInput("PCA2","Choose the PC for the y-axis:",value = ifelse(n2>n1,n2,n1),min = n1,max=nC-1,step = 1)
    })
    
    PCAdata <- reactiveValues(data = NULL,
                              var = NULL,
                              contrib = NULL)
    
    observeEvent(input$aPCA,{
      cDataPSI <- cDataPSI$data
      m <- as.matrix(cDataPSI[,PSIcol])
      rownames(m) <- cDataPSI$ID
      if(!is.null(res$k2rgFile)) {
        rownames(m) <- paste(cDataPSI$ID,cDataPSI$GeneName,cDataPSI$EventPosition,cDataPSI$EventType,sep="@")
      }
      if(input$PCAtype=="all") {
        n <- 0
      } else {
        n <- input$PCAn
        if(n<2) {
          PCAdata$data <- NULL
          PCAdata$var <- NULL
          PCAdata$contrib <- NULL
          return()
        }
      }
      pca <- .ShinyPCA(m,1,nC-1,n)
      data <- pca$li
      data$group <- conditions
      data$label <- rownames(data)
      PCAdata$data <- data
      PCAdata$var <- pca$eig/sum(pca$eig)
      PCAdata$contrib <- data.frame(round(get_pca_var(pca)$contrib*100,2)) # Contribution de chaque gene a PC
    })
    
    PCAplot <- reactive({
      if(is.null(PCAdata$data)) {
        return()
      }
      data <- PCAdata$data
      pcaVar <- PCAdata$var
      a1 <- input$PCA1
      a2 <- input$PCA2
      data <- data.frame(x=data[[paste("Axis",a1,sep="")]],
                         y=data[[paste("Axis",a2,sep="")]],
                         group=data$group,
                         label=data$label)
      g<-ggplot(data=data, aes_string(x="x", y="y", colour="group", label="label"))+
        geom_hline(yintercept = 0,size=0.3)+
        geom_vline(xintercept = 0,size=0.3)+
        geom_point(size=3,stroke=0.2)+
        labs(x=paste0("PC",a1," = ",round(pcaVar[a1]*100, 0),"%",sep=""),y=paste0("PC",a2," = ",round(pcaVar[a2]*100, 0),"%",sep=""))
      return(list(data,g))
    })
    
    output$PCAplot <- renderPlot({
      PCAplot()[[2]]
    },
    height = 1000)
    
    output$selectedPointsPCA <- DT::renderDataTable({
      data <- PCAplot()[[1]]
      dataToShow <- data
      if (!is.null(input$plotBrushPCA)){
        dataToShow <- brushedPoints(data, input$plotBrushPCA)
      } else if (!is.null(input$plotClickPCA)){
        dataToShow <- nearPoints(data, input$plotClickPCA)
      }
      DT::datatable(dataToShow, rownames = FALSE,
                    options = list(scrollX = TRUE),
                    caption = 'Select points on the plot to show them in this table.')
    })
    
    ## Plot
    plotData <- reactiveValues(data = NULL)
    
    observeEvent(input$aPlot,{
      data <- cDataDiff$data[input[["diffTable_rows_all"]],]
      if(!is.null(res$k2rgFile)) {
        data$EventType[grep(",",data$EventType)] <- unlist(lapply(strsplit(data$EventType[grep(",",data$EventType)],","),"[[",1))
        data$Biotype[grep("protein_coding",data$Biotype)] <- "protein_coding"
        data$Biotype[grep(",",data$Biotype)] <- unlist(lapply(strsplit(data$Biotype[grep(",",data$Biotype)],","),"[[",1))
        data$Strand[grep(",",data$Strand)] <- unlist(lapply(strsplit(data$Strand[grep(",",data$Strand)],","),"[[",1))
      }
      plotData$data <- data
    })
    
    Plot <- reactive({
      if(is.null(plotData$data)) {
        return()
      }
      cdata <- plotData$data
      if(input$plotDensity) {
        xName <- "x"
        yName <- "y"
        n <- length(input$plotXdensity)
        if(!is.null(res$k2rgFile)) {
          data <- data.frame(ID=rep(cdata$ID,n),
                             GeneName=rep(cdata$GeneName,n),
                             EventPosition=rep(cdata$EventPosition,n),
                             EventType=rep(cdata$EventType,n))
        } else {
          data <- data.frame(ID=rep(cdata$ID,n))
        }
        data$x <- rep(c(input$plotXdensity),rep(nrow(cdata),n))
        data$y <- unlist(cdata[,input$plotXdensity])
        if(input$plotYtrans!="None") {
          minVal <- min(data$y[data$y!=0])
          data$y[data$y==0]=minVal
          data$y <- log10(data$y)
          if(input$plotXtrans=="-log10") {
            data$y <- -data$y
          }
        }
        if(input$plotXgroup!="None") {
          data$wrapX <- unlist(cdata[,input$plotXgroup])
        }
        if(input$plotYgroup!="None") {
          data$wrapY <- unlist(cdata[,input$plotYgroup])
        }
        g<-ggplot(data = data, aes_string(x="x",y="y",fill="x"))+
          theme_bw()+
          theme(strip.background =element_rect(fill="white"),strip.text.y = element_text(angle = 0),legend.position = "none")+
          geom_violin(scale = "width") +
          geom_boxplot(aes(alpha=100),width=0.1,outlier.shape = NA) +
          stat_summary(fun=mean, geom="point", color="black", size=2)+
          stat_summary(fun=median, geom="point", fill="white", shape=23, size=1)
      } else {
        data <- cdata
        doSize <- F
        if(input$plotSize & min(data[["EventCoverageMean"]])<100) {
          doSize <- T
          alphaEstimation <- data[["EventCoverageMean"]]/100
          data$alpha <- ifelse(alphaEstimation>1,1,alphaEstimation)
          data$size <- data$alpha*2
        }
        data$pch <- ifelse(data$Adjusted_pvalue<=input$plotFDR,"19","1")
        data$col <- ifelse(data$pch=="1","black",
                           ifelse(data$DeltaPSI<=(-input$plotDPSI),"blue",
                                  ifelse(data$DeltaPSI>=input$plotDPSI,"red","black")
                           )
        )
        cols <- c("black"="black","blue"="blue","red"="red")
        pchs <- c("1"=1,"19"=19)
        xName <- input$plotX
        yName <- input$plotY
        if(input$plotXtrans!="None") {
          xName <- paste("log10(",xName,")",sep="")
          dataX <- data[[input$plotX]]
          minVal <- min(dataX[dataX!=0])
          dataX[dataX==0] <- minVal
          dataX <- log10(dataX)
          if(input$plotXtrans=="-log10") {
            xName <- paste("-",xName,sep="")
            dataX <- (-1)*dataX
          }
          data$x <- dataX
        } else {
          data$x <- data[[xName]]
        }
        if(input$plotYtrans!="None") {
          yName <- paste("log10(",yName,")",sep="")
          dataY <- data[[input$plotY]]
          minVal <- min(dataY[dataY!=0])
          dataY[dataY==0] <- minVal
          dataY <- log10(dataY)
          if(input$plotYtrans=="-log10") {
            yName <- paste("-",yName,sep="")
            dataY <- (-1)*dataY
          }
          data$y <- dataY
        } else {
          data$y <- data[[yName]]
        }
        if(input$plotXgroup!="None") {
          data$wrapX <- data[[input$plotXgroup]]
        }
        if(input$plotYgroup!="None") {
          data$wrapY <- data[[input$plotYgroup]]
        }
        if(input$plotX%in%colDiffDiscrete) {
          g <- ggplot(data = data, aes_string(x="x",y="y",fill=input$plotX))+
            theme_bw()+
            theme(strip.background =element_rect(fill="white"),strip.text.y = element_text(angle = 0),legend.position = "none")+
            geom_violin(scale = "width") +
            geom_boxplot(aes(alpha=100),width=0.1,outlier.shape = NA) +
            stat_summary(fun=mean, geom="point", color="black", size=2)+
            xlab(xName)+
            ylab(yName)+
            stat_summary(fun=median, geom="point", fill="white", shape=23, size=1)
        } else {
          if(doSize) {
            g <- ggplot(data = data, aes_string(x="x",y="y",shape="pch",color="col",alpha="alpha",size="size"))+
              theme_bw()+
              geom_point(stroke=0.15)+
              scale_color_manual(values = cols)+
              scale_shape_manual(values=pchs)+
              scale_alpha(range = c(0,1))+
              scale_size(range=c(0,2))+
              xlab(xName)+
              ylab(yName)+
              guides(color="none",shape="none",alpha="none",size="none")
          } else {
            g <- ggplot(data = data, aes_string(x="x",y="y",shape="pch",color="col"))+
              theme_bw()+
              geom_point(stroke=0.15)+
              scale_color_manual(values = cols)+
              scale_shape_manual(values=pchs)+
              xlab(xName)+
              ylab(yName)+
              guides(color="none",shape="none")
          }
        }
      }
      if(input$plotXgroup!="None") {
        if(input$plotYgroup!="None") {
          g <- g+
            facet_wrap(wrapY~wrapX)
        } else {
          g <- g+
            facet_wrap(.~wrapX)
        }  
      } else if(input$plotYgroup!="None") {
        g <- g+
          facet_wrap(wrapY~.)
      }
      list(data,xName,yName,g)
    })
    
    output$Plot <- renderPlot({
      Plot()[[4]]
    },height = 1000)
    
    output$selectedPoints <- DT::renderDataTable({
      lData <- Plot()
      data <- lData[[1]]
      xName <- lData[[2]]
      yName <- lData[[3]]
      dataToShow <- data
      if (!is.null(input$plotBrush)){
        dataToShow <- brushedPoints(data, input$plotBrush)
      } else if (!is.null(input$plotClick)){
        dataToShow <- nearPoints(data, input$plotClick)
      }
      if(!input$plotDensity) {
        dataToShow <- dataToShow[,which(colnames(dataToShow)%in%c(input$fCol))]
      }
      DT::datatable(dataToShow, rownames = FALSE,
                    extensions = c('Buttons','ColReorder'),
                    options = list(scrollX = TRUE,
                                   dom = 'Bfrtip',
                                   colReorder = TRUE,
                                   buttons =c('copy','csv','excel')),
                    caption = 'Select points on the plot to show them in this table.')
    })
  }

  shinyApp(ui, server, options = list(port=3838, host='0.0.0.0'))
}

