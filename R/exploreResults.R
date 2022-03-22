exploreResults <- function(rdsFile, k2rgRes=NA) {
  ## rdsFile is outputed by "writeOutputKissDE" with a ".rds" extension
  ## k2rgRes is a kissplice2refgenome file and is not mandatory. By default, we will search the k2rg file in the k2rg field in the rdsFile.
  
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
  
  if(!is.na(k2rgRes)) {
    if(!is.character(k2rgRes)) {
      stop("Input error: 'k2rgRes' must be a character.")
    }
    
    if(!file.exists(k2rgRes)) {
      stop(paste("Input error: 'k2rgRes' \"",k2rgRes, "\" does not exists. If you have moved it, please enter its new localisation with the 'k2rgRes' parameter.", sep=""))
    }
  }
  
  ## Read rds file
  res <- readRDS(rdsFile)
  
  ## Rds file check
  if(!all(names(res)%in%c("finalTable", "correctedPVal", "uncorrectedPVal", "resultFitNBglmModel", "f/psiTable", "k2rgFile", "k2rgRes","dfSS"))) {
    stop(paste("Input error: 'rdsFile' \"",rdsFile, "\" does not correspond to the expected format. Is that a kissDE result rds file?", sep=""))
  }
  
  if(!is.na(k2rgRes)) {
    if(!is.na(res$k2rgRes)) {
      warning(paste("Replacing kissDE+kissplice2refgenome result file \"",res$k2rgRes,"\" by user defined value file \"",k2rgRes,"\". Are you sure that this later file correspond to the \"",rdsFile,"\" rds file?",sep=""))
    } else {
      stop(paste("Input error: a k2rgRes file is provided but 'rdsFile' \"",rdsFile,"\" result from a kissDE run without a k2rg file. You may want to run kissplice2refgenome prior to kissDE to annotate the alternative splicing events."))
    }
    res$k2rgRes <- k2rgRes
    resK2RG <- read.table(res$k2rgRes, sep="\t", comment.char = "", header = T)
  } else {
    if(is.na(res$k2rgRes)) {
      warning(paste("rds file \"",rdsFile,"\" have no kissplice2refgenome file associated with it. You may want to run kissplice2refgenome prior to kissDE to annotate the alternative splicing events.",sep=""))
      resK2RG <- NA
    } else {
      resK2RG <- read.table(res$k2rgRes, sep="\t", comment.char = "", header = T)
    }
  }
  
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
  meanPSI[[C1lab]] <- ifelse(okC1,
                                                   rowMeans(PSItable[,grep(condRepl[1],colnames(PSItable))],na.rm = T),
                                                   NA)
  meanPSI[[C2lab]] <- ifelse(okC2,
                                                   rowMeans(PSItable[,grep(condRepl[2],colnames(PSItable))],na.rm = T),
                                                   NA)
  ### FDR/dPSIs
  diffTable <- res$finalTable[,c("ID","Adjusted_pvalue","Deltaf/DeltaPSI","lowcounts")]
  colnames(diffTable) <- c("ID","Adjusted_pvalue","DeltaPSI","lowcounts")
  rownames(diffTable) <- NULL
  diffTable[,2] <- as.numeric(formatC(diffTable[,2],digits = 1, format = "e"))
  diffTable[,3] <- round(diffTable[,3]*100,1)
  ## Add the meanPSIs
  o <- diffTable$ID # order for after the merge
  diffTable <- merge(meanPSI, diffTable, by=1, all.x=F, all.y=T)
  hideCol <- c("lowcounts")
  
  filterPanelEvents <- ""
  filterPanelBiotypes <- ""
  keepPanelEvents <- ""
  keepPanelBiotypes <- ""
  filterPanelRepeats <- ""
  filterPanelEventsDiff <- ""
  filterPanelBiotypesDiff <- ""
  keepPanelEventsDiff <- ""
  keepPanelBiotypesDiff <- ""
  filterPanelRepeatsDiff <- ""
  
  ### ADDITIONAL INFORMATIONS given by k2rg
  if(!is.na(res$k2rgRes)) {
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
        cAltBloc <- paste(chrom,":",cBloc[1],"-",cBloc[2],sep="")
        j=3
        while(j<length(cBloc)) {
          cAltBloc <- c(cAltBloc, paste(chrom,":",cBloc[j],"-",cBloc[j+1],sep=""))
          j=j+2
        }
        dfSS$alternativeBlocs[i] <- paste(cAltBloc,collapse = "_")
      }
      dfSS <- dfSS[,c("ID","constitutiveBlocs","alternativeBlocs")]
      res$dfSS <- dfSS
      saveRDS(res,rdsFile)
    }
    
    ## Less important info
    dfAddInfo <- resK2RG[,c(16,14,6:8,13,12,18,10,22,23)]
    colnames(dfAddInfo) <- c("ID","ComplexEvent","VariablePartLength","Frameshift","inCDS","Paralogs","upperPathSS","lowerPathSS","unkownSS","SS_IR","NormalisedCounts")
    asNumCompEvents <- as.numeric(dfAddInfo$ComplexEvent[dfAddInfo$ComplexEvent!="-"])
    dfAddInfo$ComplexEvent <- factor(dfAddInfo$ComplexEvent, levels = c("-",as.character(unique(sort(asNumCompEvents)))))
    ## Merge
    PSItable <- merge(dfInfo,PSItable,by=1,all.x=F,all.y=T)
    diffTable <- merge(dfAddInfo,diffTable,by=1,all.x=F,all.y=T)
    diffTable <- merge(dfSS,diffTable,by=1,all.x=F,all.y=T)
    diffTable <- merge(dfInfo,diffTable,by=1,all.x=F,all.y=T)
    hideCol <- colnames(diffTable)[c(5,7:19,24)]
    
    events <- unique(PSItable$EventType)
    events <- events[-grep(",",events)]
    biotypes <- unique(PSItable$Biotype)
    biotypes <- biotypes[-grep(",",biotypes)]
    eventsDiff <- unique(diffTable$EventType)
    eventsDiff <- eventsDiff[-grep(",",eventsDiff)]
    biotypesDiff <- unique(diffTable$Biotype)
    biotypesDiff <- biotypesDiff[-grep(",",biotypesDiff)]
    filterPanelEvents <- selectizeInput("fEvents","Filter this type of event:",choices = events,selected = NULL,multiple=T)
    filterPanelBiotypes <- selectizeInput("fBio","Filter this type of biotype:",choices = biotypes,selected = NULL,multiple=T)
    keepPanelEvents <- selectizeInput("fEventsK","Keep this type of event:",choices = events,selected = NULL,multiple=T)
    keepPanelBiotypes <- selectizeInput("fBioK","Keep this type of biotype:",choices = biotypes,selected = NULL,multiple=T)
    filterPanelRepeats <- checkboxInput("fRepeats","Filter repeat-linked event",F)
    filterPanelEventsDiff <- selectizeInput("fEventsDiff","Filter this type of event:",choices = eventsDiff,selected = NULL,multiple=T)
    filterPanelBiotypesDiff <- selectizeInput("fBioDiff","Filter this type of biotype:",choices = biotypesDiff,selected = NULL,multiple=T)
    keepPanelEventsDiff <- selectizeInput("fEventsDiffK","Keep this type of event:",choices = eventsDiff,selected = NULL,multiple=T)
    keepPanelBiotypesDiff <- selectizeInput("fBioDiffK","Keep this type of biotype:",choices = biotypesDiff,selected = NULL,multiple=T)
    filterPanelRepeatsDiff <- checkboxInput("fRepeatsDiff","Filter repeat-linked event",F)
  }
  diffTable <- diffTable[match(o, diffTable$ID),]
  colDiff <- colnames(diffTable)[-grep("ID|Name|Position|Blocs|Counts|SS$",colnames(diffTable))]
  colDiffContinuous <- colDiff[grep("meanPSI|pvalue|DeltaPSI|Length",colDiff)]
  colDiffDiscrete <- colDiff[-grep("meanPSI|pvalue|DeltaPSI|Length",colDiff)]
  colMeanPSI <- colDiff[grep("meanPSI",colDiff)]
  ## Shiny interface
  # ui
  ui <- fluidPage(
    title="KissDE results",
    
    navbarPage(paste(C1," vs ",C2,sep=""),
               tabPanel("Differential anakysis",
                        fluidRow(
                          column(width=12,
                                 h1(strong(paste("kissDE results (",C1," vs ",C2,")",sep="")),
                                 ),
                                 sidebarLayout(position="right",
                                               sidebarPanel(h2("Filter Events"),
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
                                                            keepPanelEventsDiff,
                                                            filterPanelEventsDiff,
                                                            keepPanelBiotypesDiff,
                                                            filterPanelBiotypesDiff,
                                                            filterPanelRepeatsDiff
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
                                          checkboxGroupInput("fCol", "Hide these columns:", choices = unique(colnames(diffTable)), selected = hideCol,inline = T)
                                        )
                                 ),
                                 column(width=6,
                                        downloadButton("dlDiff","Download the printed Dataset")
                                 )
                          ),
                          column(width=12,
                                 h1("Plot"),
                                 sidebarLayout(position="right",
                                               sidebarPanel(h2("Plot parameters"),
                                                            p("By default, plot a Volcano plot seperated by event type."),
                                                            p("x-axis: discrete or continuous values are allowed. If the 'Plot density' checbox is selected, only continuous values are allowed and multiple choices are allowed."),
                                                            p("y-axis: only continuous values are allowed. If coupled with a discrete x-axis, the density will be printed (box-violin plot)"),
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
                                                            selectInput("plotXgroup","Horizontal wrap:",choices = c("None",colDiffDiscrete), selected = "EventType"),
                                                            selectInput("plotYgroup","Vertical wrap:",choices = c("None",colDiffDiscrete))
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
                                                            filterPanelRepeats
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
      if(!is.na(res$k2rgRes)) {
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
      if(!is.na(res$k2rgRes)) {
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
      }
      cDataDiff$data <- data
      data <- data[,-which(colnames(data)%in%input$fCol)]
      data
    })
    
    
    # display the PSI table
    output$PSItable <- DT::renderDataTable(DT::datatable(
      dataPSI(),
      rownames = FALSE,
      filter = "top",
      options = list(
        dom = 'Bfrtip',
        pageLength = 15,
        scrollX = TRUE
      )
    ))
    
    # display the diff table
    output$diffTable <- DT::renderDataTable(DT::datatable(
      dataDiff(),
      rownames = FALSE,
      filter = "top",
      options = list(
        dom = 'Bfrtip',
        pageLength = 15,
        scrollX = TRUE
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
      if(!is.na(res$k2rgRes)) {
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
      if(!is.na(res$k2rgRes)) {
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
        xLab <- "x"
        yLab <- "y"
        n <- length(input$plotXdensity)
        if(!is.na(res$k2rgRes)) {
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
        data$pch <- ifelse(data$Adjusted_pvalue<=input$plotFDR,"19","1")
        data$col <- ifelse(data$pch=="1","black",
                           ifelse(data$DeltaPSI<=(-input$plotDPSI),"blue",
                                  ifelse(data$DeltaPSI>=input$plotDPSI,"red","black")
                           )
        )
        cols <- c("black"="black","blue"="blue","red"="red")
        pchs <- c("1"=1,"19"=19)
        xLab <- input$plotX
        yLab <- input$plotY
        if(input$plotXtrans!="None") {
          minVal <- min(data[[input$plotX]][data[[input$plotX]]!=0])
          data[[input$plotX]][data[[input$plotX]]==0]=minVal
          data[[input$plotX]] <- log10(data[[input$plotX]])
          if(input$plotXtrans=="-log10") {
            data[[input$plotX]] <- -data[[input$plotX]]
          }
        }
        if(input$plotYtrans!="None") {
          minVal <- min(data[[input$plotY]][data[[input$plotY]]!=0])
          data[[input$plotY]][data[[input$plotY]]==0]=minVal
          data[[input$plotY]] <- log10(data[[input$plotY]])
          if(input$plotYtrans=="-log10") {
            data[[input$plotY]] <- -data[[input$plotY]]
          }
        }
        if(input$plotXgroup!="None") {
          data$wrapX <- data[[input$plotXgroup]]
        }
        if(input$plotYgroup!="None") {
          data$wrapY <- data[[input$plotYgroup]]
        }
        if(input$plotX%in%colDiffDiscrete) {
          g <- ggplot(data = data, aes_string(x=input$plotX,y=input$plotY,fill=input$plotX))+
            theme_bw()+
            theme(strip.background =element_rect(fill="white"),strip.text.y = element_text(angle = 0),legend.position = "none")+
            geom_violin(scale = "width") +
            geom_boxplot(aes(alpha=100),width=0.1,outlier.shape = NA) +
            stat_summary(fun=mean, geom="point", color="black", size=2)+
            stat_summary(fun=median, geom="point", fill="white", shape=23, size=1)
        } else {
          g <- ggplot(data = data, aes_string(x=input$plotX,y=input$plotY,shape="pch",color="col"))+
            theme_bw()+
            geom_point(stroke=0.15)+
            scale_color_manual(values = cols)+
            scale_shape_manual(values=pchs)+
            guides(color="none",shape="none")
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
      list(data,xLab,yLab,g)
    })
    
    output$Plot <- renderPlot({
      Plot()[[4]]
    },height = 1000)
    
    output$selectedPoints <- DT::renderDataTable({
      lData <- Plot()
      data <- lData[[1]]
      #xLab <- lData[[2]]
      #yLab <- lData[[3]]
      dataToShow <- data
      if (!is.null(input$plotBrush)){
        dataToShow <- brushedPoints(data, input$plotBrush)
      } else if (!is.null(input$plotClick)){
        dataToShow <- nearPoints(data, input$plotClick)
      }
      if(!input$plotDensity) {
        dataToShow <- dataToShow[,-which(colnames(dataToShow)%in%input$fCol)]
      }
      DT::datatable(dataToShow, rownames = FALSE,
                    options = list(scrollX = TRUE),
                    caption = 'Select points on the plot to show them in this table.')
    })
  }

  shinyApp(ui, server, options = list(port=3838, host='0.0.0.0'))
}

