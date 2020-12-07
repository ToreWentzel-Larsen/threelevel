# background functions used in the shiny app
cheungdek <- function(metframe, rmaobj) {
  n <- length(metframe$v)
  list.inverse.variances <- 1 / (metframe$v)
  sum.inverse.variances <- sum(list.inverse.variances)
  squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2
  list.inverse.variances.square <- 1 / (metframe$v^2)
  sum.inverse.variances.square <-
    sum(list.inverse.variances.square)
  numerator <- (n - 1) * sum.inverse.variances
  denominator <- squared.sum.inverse.variances -
    sum.inverse.variances.square
  estimated.sampling.variance <- numerator / denominator
  I2_1 <- (estimated.sampling.variance) / (rmaobj$sigma2[1]
                                           + rmaobj$sigma2[2] + estimated.sampling.variance)
  I2_2 <- (rmaobj$sigma2[1]) / (rmaobj$sigma2[1]
                                + rmaobj$sigma2[2] + estimated.sampling.variance)
  I2_3 <- (rmaobj$sigma2[2]) / (rmaobj$sigma2[1]
                                + rmaobj$sigma2[2] + estimated.sampling.variance)
  amountvariancelevel1 <- I2_1 * 100
  amountvariancelevel2 <- I2_2 * 100
  amountvariancelevel3 <- I2_3 * 100
  res <- c(amountvariancelevel1,amountvariancelevel2,amountvariancelevel3)
  names(res) <- c("level1","level2","level3")
  return(list(rma=rmaobj, varcomp=res))
} # end function cheungdek

awf <- function(frame1, dec1=2) {
  overall <- rma.mv(y, v, random = list(~ 1 | effectsizeID, ~ 1 | studyID), 
                    tdist=TRUE, data=frame1)
  oeff <- c(overall$b,overall$ci.lb,overall$ci.ub,overall$pval)
  names(oeff) <- c("effect","lower",
                   "upper","p-value")
  oeffd <- oeff # choose number of decimals
  oeffd[1:3] <- round(oeff[1:3],dec1)
  oeffd[4] <- round(oeff[4],3)
  oeffd <- as.character(oeffd)
  names(oeffd) <- names(oeff)
  if (oeff[4]<.001) oeffd[4] <- "<.001"
  seo <- sqrt(overall$sigma2) # standard error
  owithin <- update(overall, sigma2=c(0,NA))
  obetween <- update(overall, sigma2=c(NA,0))
  pwithin <- anova(overall,owithin)$pval
  pbetween <- anova(overall,obetween)$pval
  varcomp <- data.frame(standarderror=seo, plikrat=c(pwithin,pbetween),
                        row.names=c("within","between"))
  varcompd <- data.frame(domail=c("within","between"), 
                         standarderror=varcomp$standarderror, 
                         plikrat=varcomp$plikrat)
  varcompd$standarderror <- round(varcomp$standarderror,dec1)
  varcompd$plikrat <- round(varcomp$plikrat,3)
  varcompd$standarderror <- as.character(varcompd$standarderror) 
  varcompd$plikrat <- as.character(varcompd$plikrat)
  if (varcomp[1,2]<.001) varcompd[1,3] <- "<.001" 
  if (varcomp[2,2]<.001) varcompd[2,3] <- "<.001" 
  acheung <- cheungdek(frame1, overall)$varcomp
  names(acheung) <- c("level 1","level 2","level 3")
  acheungd <- as.character(round(acheung,dec1))
  names(acheungd) <- names(acheung)
  return(list(oeff=oeff, oeffd=oeffd, varcomp=varcomp,
              varcompd=varcompd, cheung=acheung, cheungd=acheungd))
} # end function awf

# moderator analysis, one categorical or continuous variable
mod1f <- function(frame1, covnr=5, dec1=2, centr1=0) {
  if (covnr<5) stop("the moderators begin at nr 5")
  names1 <- names(frame1)[covnr]
  names(frame1)[covnr] <- "covact" # temporary name
  classc <- class(frame1[,covnr])
  if (!(classc %in% c("factor","numeric","integer"))) 
    stop("the variable has to be categorical or continuous")
  if (classc=="factor") {
    catc <- levels(frame1$covact)
    ncatc <- length(catc)
    nest <- ncatc*(ncatc+1)/2
    estnr <- 1:nest
    estnames <- catc
    estgr <- rep(0,ncatc)
    for (levnr in 1:(ncatc-1)) {
      estnames <- c(estnames, paste0(catc[(levnr+1):ncatc]," - ", catc[levnr]))
      estgr <- c(estgr,rep(levnr,ncatc-levnr))
    }
    estframe <- data.frame(name=estnames)
    estframe$estimate <- rep(NA,nest)
    estframe$standarderror <- rep(NA,nest)
    estframe$lower <- rep(NA,nest)
    estframe$upper <- rep(NA,nest)
    estframe$p.value <- rep(NA,nest)
    msemod <- rep(NA,2)
    for (catnr in 1:ncatc) {
      levc <- catc[catnr]
      modc <- rma.mv(y, v, mods = ~ relevel(covact, ref=levc), tdist=TRUE, data=frame1,
                     random = list(~ 1 | effectsizeID, ~ 1 | studyID))
      if (catnr==1) msemod <- sqrt(modc$sigma2)
      estframe[catnr,-1] <- data.frame(modc$b[1], modc$se[1], modc$ci.lb[1], 
                                       modc$ci.ub[1], modc$pval[1])
      moduse <- (catnr+1):ncatc
      if (catnr<ncatc) estframe[estgr==catnr,-1] <-  
        data.frame(modc$b[moduse], modc$se[moduse], modc$ci.lb[moduse], 
                   modc$ci.ub[moduse], modc$pval[moduse])
    }
    modr <- modc
  } # end computations for categorical variable
  if (classc!="factor") {
    frame1$contvar <- frame1[,covnr] - centr1
    modk <- rma.mv(y, v, mods = ~ I(contvar-centr1), tdist=TRUE, data=frame1,
                   random = list(~ 1 | effectsizeID, ~ 1 | studyID))
    msemod <- sqrt(modk$sigma2)
    estframe <- data.frame(
      navn=c("predik (centr)","slope"),
      estimate=modk$b, standarderror=modk$se, lower=modk$ci.lb,
      upper=modk$ci.ub, p.value=modk$pval)
    modr <- modk
  } # end computations for continuous variable
  names(msemod) <-c("within","between")
  msemodd <- as.character(round(msemod,dec1))
  names(msemodd) <- names(msemod)
  hetonlyp <- 1
  if (classc=="factor") if (ncatc>2) hetonlyp <- 0
  if (hetonlyp==1) {
    mpval <- c(modr$QEp)
    names(mpval) <- c("heterogeneity")
    mpvald <- as.character(round(mpval,3))
    names(mpvald) <- c("heterogeneity")
    mpvald[mpval<.001] <- "<0.001"
  } else {
    mpval <- c(modr$QEp, modr$QMp)
    names(mpval) <- c("heterogen","moderator")
    mpvald <- as.character(round(mpval,3))
    names(mpvald) <- c("heterogeneity","moderator")
    mpvald[mpval<.001] <- "<0.001"
  } # end p-values for heterogeneity and the moderator
  names(estframe)[1] <- names1
  estframed <- estframe
  estframed$estimate <- as.character(round(estframe$estimate, dec1))
  estframed$standarderror <- as.character(round(estframe$standarderror, dec1))
  estframed$lower <- as.character(round(estframe$lower, dec1))
  estframed$upper <- as.character(round(estframe$upper, dec1))
  estframed$p.value <- as.character(round(estframe$p.value, 3))
  estframed$p.value[estframe$p.value<.001] <- "<0.001"
  return(list(semod=msemod, semodd=msemodd,pval=mpval, pvald=mpvald,
              estframe=estframe, estframed=estframed))
} # end function mod1f

# moderator analysis, more than one categorical or continuous variable
modmf <- function(frame1, covn="pyear,type,pstat", dec1=2, centr1="0", ref1="gen,no") {
  nframe1 <- dim(frame1)[1]
  nvar <- dim(frame1)[2]
  covnam <- unlist(strsplit(covn, split=","))
  ncov <- length(covnam)
  covnamall <- names(frame1)[-c(1:4)]
  if (sum(covnam %in% covnamall)<ncov) 
    stop("Choose a name for a variable, from variable 5 onwards")
  centrn <- as.numeric(unlist(strsplit(centr1, split=",")))
  refc <- unlist(strsplit(ref1, split=","))
  ncov <- length(covnam)
  ncentr <- length(centrn)
  classc <- rep("unknown",ncov)
  nref <- length(refc)
  for (ii in 1:ncov) {
    varnri <- (1:nvar)[names(frame1)==covnam[ii]]
    vari <- frame1[,varnri]
    classi <- class(vari)
    classc[ii] <- classi
  } # end determine classes
  if (sum(classc %in% c("numeric","integer","factor"))<ncov) 
    stop("Moderators have to be numerical or categorical")
  catmult <- rep(0, ncov) # initializing flags for factors with >2 categories
  for (ii in 1:ncov) {
    vari <- frame1[,covnam[ii]]
    if (class(vari)=="factor") {
      levi <- levels(vari)
      nlevi <- length(levi)
      if (nlevi>2) catmult[ii] <- 1
    }
  } # end set flags for factors with >2 categories
  ncatmult <- sum(catmult)
  nclassnum <- sum(classc %in% c("numeric","integer"))
  nclasscat <- sum(classc %in% c("factor"))
  if (nclassnum>0&nclassnum!=ncentr) 
    stop("The number of centering values has to be the same as the mumber of numerical variables")
  if (nclasscat!=nref)
    stop("The number of reference categories has to be the same as the number of categorival variables")
  if (nclasscat>0) {
    covnamcat <- covnam[classc=="factor"]
    for (ii in 1:nclasscat) {
      varnri <- (1:nvar)[names(frame1)==covnamcat[ii]]
      vari <- frame1[,varnri]
      vari <- relevel(vari, ref=refc[ii])
      frame1[,varnri] <- vari
    } 
  } # end relevel
  if (nclassnum>0) {
    covnamnum <- covnam[classc %in% c("numeric","integer")]
    for (ii in 1:nclassnum) {
      varnri <- (1:nvar)[names(frame1)==covnamnum[ii]]
      vari <- frame1[,varnri]
      vari <- vari-centrn[ii]
      frame1[,varnri] <- vari
    } 
  } # end centering
  formm <- as.formula(paste0("~ ", paste(covnam, collapse=" + ")))
  modm <- rma.mv(y, v, mods= formm, tdist=TRUE, data=frame1,
                 random = list(~ 1 | effectsizeID, ~ 1 | studyID))
  msemod <- sqrt(modm$sigma2)
  names(msemod) <-c("within","between")
  msemodd <- as.character(round(msemod,dec1))
  names(msemodd) <- names(msemod)
  mpval <- c(modm$QEp, modm$QMp, rep(NA,length(covnam[catmult==1])))
  names(mpval) <- c("heterogeneity","moderators",covnam[catmult==1])
  if (ncatmult>0)  {
    covnammult <- covnam[catmult==1]
    for (covnamnr in 1:length(covnammult)) {
      anovanr <- anova(modm, btt=covnammult[covnamnr])
      mpval[2+covnamnr] <- anovanr$QMp
    }
  }
  mpvald <- as.character(round(mpval,3))
  names(mpvald) <-names(mpval)
  mpvald[mpval<.001] <- "<0.001"
  estnames <- row.names(modm$b)
  estnames[1] <-"pred.ref.centr"
  estframe <- data.frame(name=estnames)
  estframe$estimate <- modm$b
  estframe$standarderror <- modm$se
  estframe$lower <- modm$ci.lb
  estframe$upper <- modm$ci.ub
  estframe$p.value <- modm$pval
  estframed <- estframe
  estframed$estimate <- as.character(round(estframe$estimate, dec1))
  estframed$standarderror <- as.character(round(estframe$standarderror, dec1))
  estframed$lower <- as.character(round(estframe$lower, dec1))
  estframed$upper <- as.character(round(estframe$upper, dec1))
  estframed$p.value <- as.character(round(estframe$p.value, 3))
  estframed$p.value[estframe$p.value<.001] <- "<0.001"
  return(list(semod=msemod, semodd=msemodd, pval=mpval, pvald=mpvald,
              estframe=estframe, estframed=estframed))
} # end function modmf

# user interface
# shiny app, user input
ui <- fluidPage(
  # App title ----
  titlePanel("Meta regression"),  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      fileInput("fil1", "Retrieve the results from the studies:",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
 
      radioButtons("sep", "Separating character for the studies:",
                   choices = c(comma = ",",
                               semicolon = ";",
                               tabulator = "\t"),
                   selected = ";"),
      
      radioButtons("desfil", "Decimal character for the studies:",
                   choices = c(comma = ",",
                               point="."),
                   selected = "."),
      
      radioButtons("studies", "Show the studies:",
                   choices = c(none="none",
                               selection="selection",
                               all="all"),
                   selected = "none"),
      
      numericInput(inputId = "nselection",
                   label = "If selection, the number of studies:",
                   value = 3),
      
      numericInput(inputId = "ndecr",
                   label = "The number of decimales in the results:",
                   value = 3),
 
      numericInput(inputId = "onemod",
                   label = "Choose a moderator:",
                   value = 5),
      
      numericInput(inputId = "centr",
                   label = "If continuous, choose a centering value:",
                   value = 0),
 
      textInput(inputId = "multmod",
                label = "Choose multiple moderators, names with comma between:",
                value = "pstat,type,pyear"),
      
      textInput(inputId = "refcat",
                label = "Choose references for the categoricals:",
                value = "no,gen"),
      
      textInput(inputId = "centrnum",
                label = "Choose centrering values for the continuous:",
                value = "0"),
      
      br(),
      br(),
      
    ),
    
    # names for the stuff to show
    mainPanel(
      h2(textOutput("basetext")),
      tableOutput("raw"),
      h3(textOutput("overalltext")),
      tableOutput("overall"),
      h3(textOutput("varcomptext")),
      tableOutput("varcomp"),
      h3(textOutput("cheungtext1")),
      h4(textOutput("cheungtext2")),
      tableOutput("cheung"),
      h3(textOutput("modtext")),
      h4(textOutput("modtextp")),
      tableOutput("modp"),
      h4(textOutput("modtextse")),
      tableOutput("modvarc"),
      h4(textOutput("modtextcoef")),
      tableOutput("modcoef"),
      h3(textOutput("multtext")),
      h4(textOutput("multmodtextp")),
      tableOutput("multmodp"),
      h4(textOutput("multmodtextse")),
      tableOutput("multmodvarc"),
      h4(textOutput("multmodtextcoef")),
      tableOutput("multmodcoef")
    )
  )
)

################################
# shiny app, server instructions

server <- function(input, output) {
  
  df <- reactive({
    read.csv(input$fil1$datapath,
             header = TRUE,
             sep = input$sep,
             quote = "",
             dec=input$desfil)
  })
  
  aw <- reactive({
    req(input$fil1)
    da <- df()
    awa <- awf(df(), dec1=input$ndecr)
    return(awa)
    })
  
  mod <- reactive({
    req(input$fil1)
    da <- df()
    return(mod1f(da,covnr=input$onemod,centr1=input$centr,dec1=input$ndecr))
  })
  
  multmod <- reactive({
    req(input$fil1)
    da <- df()
    return(modmf(da,covn=input$multmod, dec1=input$ndecr,
                 centr1=input$centrnum, ref1=input$refcat))
  })
  
  output$raw <- renderTable({
    req(input$fil1)
    da <- df()
    stud <- input$studies
    nsel <- input$nselection
    show <- NULL
    if (stud=="all") show <- da else
      if (stud=="selection") {
        id <- da[,1]
        uid <- unique(id)
        include <- sort(sample(uid,nsel))
        show <- da[id %in% include,]
      }
    return(show)
  })
  
  output$basetext <- renderText({
    req(input$fil1)
    stud <- input$studies
    res <- ""
    if (stud!="none") res <- "Data from the studies"
    return(res)
  })

  output$overalltext <- renderText({
    req(input$fil1)
    return("Overall meta regression")
  })

  output$varcomptext <- renderText({
    req(input$fil1)
    return("Variance components with likelihood ratio tests")
  })
  
  output$cheungtext1 <- renderText({
    req(input$fil1)
    return("Total variance, percent distribution by levels")
  })
  
  output$cheungtext2 <- renderText({
    req(input$fil1)
    return("(Cheung's procedure)")
  })
  
  output$modtext <- renderText({
    req(input$fil1)
    return("Analysis for one moderator")
  })
  
  output$modtextp <- renderText({
    req(input$fil1)
    return("Summary p-values")
  })  
  
  output$modtekstse <- renderText({
    req(input$fil1)
    return("Variance components, standard deviation")
  })
  
  output$multtext <- renderText({
    req(input$fil1)
    return("Model with multiple moderators")
  })
  
  output$multmodtextp <- renderText({
    req(input$fil1)
    return("Summary p-values")
  })
  
  output$multmodtextse <- renderText({
    req(input$fil1)
    return("Variance components, standard deviation")
  })
  
  output$multmodtextcoef <- renderText({
    req(input$fil1)
    return("Coefficients")
  })
  
  output$overall <- renderTable({
    req(input$fil1)
    awa <- aw()
    awo <- awa$oeffd
    awot <- t(awo)
    return(awot)    
  })

  output$varcomp <- renderTable({
    req(input$fil1)
    awa <- aw()
    awv <- awa$varcompd
    names(awv)[3] <- "p-value"
    return(awv)    
  })
  
  output$cheung <- renderTable({
    req(input$fil1)
    awa <- aw()
    awc <- awa$cheungd
    awct <- t(awc)
    return(awct)    
  })
  
  output$modp <- renderTable({
    req(input$fil1)
    md1 <- mod()
    pval <- md1$pvald
    pvalt <- t(pval)
    return(pvalt)
  }) 
  
  output$modvark <- renderTable({
    req(input$fil1)
    md1 <- mod()
    sem <- md1$semodd
    semt <- t(sem)
    return(semt)
  })
  
  output$modcoef <- renderTable({
    req(input$fil1)
    md1 <- mod()
    return(md1$estframed)
  })
  
  output$multmodp <- renderTable({
    req(input$fil1)
    mu1 <- multmod()
    pvalu <- mu1$pvald
    pvalut <- t(pvalu)
    return(pvalut)
  })
  
  output$multmodvarc <- renderTable({
    req(input$fil1)
    mu1 <- multmod()
    semu <- mu1$semodd
    semut <- t(semu)
    return(semut)
  })
 
  output$multmodcoef <- renderTable({
    req(input$fil1)
    mu1 <- multmod()
    return(mu1$estframed)
  })

}

#############################
# shiny app, putting together
shinyApp(ui, server)
