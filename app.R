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
  names(res) <- c("within","intermediate","between")
  return(list(rma=rmaobj, varcomp=res))
} # end function cheungdek

awf <- function(frame1, dec1=2) {
  overall <- rma.mv(y, v, random = list(~ 1 | effectsizeID, ~ 1 | studyID), 
                    tdist=TRUE, data=frame1)
  megger <- rma.mv(y, v, mods=~I(sqrt(v)), random = list(~ 1 | effectsizeID, ~ 1 | studyID), 
                   tdist=TRUE, data=frame1) # modified Egger test, Marengo and Montag
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
  ohet <- c(overall$QE,overall$k - overall$p,overall$QEp)
  # changed degrees of freedom k-p, to work with metafor 3.0_0 as well, 
  # where dfs is a character input. k-p is taken from the code for rma.mv
  names(ohet) <- c("Q","df","p")
  ohetd <- as.character(round(ohet,dec1))
  names(ohetd) <- names(ohet)
  ohetd[2] <- as.character(round(ohet[2]))
  ohetd[3] <- as.character(round(ohet[3],3))
  if (ohet[3]<.001) ohetd[3] <- "<.001"
  if (ohet[3]<.001) ohetd[3] <- "<.001"
  mhet <- c(megger$QE,megger$k - megger$p,megger$QEp)
  names(mhet) <- c("Q","df","p")
  mhetd <- as.character(round(mhet,dec1))
  names(mhetd) <- names(mhet)
  mhetd[2] <- as.character(round(mhet[2]))
  mhetd[3] <- as.character(round(mhet[3],3))
  if (mhet[3]<.001) mhetd[3] <- "<.001"
  megger.res <- c(megger$b[2], megger$se[2], megger$zval[2], megger$k - megger$p,
                  megger$ci.lb[2], megger$ci.ub[2], megger$pval[2])
  names(megger.res) <- c("estimate","standarderror","t","df","lower","upper","p-value")
  meggerd <- as.character(megger.res)
  names(meggerd) <-names(megger.res)
  meggerd[c(1:3,5,6)] <- as.character(round(megger.res[c(1:3,5,6)],dec1))
  meggerd[4] <- as.character(round(megger.res[4]))
  meggerd[7] <- as.character(round(megger.res[7],3))
  if (megger.res[7]<.001) meggerd[7] <- "<.001"
  varcomp <- data.frame(standarderror=seo, plikrat=c(pwithin,pbetween),
                        row.names=c("within","between"))
  varcompd <- data.frame(domain=c("within","between"), 
                         standarderror=varcomp$standarderror, 
                         plikrat=varcomp$plikrat)
  varcompd$standarderror <- round(varcomp$standarderror,dec1)
  varcompd$plikrat <- round(varcomp$plikrat,3)
  varcompd$standarderror <- as.character(varcompd$standarderror) 
  varcompd$plikrat <- as.character(varcompd$plikrat)
  if (varcomp[1,2]<.001) varcompd[1,3] <- "<.001" 
  if (varcomp[2,2]<.001) varcompd[2,3] <- "<.001" 
  acheung <- cheungdek(frame1, overall)$varcomp
  names(acheung) <- c("within","intermediate","between")
  acheungd <- as.character(round(acheung,dec1))
  names(acheungd) <- names(acheung)
  return(list(ovall=overall, oeff=oeff, oeffd=oeffd, varcomp=varcomp,
              varcompd=varcompd, cheung=acheung, cheungd=acheungd,
              megger=megger.res, meggerd=meggerd,
              ohet=ohet,ohetd=ohetd,mhet=mhet,mhetd=mhetd))
} # end function awf

# function to extract correlations between levels, one moderator
corr1catf <- function(vb, levcat) { 
  # vb, variance-covariance matrix of the moderator 
  # levcat, levels of the moderator, in correct order
  ncat <- length(levcat)
  vbcat <- vb # initializing, should be var-cov matr for means not contrasts
  rownames(vbcat) <- levcat
  colnames(vbcat) <- levcat
  vbcat[1,-1] <- vb[1,-1] + vb[1,1] # first row correct 
  vbcat[-1,1] <- vb[-1,1] + vb[1,1] # first column correct 
  vbcat[-1,-1] <- vb[-1,-1] + vb[1,1] + # what remails correct
    matrix(rep(vb[1,-1],ncat-1), ncol=ncat-1, byrow=TRUE) +
    matrix(rep(vb[-1,1],ncat-1), nrow=ncat-1)
  vbcor <- diag(1/sqrt(diag(vbcat))) %*%
    vbcat %*% diag(1/sqrt(diag(vbcat)))
  diag(vbcor) <- sqrt(diag(vbcat)) # replace diag by se for easy check
  rownames(vbcor) <- levcat
  colnames(vbcor) <- levcat
  return(list(varcov=vb, varlev=vbcat, corr=vbcor))
} # end function corr1catf

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
      if (catnr==1) {
        msemod <- sqrt(modc$sigma2)
        modcorr <- corr1catf(vb=modc$vb, levcat=catc)$corr
        modcorrd <- matrix(as.character(round(modcorr, dec1)), ncol=ncatc)
        modcorrd <- data.frame(cbind(catc, modcorrd), row.names=catc)
        names(modcorrd) <- c("category",catc)
      }
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
    modcorr <- "no correlations for a continuous moderator"
    modcorrd <- "no correlations for a continuous moderator"
  } # end computations for continuous variable
  names(msemod) <-c("within","between")
  msemodd <- as.character(round(msemod,dec1))
  names(msemodd) <- names(msemod)
  hetonlyp <- 1
  if (classc=="factor") if (ncatc>2) hetonlyp <- 0
  if (hetonlyp==1) {
    hetm1 <- c(modr$QE,modr$k - modr$p,modr$QEp)
    names(hetm1) <- c("Q","df","p")
    hetm1d <- as.character(hetm1)
    names(hetm1d) <- names(hetm1)
    hetm1d[1] <- as.character(round(hetm1[1],dec1))
    hetm1d[2] <- as.character(round(hetm1[2]))
    hetm1d[3] <- as.character(round(hetm1[3],3))
    hetm1d[3][hetm1[3]<.001] <- "<0.001"
  } else {
    hetm1 <- c(modr$QE,modr$k - modr$p,modr$QEp, modr$QM, modr$QMdf, modr$QMp)
    names(hetm1) <- c("Q","df","p","Qm","df1m","df2m","pm")
    hetm1d <- hetm1
    names(hetm1d) <- names(hetm1)
    hetm1d <- as.character(hetm1d)
    hetm1d[1] <- as.character(round(hetm1[1],dec1))
    hetm1d[2] <- as.character(round(hetm1[2]))
    hetm1d[3] <- as.character(round(hetm1[3],3))
    hetm1d[3][hetm1[3]<.001] <- "<0.001"
    hetm1d[4] <- as.character(round(hetm1[4],dec1))
    hetm1d[5] <- as.character(round(hetm1[5]))
    hetm1d[6] <- as.character(round(hetm1[6]))
    hetm1d[7] <- as.character(round(hetm1[7],3))
    hetm1d[7][hetm1[7]<.001] <- "<0.001"   
  } # end p-values for heterogeneity and the moderator
  names(estframe)[1] <- names1
  estframed <- estframe
  estframed$estimate <- as.character(round(estframe$estimate, dec1))
  estframed$standarderror <- as.character(round(estframe$standarderror, dec1))
  estframed$lower <- as.character(round(estframe$lower, dec1))
  estframed$upper <- as.character(round(estframe$upper, dec1))
  estframed$p.value <- as.character(round(estframe$p.value, 3))
  estframed$p.value[estframe$p.value<.001] <- "<0.001"
  return(list(semod=msemod, semodd=msemodd, 
              estframe=estframe, estframed=estframed,
              hetm1=hetm1d, corrlevel=modcorr, corrleveld=modcorrd))
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
  hetmult <- c(modm$QE, modm$k - modm$p, modm$QEp, 
               modm$QM, modm$QMdf, modm$QMp)
  names(hetmult) <- c("Q","df","p","Qm","df1m","df2m","pm")
  hetmultd <- hetmult
  names(hetmultd) <- names(hetmult)
  hetmultd <- as.character(hetmultd)
  hetmultd[1] <- as.character(round(hetmult[1],dec1))
  hetmultd[2] <- as.character(round(hetmult[2]))
  hetmultd[3] <- as.character(round(hetmult[3],3))
  hetmultd[3][hetmult[3]<.001] <- "<0.001"
  hetmultd[4] <- as.character(round(hetmult[4],dec1))
  hetmultd[5] <- as.character(round(hetmult[5]))
  hetmultd[6] <- as.character(round(hetmult[6]))
  hetmultd[7] <- as.character(round(hetmult[7],3))
  hetmultd[7][hetmult[7]<.001] <- "<0.001"   
  mpval <- rep(NA,length(covnam[catmult==1]))
  names(mpval) <- covnam[catmult==1]
  if (ncatmult>0)  {
    covnammult <- covnam[catmult==1]
    for (covnamnr in 1:length(covnammult)) {
      anovanr <- anova(modm, btt=covnammult[covnamnr])
      mpval[covnamnr] <- anovanr$QMp
    }
  }
  mpvald <- as.character(round(mpval,3))
  names(mpvald) <- names(mpval)
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
  return(list(semod=msemod, semodd=msemodd, mpval=mpval, mpvald=mpvald,
              hetmult=hetmultd, estframe=estframe, estframed=estframed))
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
      radioButtons("ls", "Language for results:",
                   choices = c(English = "English",
                               norsk = "norsk"),
                   selected = "English"),
      
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
                   label = "The number of decimals in the results:",
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
      h3(textOutput("overallhettext")),
      tableOutput("overallhet"),
      h3(textOutput("varcomptext")),
      tableOutput("varcomp"),
      h3(textOutput("cheungtext1")),
      h4(textOutput("cheungtext2")),
      tableOutput("cheung"),
      h4(textOutput("meggertext")),
      tableOutput("megger"),
      h4(textOutput("meggerhettext")),
      tableOutput("meggerhet"),
      plotOutput(outputId = "funnel"),
      h3(textOutput("modtext")),
      h4(textOutput("modtexttest")),
      h5(textOutput("modtexttest2")),
      tableOutput("modtest"),
      h4(textOutput("modtextp")),
      tableOutput("modp"),
      h4(textOutput("modtextse")),
      tableOutput("modvarc"),
      h4(textOutput("modtextcoef")),
      tableOutput("modcoef"),
      h4(textOutput("modtextcorr")),
      tableOutput("modcorr"),
      h3(textOutput("multtext")),
      h4(textOutput("multmodtexthet")),
      h5(textOutput("multmodtexthet2")),
      tableOutput("multmodhet"),
      h4(textOutput("multmodtextp")),
      h5(textOutput("multmodtextp2")),
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
             dec=input$desfil,
             stringsAsFactors=TRUE)
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
    lasp <- input$ls
    res <- ""
    if (stud!="none"&lasp=="English") res <- "Data from the studies"
    if (stud!="none"&lasp=="norsk") res <- "Grunnlagsdata"
    return(res)
  })

  output$overalltext <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Overall meta regression"
    if (lasp=="norsk") res <- "Overall meta-regresjon"
    return(res)
  })
  
  output$overallhettext <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Overall heterogeneity"
    if (lasp=="norsk") res <- "Overall heterogenitet"
    return(res)
  })

  output$varcomptext <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Variance components with likelihood ratio tests"
    if (lasp=="norsk") res <- "Varianskomponenter med likelihood ratio-tester"
    return(res)
  })
  
  output$cheungtext1 <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Total variance, percent distribution by levels"
    if (lasp=="norsk") res <- "Total varians, prosentfordeling på nivåene"
    return(res)
  })
  
  output$cheungtext2 <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "(Cheung's procedure)"
    if (lasp=="norsk") res <- "(Cheungs framgangsmåte)"
    return(res)
  })
  
  output$meggertext <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Modified Egger test (Marengo and Montag)"
    if (lasp=="norsk") res <- "Modifisert Egger-test (Marengo og Montag)"
    return(res)
  })
  
  output$meggerhettext <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Modified Egger, residual heterogeneity"
    if (lasp=="norsk") res <- "Modifisert Egger, residual-heterogenitet"
    return(res)
  })
  
  output$modtext <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Analysis for one moderator"
    if (lasp=="norsk") res <- "Analyse for en moderator"
    return(res)
  })
  
  output$modtexttest <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Summary tests. First residual heterogeneity"
    if (lasp=="norsk") res <- "Samlede tester, først residual-heterogenitet"
    return(res)
  })  
  
  output$modtexttest2 <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- ".....next (suffix m) for the moderator if necessary"
    if (lasp=="norsk") res <- ".....deretter (endelse m) for moderatoren om nødvendig"
    return(res)
  })  
  
  output$modtekstse <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Variance components, standard deviation"
    if (lasp=="norsk") res <- "Varianskomponenter, standardavvik"
    return(res)
  })
  
  output$multtext <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Model with multiple moderators"
    if (lasp=="norsk") res <- "Modell med flere moderatorer"
    return(res)
  })
  
  output$multmodtextp <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Summary p-values for individual moderators"
    if (lasp=="norsk") res <- "Samlede p-verdier for de enkelte moderatorene"
    return(res)
  })
  
  output$multmodtextp2 <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "... if categorical with more than two caregories"
    if (lasp=="norsk") res <- "... for de kategoriske med mer enn to kategorier"
    return(res)
  })
  
  
  output$multmodtexthet <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Tests of residual heterogeneity"
    if (lasp=="norsk") res <- "tester for residual-hogenitet"
    return(res)
  })
  
  output$multmodtexthet2 <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "... and (suffix m) of moderators together"
    if (lasp=="norsk") res <- "... og (endelse m) for moderatorene samlet"
    return(res)
  })
  
  output$multmodtextp2 <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "... categorical, more than two categories"
    if (lasp=="norsk") res <- "... kategoriske, med mer enn to kategorier"
    return(res)
  })
  
  output$multmodtextse <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Variance components, standard deviation"
    if (lasp=="norsk") res <- "Varianskomponenter, standardfeil"
    return(res)
  })
  
  output$multmodtextcoef <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Coefficients"
    if (lasp=="norsk") res <- "Koeffisienter"
    return(res)
  })
  
  output$modtextcorr <- renderText({
    req(input$fil1)
    lasp <- input$ls
    res <- ""
    if (lasp=="English") res <- "Correlations (diagonal, standard errors)"
    if (lasp=="norsk") res <- "Korrelasjoner (standardfeil langs diagonalen)"
    return(res)
  })
  
  output$overall <- renderTable({
    req(input$fil1)
    lasp <- input$ls
    awa <- aw()
    awo <- awa$oeffd
    if (lasp=="norsk") names(awo) <- c("effekt","nedre","øvre","p-verdi")
    awot <- t(awo)
    return(awot)    
  })
  
  output$overallhet <- renderTable({
    req(input$fil1)
    awa <- aw()
    heto <- awa$ohetd
    hetot <- t(heto)
    return(hetot)    
  }) 

  output$varcomp <- renderTable({
    req(input$fil1)
    lasp <- input$ls
    awa <- aw()
    awv <- awa$varcompd
    if (lasp=="English") names(awv)[3] <- "p-value"
    if (lasp=="norsk") names(awv) <- c("område","standardfeil","p-verdi")
    if (lasp=="norsk") awv[,1] <- c("innen","mellom")
    return(awv)    
  })
  
  output$cheung <- renderTable({
    req(input$fil1)
    lasp <- input$ls
    awa <- aw()
    awc <- awa$cheungd
    if (lasp=="norsk") names(awc) <- c("inni","midt","mellom")
    awct <- t(awc)
    return(awct)    
  })
  
  output$megger <- renderTable({
    req(input$fil1)
    lasp <- input$ls
    awa <- aw()
    awm <- awa$meggerd
    if (lasp=="norsk") names(awm) <- 
      c("estimat","standardfeil","t","df","nedre","øvre","p-verdi")
    awmt <- t(awm)
    return(awmt)    
  }) 
  
  output$meggerhet <- renderTable({
    req(input$fil1)
    awa <- aw()
    mego <- awa$mhetd
    megot <- t(mego)
    return(megot)    
  })  
  
  output$funnel <- renderPlot({
    req(input$fil1)
    awa <- aw()
    ova <- awa$ovall
    return(funnel(ova, main="Funnel diagram"))
  })
  
  output$modtest <- renderTable({
    req(input$fil1)
    md1 <- mod()
    testval <- md1$hetm1
    len <- length(testval)
    if (len==7) names(testval) <-
      c("Q","df","p","Qm","df1m","df2m","pm")
    testvalt <- t(testval)
    return(testvalt)
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
    lasp <- input$ls
    md1 <- mod()
    md1e <- md1$estframed
    if (lasp=="norsk") 
      names(md1e)[-1] <- c("estimat","standardfeil","nedre","øvre","p-verdi")
    return(md1e)
  })
  
  output$modcorr <- renderTable({
    req(input$fil1)
    lasp <- input$ls
    md1 <- mod()
    md1c <- md1$corrleveld
    return(md1c)
  })
  
  output$multmodp <- renderTable({
    req(input$fil1)
    mu1 <- multmod()
    mpvalu <- mu1$mpvald
    mpvalut <- t(mpvalu)
    return(mpvalut)
  })
  
  output$multmodhet <- renderTable({
    req(input$fil1)
    mu1 <- multmod()
    mshet <- mu1$hetmult
    names(mshet) <- c("Q","df","p","Qm","df1m","df2m","pm")
    mshett <- t(mshet)
    return(mshett)
  })
  
  output$multmodvarc <- renderTable({
    req(input$fil1)
    lasp <- input$ls
    mu1 <- multmod()
    semu <- mu1$semodd
    if (lasp=="norsk") names(semu) <- c("inni","mellom")
    semut <- t(semu)
    return(semut)
  })
 
  output$multmodcoef <- renderTable({
    req(input$fil1)
    lasp <- input$ls
    mu1 <- multmod()
    mefd <- mu1$estframed
    if (lasp=="norsk") names(mefd) <- 
      c("navn","estimat","standardfeil","nedre","øvre","p-verdi")
    return(mefd)
  })

}

#############################
# shiny app, putting together
shinyApp(ui, server)
