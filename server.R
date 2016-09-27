###### App for Basket trials
### Version 1, Feb 2014
### Version 3, Apr 2014: Approved by Rich. Hosted in glimmer
### Version 3, May 2014: Moved to shinyapps.io
### Version 4, Oct 2015: From version 3, changed the 'Introduction' page a bit (i.e., documentation) 
###               in ui.R. 
###               Updated to shinyIncubator ver 0.2.2 and shiny version 0.12.2 (latest versions)
###               Column headers now fail to display in matrixInput. Sent a question to google groups.
###               Changed label text and strata names for now. 
###               Added a tooltip to input 'delta' and made it into a slider. 
###               Probability values move at intervals of 0.01 (gamma, lambda, p_lo, p_hi).
###               8-Oct-2015: Modified the matrixInput function from tableinput.R 
###               in the shinyIncubator 0.2.2 package
######

###### Load required libraries

library(shiny)
library(shinyIncubator)
options(warn = -1) ### supress warning messages
#options(warn = 2) ### display warnings as errors for debugging

shinyServer(function(input, output) {
  ##### Slider for plo
  output$plo_slider <- renderUI({
    sliderInput("plo", p("Threshold probability for inactivity: p_lo"), 
                min=0, max=input$phi, value=0.05, step=0.01)
  })
  
  #### Table to enter interim data
  output$interim <- renderUI({
    kk <- input$kk
    val <- data.frame(Strata=paste("Strata",c(1:kk),sep="_"), Number_Evaluable="0", Number_Responses="0")
    matrixInputNew(inputId="interimdata", label=h5("Enter number evaluable and number of responders at interim analysis"), data=val)
  })  
  
  #### Table to enter prevalence
  output$prevalence <- renderUI({
    kk <- input$kk
    val <- data.frame(Strata=paste("Strata",c(1:kk),sep="_"), Prevalence=0)
    matrixInputNew(inputId="prev", label=h5("Enter the prevalence (proportion) in each strata:"), data=val)
  }) 
  
  ##### Computation, results, output
  
  ### Calculate posterior prob
  postk <- reactive({
    kk <- input$kk
    nr <- input$interimdata
    if(is.null(nr) || nrow(nr) != kk) {
      return(NULL)
    } else {
    lambda <- input$lambda
    gamma <- input$gamma
    plo <- input$plo
    phi <- input$phi
    n <- nr[,2]
    r <- nr[,3]
 
    postprob <- find_postk(lambda, gamma, r, n, plo, phi, kk)
    rown <- sapply(1:kk, function(i) {
      paste0("Strata_", i)
    })
    data.frame(Number_Evaluable=n, Number_Responses=r, Posterior_Prob_Activity=postprob, row.names=rown) 
    }
  })
  
  ### Display posterior probabilities as a table
  output$displayPostprob <- renderTable({
    res <- postk()
    data.frame(res)
  },digits=c(0,0,0,3))
  
  
  ### Simulations for sample size planning
  ssplan <- reactive({
    lambda <- input$lambda
    gamma <- input$gamma
    plo <- input$plo
    phi <- input$phi
    kk <- input$kk
    nrep <- input$nrep
    ifelse(input$equal, prev <- rep(1,kk), prev <- input$prev[,2]) 
    #print(prev)
    if(!input$equal) {
      if(sum(prev) != 1) (stop(paste("Entered prevalences do not add to 1.","\n", 
                                     "Entered values to be > 0 and < 1 and sum to 1.",sep="")))
      if(any(prev <= 0)) (stop(paste("Entered prevalences contains non-positive values.","\n", 
                               "Entered values to be > 0 and < 1 and sum to 1.",sep="")))
    }
    delta <- input$delta
    maxs <- input$maxs
    nblock <- input$nblock
    if (maxs < nblock) stop("Max sample size should be greater than subjects for interim analysis")
    #ns <- seq(mins,maxs,10)
    
    res <- samplesize(lambda, gamma, plo, phi, kk, nrep, maxs, nblock, prev, delta)
    res
  })
  

  ### Display text output from simulation
  output$simResults <- renderUI({
    if (input$sstart == 0)
      return()
    isolate({simres <- ssplan()})
    div(span(h5("Simulation Results:"),style="color:darkred"),
        p("Expected sample size:",round(simres$essav,2)),
        p("True positive rate (Expected sensitivity):",round(simres$trueposav,2)),
        p("False negative rate:",round(simres$falsenegav,2)),
        p("True negative rate (Expected specificity):",round(simres$truenegav,2)),
        p("False positive rate:",round(simres$falseposav,2)),
        #p("Expected number of discoveries:",round(simres$expDiscovery,2)),   ## no need as per Rich, Sep-26-2016
        p("Expected number of indeterminate strata:",round(simres$indeterav,2)),
        p("Probability of no false positives when all strata are negative:", round(simres$probNoFPallNeg,2)),
        div(HTML('<HR SIZE=5 WIDTH=\"100%\" NOSHADE>'),
            p(span(em("Analyzed on:",date()),style="color:blue"))
            )
      )
   })
  
  ### Display date and time of analysis
  output$dtime <- renderUI({
    div(HTML('<HR SIZE=\"6\" WIDTH=\"100%\" NOSHADE>'),
    p(span(em("Analyzed on:",date()),style="color:blue"))
    )
  })
})
