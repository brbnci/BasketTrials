
########################################################################
#								                                                       #
# Main R functions for finding the posterior probabilities at interim  #
# analysis for basket trials                                           #
#                                                                      #
########################################################################


################ Functions ###############
####
#### From Dr. Simon, Feb 2014
#### Version 1, Feb 2014
#### Modified for unequal prevalences (Jyothi, Mar. 2014)
#### Version 3, Apr 2014: Approved by Rich. Hosted in glimmer
#### Version 3, May 2014: Moved to shinyapps.io
### Version 4, Oct 2015: From version 3, changed the 'Introduction' page a bit (i.e., documentation) 
###               in ui.R. 
###               Updated to shinyIncubator ver 0.2.2 and shiny version 0.12.2 (latest versions)
###               Column headers now fail to display in matrixInput. Sent a question to google groups.
###               Changed label text and strata names for now. 
###               Added a tooltip to input 'delta' and made it into a slider. 
###               Probability values move at intervals of 0.01 (gamma, lambda, p_lo, p_hi)
###               8-Oct-2015: Modified the matrixInput function from tableinput.R 
###               in the shinyIncubator 0.2.2 package to matrixInputNew. See below.
### Version 5, Sep 2016: Corrected true positive and negative calculations. Sensitivity and specificity 
###               are just the true positive and true negative rate respectively.
###               Also calculate and output the expected number of indeterminate strata
### Version 6, Oct 2016: Used the new fast formulae for calculation pf posterior probabilities 
###               (after verification and comparing with results of existing approach), increased max 
###               number of strata to 12, and improving formatting of post prob output table.
##########################################

number2binary <- function(number,noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)){return(binary_vector)}
  else{binary_vector[-(1:(length(binary_vector)-noBits))]}
}


find_postk <- function(lambda, gamma, r, n, plo, phi){
  ## P0
  P0.num <- (1-gamma)*prod(dbinom(r, n, plo))
  P0.deno <- gamma*prod(dbinom(r, n, phi))
  P0 <- 1/(1 + (P0.num/P0.deno))
  
  ## Pk_ind
  Pk.num <- gamma*dbinom(r, n, phi)
  Pk.deno <- gamma*dbinom(r, n, phi) + (1-gamma)*dbinom(r, n, plo)
  Pk <- Pk.num/Pk.deno
  
  ## P[pi = 1|data]
  pi1.num <- (1-lambda)*prod(gamma*dbinom(r, n, phi) + (1-gamma)*dbinom(r, n, plo))
  pi1.deno <- lambda*(gamma*prod(dbinom(r, n, phi)) + (1-gamma)*prod(dbinom(r, n, plo)))
  pi1 <- 1/(1+(pi1.num/pi1.deno))
  
  ## Finally, post prob of activity in each strata
  postk <- P0*pi1 + Pk*(1 - pi1)
  postk
}


samplesize <- function(lambda, gamma, plo, phi, kk, nrep, ntot, nblock, prev=rep(1,kk), delta) {
  pnull <- 1 - gamma
  errate <- matrix(0,nrow=kk+1,ncol=7)
  
  for (npos in 0:kk) {  ### for each number of active strata: 0, 1, 2, ... or all kk active
    truepos <- falsepos <- falseneg <- trueneg <- 0
    ess <- 0	### expected sample size
    fpAtAllNeg <- 0	### To calculate P[no FPs | all strata -ve]
    expDiscovery <- 0	### To calculate expected no. of discoveries
    indeter <- 0  ### To calculate expected no. of indeterminate strata
    #print(npos)     
    for (irep in seq(nrep)) {	### simulate for number of 'nrep' reps
      pval <- plo + 0*seq(kk)	### all inactive strata
      if (npos > 0) {
        ind <- seq(npos)
        pval[ind]<- pval[ind]+(phi-plo)	### 1:npos active strata
      }
      kknow <-  kk
      activevec <- activevecnew <- 1+0*seq(kk)  ### All strata open initially
      finalprob <- rep(NA,kk)  ## Post prob for all strata initialized to missing 
      ncum <- rcum <-  0*seq(kk)	### evaluables and responders initialized to zero
      prevnew <- prev
      for (npat in seq(nblock,ntot,nblock)) {	### each interim analysis after 5 evaluable subjects
        n <- rmultinom(1, nblock, prevnew)      ### simulate number of evaluables in each strata    
        ncum <- ncum + n        
        r <- rbinom(kk, n, pval) 	### simulate number of responses in each strata (shd it be kknow?) (kk also ok..)
        r[activevec==0] <- 0
        rcum <- rcum + r     
        val <- find_postk(lambda, gamma, rcum, ncum, plo, phi)    ###(shd it be kknow?) (kknow gives NAs results)
        ### Thats bcoz, the prior is still set as if there are kk strata
        ### I think kk is ok, bcoz only accrual to the strata is stopped, the parameters set initially are still the same
        activevecnew <- activevec * (val < delta) * (val > 1-delta)	### strata getting closed (> d or < 1-d) will become 0
        ind <- (activevec == 1 & activevecnew == 0)	### strata getting closed in this iteration
        finalprob[ind] <- val[ind]
        #print(finalprob)
        kknow <- sum(activevecnew)	### current number of active strata
        activevec <- activevecnew 	### new activevec
        prevnew[activevec==0] <- 0
        if (kknow == 0) break
     }                               
      if (npos > 0) {
        truepos <- truepos + sum(finalprob[1:npos] > delta, na.rm=TRUE)
        falseneg <- falseneg + sum(finalprob[1:npos] < 1-delta, na.rm=TRUE)
      }
      if (npos < kk) {
        falsepos <- falsepos + sum(finalprob[(npos+1):kk] > delta, na.rm=TRUE)
        trueneg <- trueneg + sum(finalprob[(npos+1):kk] < 1-delta, na.rm=TRUE)
      }
      ess <- ess + sum(ncum) ### cexpected sample size
      indeter <- indeter + sum(is.na(finalprob))   ## the posterior prob for indeterminate strata will remain NA
      if (npos == 0) {	### To calculate P[no FPs | all strata -ve]
        if (sum(finalprob > delta, na.rm=TRUE) > 0)
          fpAtAllNeg <- fpAtAllNeg+1
      }
    }	### nrep simulations

    if (npos == 0) {probNoFPallNeg <- (1-fpAtAllNeg/nrep)}	### P[no FPs | all strata -ve]
    expDiscovery <- (truepos+falsepos)/nrep
    truepos <- truepos/nrep
    falsepos <- falsepos/nrep
    falseneg <- falseneg/nrep
    trueneg <- trueneg/nrep
    ess <- ess/nrep
    indeter <- indeter/nrep
    
    if(npos > 0){truepos <- truepos/npos; falseneg <- falseneg/npos}
    if(npos < kk){falsepos <- falsepos/(kk-npos); trueneg <- trueneg/(kk-npos)}
    
    errate[npos+1,1] <- truepos
    errate[npos+1,2] <- falsepos
    errate[npos+1,3] <- falseneg
    errate[npos+1,4] <- trueneg
    errate[npos+1,5] <- ess
    errate[npos+1,6] <- expDiscovery
    errate[npos+1,7] <- indeter
  }   ### over for nos. of active strata = 0, 1, 2,...kk
  #print(errate)
  ### Multiply errate with prior prob of active strata being 0, 1, 2,...kk
  prob <- (1-lambda)*dbinom(seq(0,kk,1),kk,pnull)
  prob[1] <- lambda*(1-gamma) + prob[1]
  prob[kk+1] <- lambda*gamma + prob[kk+1]

  valpos <- 0:kk
  valneg <- kk - valpos
  
  falseposav <- sum(prob*errate[,2]*valneg)/sum(prob*valneg)
  trueposav <-  sum(prob*errate[,1]*valpos)/sum(prob*valpos)
  falsenegav <-  sum(prob*errate[,3]*valpos)/sum(prob*valpos)
  truenegav <-   sum(prob*errate[,4]*valneg)/sum(prob*valneg)
  essav <- sum(prob*errate[,5])
  expDiscovery <- sum(prob*errate[,6])
  indeterav <- sum(prob*errate[,7])

  cat(trueposav,falsenegav,falseposav,truenegav,essav,expDiscovery,indeterav,probNoFPallNeg,"\n")
  list(essav=essav, expDiscovery=expDiscovery, indeterav=indeterav, probNoFPallNeg=probNoFPallNeg,
       trueposav=trueposav,falsenegav=falsenegav,falseposav=falseposav,truenegav=truenegav)
}

### Jyothi, 8-Oct-2015: Modified the matrixInput function from tableinput.R in the shinyIncubator 0.2.2 package.
### Copied from https://github.com/rstudio/shiny-incubator/blob/master/R/tableinput.R
### Lines pertining to +/- buttons and table headings hide are commented off and nothing is deleted
### These changes make it not possible for the user to increase the number of strata rows in the table without 
### actually increasing 'kk'. Table column headings are also now visible. Tested now with shiny 0.12.2.
### We still need to include the shinyincubator package for the Javascript functions

matrixInputNew <- function(inputId, label, data) {
  addResourcePath(
    prefix='tableinput', 
    directoryPath=system.file('tableinput', 
                              package='shinyIncubator'))
  
  tagList(
    singleton(
      tags$head(
        tags$link(rel = 'stylesheet',
                  type = 'text/css',
                  href = 'tableinput/tableinput.css'),
        tags$script(src = 'tableinput/tableinput.js')
      )
    ),
    
    tags$div(
      class = 'control-group tableinput-container',
      tags$label(
        class = "control-label",
        label#,
        #        tags$div(
        #          class = 'tableinput-buttons',
        #          tags$button(
        #            type = 'button', class = 'btn btn-mini tableinput-settings hide',
        #            tags$i(class = 'glyphicon glyphicon-cog icon-cog')
        #          ),
        #          HTML('<a href="#" class="tableinput-plusrow"><i class="glyphicon glyphicon-plus-sign icon-plus-sign"></i></a>'),
        #          HTML('<a href="#" class="tableinput-minusrow"><i class="glyphicon glyphicon-minus-sign icon-minus-sign"></i></a>')
        #        )
      ),
      tags$table(
        id = inputId,
        class = 'tableinput data table table-bordered table-condensed',
        tags$colgroup(
          lapply(names(data), function(name) {
            tags$col('data-name' = name,
                     'data-field' = name,
                     'data-type' = 'numeric')
          })
        ),
        tags$thead(
          #          class = 'hide',
          tags$tr(
            lapply(names(data), function(name) {
              tags$th(name)
            })
          )
        ),
        tags$tbody(
          lapply(1:nrow(data), function(i) {
            tags$tr(
              lapply(names(data), function(name) {
                tags$td(
                  div(tabindex=0, as.character(data[i,name]))
                )
              })
            )
          })
        )
      ),
      tags$div(
        class = 'tableinput-editor modal hide fade',
        tags$div(
          class = 'modal-header',
          HTML('<button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button>'),
          tags$h3(label)
        ),
        tags$div(
          class = 'modal-body',
          
          HTML('
               <form class="form-horizontal">
               <div class="control-group">
               <label class="control-label">Rows</label>
               <div class="controls">
               <input type="number" class="tableinput-rowcount">
               </div>
               </div>
               <div class="control-group">
               <label class="control-label">Columns</label>
               <div class="controls">
               <input type="number" class="tableinput-colcount">
               </div>
               </div>
               </form>'
          )
          ),
        tags$div(
          class = 'modal-footer',
          tags$a(href = '#', class = 'btn btn-primary tableinput-edit', 'OK'),
          tags$a(href = '#',
                 class = 'btn',
                 'data-dismiss' = 'modal',
                 'aria-hidden' = 'true',
                 'Cancel')
        )
          )
      )
    )
}
