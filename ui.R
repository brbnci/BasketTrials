library(shiny)

### Define UI for Basket Trial application
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
### Version 5, Sep 2016: Corrected true positive and negative calculations. Sensitivity and specificity 
###               are just the true positive and true negative rate respectively.
###               Also calculate and output the expected number of indeterminate strata
### Version 6, Oct 2016: Used the new fast formulae for calculation pf posterior probabilities 
###               (after verification and comparing with results of existing approach), increased max 
###               number of strata to 12, and improving formatting of post prob output table.
######

shinyUI(pageWithSidebar(
  headerPanel(list(HTML('<A href="http://www.cancer.gov" target="_blank"><img src="nci_new.PNG", style="float:center"></A>'),
                   HTML('<A href="http://brb.nci.nih.gov" target="_blank"><img src="brp.PNG", style="float:center"></A>'),
                      #span(h1("Biometric Research Branch", align = "center"),style="color:#00aaff"),
                      span(h2("Design, Monitoring and Analysis of Bayesian Basket Discovery Trials", align="center"))),
              "Basket Trials"),
  sidebarPanel(
    wellPanel(
      selectInput("Navigation", strong("Select:"),
                 choices=list("Introduction" = "About",
                              "Sample Size Planning" = "Sample",
                              "Interim Analysis" = "Interim"
                 ), selected="About")
    ),
    conditionalPanel(
      condition = ("input.Navigation != 'About'"),
        wellPanel(
          p(strong("Change/Set Parameters for Prior Probability:")),
          sliderInput("kk", p("Number of Strata"),
                      min=1, max=12, value = 3, step = 1),
          br(),
          sliderInput("lambda", p("Probability that activity is correlated among strata (lambda)"),
                      min=0, max=1, value = 0.5, step = 0.01),    
          br(),
          sliderInput("gamma", p("Per-stratum probability of drug activity (gamma)"),
                      min=0, max=1, value = 0.4, step = 0.01),
          br(),
          sliderInput("phi", p("Threshold probability for activity: p_hi"),
                      min=0, max=1, value = 0.25, step = 0.01),
          br(),
          uiOutput("plo_slider")
        )
      )
   ),
  
  mainPanel(
    conditionalPanel(
      condition = "input.Navigation == 'About'",
      tabsetPanel(
        tabPanel(h4("Introduction"),
                 span(h5("Basket Trials"),style="color:darkred"),
                 HTML('<P><em>Basket trials</em> are early phase II cancer clinical trials of one drug 
                      in a population of patients with histologic types and/or genomic 
                      variants thought to make them responsive to the drug. Basket trials are 
                      <em>discovery trials</em> rather than <em> hypothesis testing trials</em> where the 
                      objective is to discover rather then
                      test hypothesis on which histological type or genomic alterations sensitize the tumor to 
                      the drug. Each histological or genomic variant is referred to as a <em> strata</em>.
                      This web-based application implements a new method based on Bayesian principles for 
                      the design, monitoring and analysis of Basket trials.</p>'),
                  HTML('<p>This application currently supports <b> Sample size planning</b> and 
                    <b> Interim analysis</b> for Basket trials. The user specifies the number of strata and the 
                    prior probability that the drug is active 
                    in any particular stratum. Also, since, in general, there would be uncertainty on whether 
                    these strata are completely correlated 
                    or independent with regard to the distribution of drug activity, the user also
                    specifies a prior probability corresponding to correlation of drug activity among strata.
                    The more general situation where there are asymmetries among the strata with regard to the 
                    prior probability of drug activity will be implemented in the next version of this application</P>'),
                    
                 span(h5("Sample size planning for Basket trials"),style="color:darkred"),
                 HTML('<p> Sample size planning for Basket trials is guided through simulations. Patient accrual is  
                  simulated assuming user specified prevalence in strata (default is equal prevalence in all strata).
                  It is further assumed that interim 
                  analysis is conducted after accrual of a user specified number of evaluable subjects. 
                  At interim analysis, 
                  the posterior probability of activity for each stratum is computed using the pre-specified 
                  prior probability parameters. Accrual to a stratum is terminated if the posterior probability 
                  of activity becomes > threshold or < 1-threshold.This threshold can be set by the user at a 
                  high value, for example, the default is 0.8. </p>'),
                 HTML('<p> Simulations are done for each possible number of active strata (no strata active to all
                  strata active). The final output from simulations are the expected values of the sample size, 
                  true positve rate, false negative rate, false positive rate and true negative rate. 
                  Here, the false positive rate is defined as the expected number of false positive strata divided by 
                  the expected number of true negative strata, where expected values are taken with regard to 
                  sampling the data and with regard to the prior on the number of positive strata. The other values are
                  analogously defined.</p>'), 
                 
                span(h5("Interim analysis for Basket trials"),style="color:darkred"),
                HTML('<P>Interim analysis for Basket trials amounts to just computing the posterior probability
                of activity for each stratum at the interim analysis time point. This computation is based on the
                number of evaluable subjects and the number of responses for each stratum at the interim 
                analysis time.</p>'),
                                
                span(h5('Basic instructions for using this application'),style="color:darkred"), 
                HTML('<ol><li>This application works best with latest versions of Firefox or Chrome, 
                it has not been tested with other browsers.</li>
                <li>For sample size planning, select the <b>Sample Size Planning</b> option on the  
                panel on the left. Enter the parameters for the calculation of the prior probabilities 
                as well as the simulation related parameters. Click on the <b>Click to start simulations </b>
                button. Depending on the parameters set, the simulations may take a few minutes to complete. 
                The simulation time increases greatly with the increase in the number of strata. </li>
               <li>For interim analysis, select the <b>Interim Analysis</b> option on the  
                panel on the left. Enter the parameters for the calculation of the prior probabilities 
                as well as the number of subjects evaluable, and number of responses at the interim analysis 
                time. The posterior probabilites for each strata are computed from these data and are tabulated.</li>
               </ol>'),
                span(h5('About this application'),style="color:darkred"),
                HTML('<P>This web application is developed using the <A href="http://www.rstudio.com/shiny" target="_blank">Shiny package</A> 
                from <A href="http://www.rstudio.com/shiny" target="_blank">RStudio</A>.</p>'),
                HTML('<p>For questions and comments on this application, please contact 
                <a href="mailto:rsimon@mail.nih.gov">Dr. Richard Simon</a> 
                at the <a href="http://brb.nci.nih.gov" target="_blank">Biometric Research Program</a></p>'),
              value=1),id="panel1")
      ),
    
    conditionalPanel(
      condition = "input.Navigation == 'Interim'",
      tabsetPanel(
        tabPanel(span(h4("Interim Analysis"),style="color:darkred"),
          wellPanel(  
          uiOutput("interim"),
          #tags$th(type="text/css", "#interim th {display: table-header-group}"),
          helpText(span(h6("Select and overwrite to change the values in the cells"),style="color:red"))
          ),
          wellPanel(span(h5("Interim Analysis Results: Posterior Probability of Activity"), style="color:darkred"),
                   tableOutput("displayPostprob"),
                   htmlOutput("dtime")
          )
        )
      )),
    
    conditionalPanel(
      condition = "input.Navigation == 'Sample'",
        tabsetPanel(
          tabPanel(span(h4("Sample Size Planning"),style="color:darkred"),
          wellPanel(
            h5("Change/Set parameters for simulations for sample size planning:"),
            numericInput("maxs", p("Maximum total sample size:"), value=20),
            br(),
            numericInput("nblock", p("Interim analysis after how many evaulable subjects?"), value=5),
            br(),
            tags$div(title="This is 'delta'. Accural stops if posterior probability 
                            in a stratum becomes more than delta or < 1-delta",
              sliderInput("delta", p("Threshold posterior probability in a stratum to stop accrual"),
                    min=0.5, max=1, value = 0.8, step = 0.01)
            ),
            br(),
            numericInput("nrep", p("Number of simulations"), value=1000),
            checkboxInput("equal", "Assume equal prevalance among strata?", value=TRUE),
            conditionalPanel(
              condition = "!input.equal",
              uiOutput("prevalence"),
              helpText(span(h6("Select and overwrite to change values in the cells. Entered values to be > 0 and < 1 and sum to 1."),style="color:red"))
             ),
            actionButton(inputId="sstart", label = strong("Click to start simulations"))
          ),
        wellPanel(
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           id='progressIndicator',
                           img(src = "Process1.gif"),
                           span(h6("Working...Please Wait..."),style="color:red")
                           ),
          htmlOutput("simResults")           ### for simulation results output display 
         ),
   #tags$style(type="text/css", "h6 { color: red; }"),
   tags$head(tags$style(type="text/css",
                         ' #progressIndicator {',
                         ' position: fixed; top: 200px; right: 25px; width: 150px; height: 50px;',
                         ' padding: 10px; border: 2px solid #CCC; border-radius: 8px;',
                         '}'))
      )
    )
  )
  ) ### of mainpanel
  ))