library(shiny)
library(shinycssloaders)
library(rhandsontable)
library(dplyr)
library(Rcpp)
library(DT)
library(clinfun)
ui <- fluidPage(
  navbarPage(
   # "",
    
    tags$head(
      tags$style(
        HTML('
      .title {
        text-align: center;
        font-size: 24px;
        font-weight: bold;
      }
      .author {
        text-align: center;
        font-size: 14px; /* Changed font size to 14px */
        font-style: italic;
      }
    ')
      )
    ),
    div(class = "title", "An improved randomized two-stage Bayesian pick-the-winner
design for phase II clinical trials"),
    div(class = "author", HTML("Authors: Wanni Lei<sup>1</sup> | Maosen Peng<sup>1,2</sup> | Xi Kathy Zhou<sup>1</sup>")),
    div(class = "author", "1. Department of Population Health Sciences, Weill Cornell Medicine, New York, USA"),
    div(class = "author", "2. Department of Biostatistics and Data Science, University of Texas School of Public Health, Houston, TX, USA"),
    br(),
    tabPanel("Page 1",
             
             #titlePanel("Optimal and Minimax design"),
             sidebarPanel(
               textInput('Pa_0', 'Pa0: response rate of Arm A under H0', 0.1),
               textInput('Pb_0', 'Pb0: response rate of Arm B under H0', 0.1),
               textInput('Pa_1', 'Pa1: response rate of Arm A under H1', 0.1),
               textInput('Pb_1', 'Pb1: response rate of Arm B under H1', 0.4),
               textInput('alpha1', 'Type-I error', 0.1),
               textInput('beta1', 'Type-II error', 0.2),
               textInput("delta1","Delta: Arm B will be the winner if Pr(B >A) > δ", 0.8),
               #textInput("nmax","nmax: Maximum sample size for a single arm", 100),
               actionButton("runButton1", "Run")
             ),
             mainPanel(
               verbatimTextOutput("Parameters1"),
               verbatimTextOutput("Hypothesis1"),
               verbatimTextOutput("explain1"),
               tabsetPanel(
                 tabPanel("Optimal and Minimax design",withSpinner(DT::DTOutput('results1'))),
                 tabPanel("Reference", 
                          tags$a(href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5687612/", 
                                 "1. Chen, D. T., Huang, P. Y., Lin, H. Y., Chiappori, A. A., Gabrilovich, D. I., Haura,
E. B., ... & Gray, J. E. (2017).
A Bayesian pick-the-winner design in a randomized phase II clinical trial. Oncotarget, 8(51), 88376."),
                          br(),
                          tags$a(href = "https://www.sciencedirect.com/science/article/abs/pii/0197245689900159", 
                                 "2. Simon, R. (1989). Optimal two-stage designs for phase II clinical trials. Controlled clinical trials, 10(1), 1-10.")
                 )
               )
             )
    ),  
    tabPanel("Page 2",
             fluidPage(
               titlePanel("Trial Setting"),
               #verbatimTextOutput("Instruction"),
               # Rows 1-4 with 2 columns each
               sidebarPanel(width = 10,
              verbatimTextOutput("Instruction"),
                # Adding an editable table 
                rHandsontableOutput("tableInput"),
                 actionButton("runButton2", "Run")),
               mainPanel(
                 titlePanel("Operating Characteristics"),
                 
                 tabsetPanel(
                   tabPanel(
                            "Overall", DT::DTOutput('results2'),
                            verbatimTextOutput("Instruction2")),
                   tabPanel("Each arm under H0",
                           DT::DTOutput('results31'),
                            DT::DTOutput('results32'),
                            DT::DTOutput('results33'),
                            DT::DTOutput('results34'),
                            DT::DTOutput('results35')),
                   
                   tabPanel("Each arm under H1",
                            DT::DTOutput('results41'),
                            DT::DTOutput('results42'),
                            DT::DTOutput('results43'),
                            DT::DTOutput('results44'),
                            DT::DTOutput('results45'))
                 )
               )
               
               
               
             )
    ), ## tabPanel
    tabPanel("Page 3",
             fluidPage(
               titlePanel('Decision rule'),
               sidebarPanel(width = 8,
               selectInput("design_select", label = "Choose a design:",
                           choices = c("Optimal", "Minimax")),
               verbatimTextOutput("Instruction_page3"),
               actionButton("runButton3", "Run")),
               mainPanel(
                 DT::DTOutput("results_dr")
                 )  # mainPanel
      ) #fluidPage
    ) #tabPanel
  ) ## navbarPage
) ## fluidPage


