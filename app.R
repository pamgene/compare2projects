library(shiny)
library(shinyjs)

ui <- fluidPage(
    useShinyjs(),
    tags$div(class = "jumbotron text-center", style = "margin-bottom:0px;margin-top:0px",
           tags$h1(class = 'jumbotron-heading', style = 'margin-bottom:0px;margin-top:0px', 'Compare two projects'),
           p('Compare peptide and kinase results of repeated projects.')
    ),
    sidebarLayout(
        sidebarPanel(
            h3("Usage"),
            p("Save phosphosite and kinase analysis result files (BN or Tercen) of the two projects in separate folders.
              In each folder, the files should have the exact same naming.
              In each file, the namespaces should be exactly the same.
              Name files the same as for automated reporting."),
            h2("Parameters"),
            p("Choose the output format for your report below. 'Summary stats' provides an overview, while 'Per comparison' gives detailed results for each comparison."),
            radioButtons(
              "output_type", "Output_type", choices = c("Summary stats" = "summary", "Per comparison" = "per_comp")),
            dateInput("date", "Report Date"),
            textInput("data1_name", "Data 1 name", value = 'data1'),
            textInput("data2_name", "Data 2 name", value = 'data2'),
            radioButtons(
              "uka_version", "UKA version", choices = c("All vs all" = "all_vs_all", "Single comparison" = "single_comp")),
            numericInput('pcutoff', "P value cutoff", value = 0.05, min = 0.0001, max = 1),
            numericInput('fscorecutoff', "Final Score cutoff", value = 1.3, min = 0, max = 12),
            fluidRow(
              disabled(actionButton("go", "GO", class = "btn-lg btn-success")),
              br(),
              disabled(downloadButton("download_results", "Download results"))
            )
        ),
        mainPanel(
            br(),
            h1("Files"),
            fileInput("data1_files", "Upload .csv or .txt files from data1.", multiple = TRUE, accept = c(".csv", ".txt")),
            fileInput("data2_files", "Upload .csv or .txt files from data2.", multiple = TRUE, accept = c(".csv", ".txt")),
            h2("Data 1 files"),
            tableOutput('files1'),
            h2("Data 2 files"),
            tableOutput('files2')
        )
    )
)

server <- function(input, output, session) {
  make_datafolders <- function(folders) {
    for (folder in folders){
      if (!dir.exists(folder)) {
        dir.create(folder, recursive = T)
      }
    }
  }
  
  remove_datafolders <- function(folders) {
    for (folder in folders){
      if (dir.exists(folder)) {
        unlink(folder, recursive = T)
      }
    }
  }
  

  make_result_zip <- function() {
    file_name <- paste0(input$data2_name, " vs ", input$data1_name, "_results_", format(input$date, "%y%m%d"), ".zip")
    zip(file_name, c("results/"))
    
  }
  
  move_datafiles <- function(datafiles, folder){
    for (i in 1:nrow(datafiles)) {
      aFile <- datafiles[i,]
      file.copy(from = aFile$datapath, to = file.path(folder, aFile$name))
    }
  }

  observe({
    if (is.null(input$data1_files) & is.null(input$data2_files)) {
      return()
    } else if (!is.null(input$data1_files) & !is.null(input$data2_files)) {
      remove_datafolders(c('results', 'data1', 'data2'))
      make_datafolders(c('results', 'data1', 'data2'))
      move_datafiles(input$data1_files, "data1")
      move_datafiles(input$data2_files, "data2")
      enable("go")
      disable("download_results")
    }
  })
  
  output$files1 <- renderTable(input$data1_files)
  output$files2 <- renderTable(input$data2_files)
  
  observeEvent(input$go, {

    withProgress(message = "Working on results...", {
      data1 <<- 'data1'
      data2 <<- 'data2'
      data1_name <<- input$data1_name
      data2_name <<- input$data2_name
      fscorecutoff <<- input$fscorecutoff
      pcutoff <<- input$pcutoff
      uka_version <<- input$uka_version 
      output_type <<- input$output_type
      source("main.R")
    })
    showNotification("Results complete!", type = "message")
    make_result_zip()
    enable("download_results")
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste0(input$data2_name, " vs ", input$data1_name, "_results_", format(input$date, "%y%m%d"), ".zip")
    },
    content = function(file) {
      file.copy(from = paste0(input$data2_name, " vs ", input$data1_name, "_results_", format(input$date, "%y%m%d"), ".zip"),
                to = file)
    },
    contentType = "application/zip"
  )
}

# Run the application
shinyApp(ui = ui, server = server)
