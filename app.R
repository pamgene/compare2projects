library(shiny)
library(shinyjs)

ui <- fluidPage(
    useShinyjs(),
    fluidPage(
        sidebarLayout(
            sidebarPanel(
                h1("Parameters"),
                dateInput("date", "Report Date"),
                textInput("data1_name", "Data 1 name", value = 'data1'),
                textInput("data2_name", "Data 2 name", value = 'data2'),
                numericInput('pcutoff', "P value cutoff", value = 0.05, min = 0.0001, max = 1),
                numericInput('fscorecutoff', "Final Score cutoff", value = 1.3, min = 0, max = 12),
                fluidRow(
                  disabled(actionButton("go", "GO", class = "btn-lg btn-success")),
                  br(),
                  uiOutput("download")
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
)

server <- function(input, output, session) {
  make_datafolders <- function(folders) {
    for (folder in folders){
      if (!dir.exists(folder)) {
        dir.create(folder)
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
      make_datafolders(c('data1', 'data2', 'results'))
      move_datafiles(input$data1_files, "data1")
      move_datafiles(input$data2_files, "data2")
      
      enable("go")
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
      source("main.R")
    })
    showNotification("Results complete!", type = "message")
    make_result_zip()

    output$download <- renderUI({
      downloadButton("download_results", "Download results")
    })
    
    outputOptions(output, "download_results", suspendWhenHidden = FALSE)
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