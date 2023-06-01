library(shiny)
library(shinydisconnect)
library(shinyjs)
library(tidyverse)
library(ckanr)
library(jsonlite)
library(TissueEnrich)
library(flextable)

# Define the list of genes and annotations
genes <- readRDS("../data/genes.rds")
gene_list <- genes$symbol
gene_context <- genes$context


annotations <- c("tissue_expn",
                 "cell_type_expn",
                 "phewas",
                 "focused_phewas",
                 "maf",
                 "mendelian_phenotypes",
                 'mouse_phenotypes(enexpr(gene))',
                 "MF",
                 "eQTL",
                 "pathways",
                 "tool_compounds")

annotation_labels <- c("(Gene) Tissue-specific expression (RNA and protein)",
                      "(Gene) Cell-type-specific expression for a tissue of interest",
                      "(Variant) PheWAS readouts, general",
                      "(Variant) PheWAS readouts, tissue-focused",
                      "(Variant) MAF in various populations",
                      "(Gene and Variant) Mendelian phenotypes",
                      "(Gene) Mouse phenotypes",
                      "(Gene) Molecular function of gene",
                      "(Variant) eQTL effect on RNA levels in specific tissues",
                      "(Gene) Pathway",
                      "(Gene) Tool compounds")
availability <- rep(FALSE, length(annotation_labels))
availability[c(1,2,4,5,7,8,9)] <- TRUE

annotation_modules <- tibble(label = annotation_labels,
                             module = annotations,
                             module_index = seq_along(annotations),
                             availability = availability)

annotation_list <- annotation_modules %>%
  dplyr::filter(availability == TRUE) %>%
  dplyr::select(label, module_index) %>%
  tibble::deframe()



css <-
  "
a#maze-shiny-apps {
display: none;
}
body {
background: #fff;
}
.container {
margin: 0;
padding: 15px;
}
.green {
background: .green;
}
#gene-container .selectize-input {
  font-size: 20px;
line-height: 20px;
}

#gene-container .selectize-dropdown {
font-size: 16px;
line-height: 22px;
}

#gene-container .selectize-dropdown-content {
max-height: 225px;
padding: 0;
}

#gene-container .selectize-dropdown-content .option {
border-bottom: 1px dotted #ccc;
}

#annotation-container .selectize-input {
  font-size: 20px;
line-height: 20px;
}

#annotation-container .selectize-dropdown {
font-size: 16px;
line-height: 22px;
}

#annotation-container .selectize-dropdown-content {
max-height: 225px;
padding: 0;
}

#annotation-container .selectize-dropdown-content .option {
border-bottom: 1px dotted #ccc;
}


#submitSelection {
font-weight: bold;
font-size: 18px;
padding: 5px 20px;
}

#helpText {
font-size: 18px;
}

#btn {
font-size: 20px;
}
"


ui <- fluidPage(
  shinydisconnect::disconnectMessage2(),
  shinyjs::useShinyjs(),
  tags$style(css),
  titlePanel("Select Inputs"),
  fluidRow(
    column(
      width = 6,
      div(id = "annotation-container",
          selectInput(inputId = "context", label = "Gene:", choices = c("all", gene_context),
                      selected = 1, multiple = FALSE, selectize = TRUE, width = "100%")
      )
    )
  ),
  fluidRow(
    column(
      width = 3,
      div(id = "gene-container",
        selectInput(inputId = "gene", label = "Gene:", choices = gene_list,
                    selected = 1, multiple = FALSE, selectize = TRUE, width = "100%")
        )
    ),
    column(
      width = 9,
      div(id = "annotation-container",
          selectInput(inputId = "annotation", label = "Annotation:", choices = annotation_list,
                      selected = 1, multiple = FALSE, selectize = TRUE, width = "100%")
      )
    )
  ),
  actionButton("submitSelection", "Deliver", class = "btn-success"),
  shiny::hr(),

  htmlOutput("selection"),

)

# Define the server logic
server <- function(input, output, session) {
  observeEvent(input$context, {

    # Perform filtering based on the selected category
    if(input$context == "all")
      selected_context <- genes %>% pull(context)
    else
      selected_context <- c(input$context)

    filtered_genes <- genes %>%
      dplyr::filter(context %in% selected_context) %>%
      dplyr::select(symbol, gene_id) %>%
      tibble::deframe()

    # Update the choices in the item selectInput
    updateSelectInput(session, "gene", choices = filtered_genes)
  })


    # Define the reactive expression to generate output based on the selected gene and annotation
  selected_output <- reactive({

    gene <- input$gene
    annotation <- input$annotation

    # Example output: displaying the selected gene and annotation
    output_text <- paste("Selected Gene:", gene, "<br>")
    output_text <- paste(output_text, "Selected Annotation:", annotation, "<br>")

    if(annotation == 7) {
      isolate(
        eval(parse(text = annotation_modules$module[as.numeric(input$annotation)]))
      )
    } else {
      isolate(
        eval(parse(text = "ls()"))
      )
      output_text
    }
  })

  # Render the output as a table or a figure, based on the selected gene and annotation
  output$selection <- renderUI({
    selected_output()
  })
}

# Run the application
shinyApp(ui, server)
