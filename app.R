suppressMessages({
  library(shiny)
  library(DT)
  library(factoextra)
  library(tidyverse)
  library(patchwork)
  library(DESeq2)
  library(purrr)
  library(shinycssloaders)
  library(shinythemes)
  library(ggrepel)
  library(gridExtra)
})

# Define UI for application that draws a histogram
ui <- navbarPage("MolEvol – DGE", footer = fluidRow(column(2, tags$b("OVERALL PROGRESS:")),
                                                    column(1, uiOutput("footer1")),
                                                    column(1, uiOutput("footer2")),
                                                    column(1, uiOutput("footer3")),
                                                    column(1, uiOutput("footer4")),
                                                    column(1, uiOutput("footer5"))), 
  theme = shinytheme("yeti"),
  tabPanel(
    "Data upload",
    fluidRow(
      column(
        3,
        h3("Upload a file with counts and metadata to start"),
        wellPanel(
          fileInput("countfile", "Select counts file"),
          fileInput("metafile", "Select metadata file")
        )
      ),
      column(
        9,
        conditionalPanel(
          condition = "output.counts",
          h3("Counts table preview")
        ),
        shinycssloaders::withSpinner(DT::dataTableOutput("counts")),
        conditionalPanel(
          condition = "output.meta",
          h3("Metadata preview")
        ),
        shinycssloaders::withSpinner(DT::dataTableOutput("meta")),
      ),
    )
  ),
  tabPanel(
    "Data filtering",
    fluidRow(
      column(
        3,
        h3("Filter genes"),
        h4("Remove loci for which:"),
        wellPanel(
          numericInput("filtnum", value = NA, label = "This many samples:"),
          numericInput("filtreads", value = NA, label = "have this many or fewer mapped reads:")
        ),
        h3("Apply filter"),
        h4("Click to continue with filtered dataset"),
        actionButton("filter", "Filter now"),
        conditionalPanel(
          condition = "input.filter", h3("Download"),
          downloadButton("filtdown", "Download filtered counts table")
        )
      ),
      column(
        9, 
        shinycssloaders::withSpinner(plotOutput("filtplot", width = "1000px", height = "600px")),
        conditionalPanel(condition = "output.counts", h3("Preview of filtered dataset")),
        DT::dataTableOutput("filtlook")
      )
    )
  ),
  tabPanel("PCA", shinycssloaders::withSpinner(
    plotOutput("PCA", width = "1200px", height = "600px")
  )),
  tabPanel(
    "DGE",
    fluidRow(
      column(
        3,
        h3("Determine differentially expressed genes"),
        wellPanel(
          uiOutput("group1"), uiOutput("group2"),
          actionButton("calc.dge", "Start DESeq2 calculation")
        ),
        conditionalPanel(
          condition = "output.samples",
          h3("Download results"),
          downloadButton("dgedown", "Download all DESeq2 results"),
          hr(),
          sliderInput("pvalue", "Set pvalue cutoff", value = 0.05, min = 0, max = 0.1, step = 0.01),
          downloadButton("siggenes", "Download DGEs")
          
        )
      ),
      column(
        9,
        shinycssloaders::withSpinner(DT::dataTableOutput("samples")),
        conditionalPanel(
          condition = "output.samples",
          shinycssloaders::withSpinner(plotOutput("MA", width = "50%"))
        )
      )
    )
  ),
  tabPanel(
    "Count plots",
    fluidRow(
      column(
        3,
        h3("Plot expression levels for genes of interest"),
        wellPanel(
          textAreaInput("genelist", "Enter one gene ID per line", height = "400px"),
          actionButton("countplot", "Plot expression levels")
        ),
        conditionalPanel(
          condition = "output.count", sliderInput("plotsize", "Change plot size",
            min = 1, max = 50, value = 15, step = 1, ticks = F
          ),
          downloadButton("countdown", "Download plot as pdf")
        )
      ),
      column(
        9,
        shinycssloaders::withSpinner(plotOutput("count", width = "800px", height = "800px"))
      )
    )
  ),
  hr(),
  
)

server <- function(input, output, server) {
  options(shiny.maxRequestSize = 500 * 1024^2)

  countstab <- reactive({
    req(input$countfile$datapath)
    read.table(input$countfile$datapath, header = TRUE, sep = "\t")
  })

  output$counts <- DT::renderDataTable(
    {
      countstab()
    },
    selection = "none",
    options = list(
      searching = FALSE
    )
  )

  metatab <- reactive({
    req(input$metafile$datapath)
    read.table(input$metafile$datapath, header = TRUE, sep = "\t")
  })

  output$meta <- DT::renderDataTable(
    {
      metatab()
    },
    selection = "none",
    options = list(
      searching = FALSE
    )
  )

  filttab <- reactiveValues()

  observeEvent(countstab(), {
    filttab$filt <- countstab()
  })

  observeEvent(input$filter, {
    filttab$filt <- countstab()[rowSums(countstab() <= input$filtreads) <= input$filtnum, ]
  })

  observeEvent(input$filter, {
    showModal(modalDialog(
      title = "Filer applied",
      paste(nrow(filttab$filt), "genes remaining."),
      easyClose = TRUE,
      footer = modalButton("Dismiss")
    ))
  })

  output$filtlook <- DT::renderDataTable(
    {
      filttab$filt
    },
    caption = "Preview of filtered dataset",
    selection = "none"
  )

  output$filtdown <- downloadHandler(
    filename = "filteredcounts.csv",
    content = function(file) {
      write.table(filttab$filt,
        quote = FALSE, row.names = FALSE, sep = ",",
        file = file, append = FALSE
      )
    }
  )

  output$filtplot <- renderPlot({
    req(countstab())
    req(input$filtreads)
    req(input$filtnum)

    filttab <- countstab()[rowSums(countstab() <= input$filtreads) <= input$filtnum, ]

    pA <- countstab() %>%
      mutate(total = rowSums(across(where(is.numeric)))) %>%
      ggplot(aes(x = total)) +
      geom_density() +
      xlim(c(0, 500)) +
      ggtitle(paste("Before filtering:", nrow(countstab()), "genes")) +
      xlab("Observed counts") +
      theme_light(base_size = 20)

    pB <- filttab %>%
      mutate(total = rowSums(across(where(is.numeric)))) %>%
      ggplot(aes(x = total)) +
      geom_density() +
      xlim(c(0, 500)) +
      xlab("Observed counts") +
      ggtitle(paste("After filtering:", nrow(filttab), "genes")) +
      theme_light(base_size = 20)

    pA / pB
  })

  pcaplot <- reactive({
    req(metatab())
    req(filttab$filt)
    rownames(filttab$filt) <- c()
    tab1 <- filttab$filt %>%
      column_to_rownames(var = colnames(filttab$filt)[1])
    tab2 <- metatab() %>%
      column_to_rownames(var = colnames(metatab())[1])
    colnames(tab2) <- "condition"
    pca.res <- prcomp(tab1, scale. = TRUE)
    fviz_eig(pca.res)
    fviz_contrib(pca.res, "var", axes = c(1, 2))
    p1 <- fviz_pca_var(pca.res,
      col.var = tab2$condition,
      geom = c("text", "point"),
      addEllipses = FALSE,
      labelsize = 4,
      pointsize = 6,
      pointshape = 19,
      repel = TRUE,
      title = "PCA"
    ) +
      theme_light(base_size = 20) +
      theme(legend.title = element_blank())

    tab1.1 <- as.data.frame(cmdscale(dist(scale(t(tab1))))) %>%
      rownames_to_column(var = "sample")
    tab2.1 <- tab2 %>%
      rownames_to_column(var = "sample")
    p2 <- full_join(tab2.1, tab1.1) %>%
      ggplot(aes(x = V1, y = V2, col = condition, label = sample)) +
      geom_point(size = 6) +
      geom_text_repel(size = 4, force = 50) +
      theme_light(base_size = 20) +
      theme(
        axis.title = element_blank(),
        legend.title = element_blank()
      ) +
      ggtitle("MDS")

    p1 | p2
  })

  output$PCA <- renderPlot({
    pcaplot()
  })
  
  output$group1 <- renderUI({
    req(metatab())
    req(countstab())
    selectInput(
      inputId = "group1",
      label = "Select group to compare to",
      choices = unique(metatab()[, 2])
    )
  })

  output$group2 <- renderUI({
    req(metatab())
    req(countstab())
    selectInput(
      inputId = "group2",
      label = "Select group to compare with",
      choices = c("", unique(metatab()[, 2]))
    )
  })

  compare <- reactiveValues()
  observeEvent(input$calc.dge, {
    compare$group1 <- input$group1
    compare$group2 <- input$group2
  })

  deg.table <- reactive({
    req(metatab())
    req(countstab())
    req(compare$group1)
    req(compare$group2)
    group1 <- compare$group1
    group2 <- compare$group2
    validate(need(group1 != group2, "Please select two different groups for comparison."))
    validate(need(group1, "Please select two different groups for comparison."))
    validate(need(group2, "Please select two different groups for comparison."))
    rownames(filttab$filt) <- c()
    tab1 <- filttab$filt %>%
      column_to_rownames(var = colnames(filttab$filt)[1])
    tab2 <- metatab() %>%
      column_to_rownames(var = colnames(metatab())[1])
    colnames(tab2) <- "condition"
    tab2$condition <- factor(tab2$condition)
    tab2$condition <- relevel(tab2$condition, ref = group1)
    samples_filt <- rownames(filter(tab2, condition %in% c(group1, group2)))
    dds <- DESeqDataSetFromMatrix(
      countData = round(select(tab1, samples_filt)),
      colData = filter(tab2, condition %in% c(group1, group2)),
      design = ~condition
    )
    DESeq(dds)
  })

  restab <- reactive({
    req(deg.table())
    results(deg.table(), format = "DataFrame", tidy = TRUE) %>%
      select(row, log2FoldChange, pvalue, padj)
  })

  output$samples <- DT::renderDataTable(
    {
      DT::datatable(restab(),
                    caption = htmltools::tags$caption(
                      style = 'font-size: 24px; font-weight: bold',
                      paste(compare$group1, "vs.", compare$group2)),
                    selection = "none") %>% 
        DT::formatStyle("padj", target = "row", 
                        backgroundColor = styleInterval(input$pvalue, c("#b9db92", "white")))
    }
  )

  output$dgedown <- downloadHandler(
    filename = paste0(compare$group1, ".vs.", compare$group2, ".dge.csv"),
    content = function(file2) {
      write.table(restab(),
        quote = FALSE, row.names = FALSE, sep = ",",
        file = file2, append = FALSE
      )
    }
  )

  output$siggenes <- downloadHandler(
    filename = paste0(compare$group1, ".vs.", compare$group2, ".", input$pvalue, ".SIGN.csv"),
    content = function(file2) {
      write.table(filter(restab(), padj <= input$pvalue),
                  quote = FALSE, row.names = FALSE, sep = ",",
                  file = file2, append = FALSE
      )
    }
  )
  
  output$MA <- renderPlot({
    req(deg.table())
    plotMA(deg.table(),
      main = paste(compare$group1, "vs.", compare$group2),
      colSig = "red"
    )
  })

  genelist <- eventReactive(input$countplot, {
    genelist <- as.character(stringr::str_split(input$genelist, "\n", simplify = TRUE))
    genelist <- unique(genelist[genelist != ""])
    genelist
  })

  countP <- reactive({
    req(genelist())
    genes <- genelist()
    genelist <- list()
    for (i in 1:length(genes)) {
      df1 <- plotCounts(deg.table(),
        gene = genes[i], intgroup = "condition",
        returnData = TRUE
      )
      df1$geneID <- genes[i]
      genelist[[i]] <- df1
    }
    genelist %>%
      purrr::reduce(full_join) %>%
      ggplot(aes(x = condition, y = count, colour = condition)) +
      geom_point(size = input$plotsize / 3) +
      facet_wrap(~geneID, scales = "free_y") +
      theme_light(base_size = input$plotsize)
  })

  output$count <- renderPlot({
    countP()
  })

  output$countdown <- downloadHandler(
    filename = "counts.pdf",
    content = function(file3) {
      pdf(file3, width = 800 / 72, height = 800 / 72)
      grid.arrange(countP(), ncol = 1)
      dev.off()
    }
  )
  
  textcolors <- reactiveValues(data = "grey", 
                               filter = "gray",
                               pca = "gray",
                               dge = "gray",
                               plots = "gray")
  
  observeEvent({
    countstab() 
    metatab()},{
    textcolors$data <- "green"
  })
  
  observeEvent(input$filter,{
    textcolors$filter <- "green"
  })
  
  observeEvent(pcaplot(),{
    textcolors$pca <- "green"
  })
  
  observeEvent(restab(),{
    textcolors$dge <- "green"
  })
  
  observeEvent(input$countplot,{
    textcolors$plots <- "green"
  })
  
  output$message1 <- renderText({"| DATA > "})
  output$message2 <- renderText({"FILTERING > "})
  output$message3 <- renderText({"PCA > "})
  output$message4 <- renderText({"DGE > "})
  output$message5 <- renderText({"PLOTS |"})
  
  output$footer1 <- renderUI({
    span(textOutput("message1"), style = paste("color:", textcolors$data))
  })
  output$footer2 <- renderUI({
    span(textOutput("message2"), style = paste("color:", textcolors$filter))
  })
  output$footer3 <- renderUI({
    span(textOutput("message3"), style = paste("color:", textcolors$pca))
  })
  output$footer4 <- renderUI({
    span(textOutput("message4"), style = paste("color:", textcolors$dge))
  })
  output$footer5 <- renderUI({
    span(textOutput("message5"), style = paste("color:", textcolors$plots))
  })
}


# Run the application
shinyApp(ui = ui, server = server)
