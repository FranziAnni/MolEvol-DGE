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
  library(ggpubr)
  library(ggheatmap)
  library(ggplotify)
})

# Define UI for application that draws a histogram
ui <- navbarPage("MolEvol – DGE",
  footer = fluidRow(
    hr(),
    column(2, HTML("<center><b>OVERALL<br>PROGRESS</b></center>")),
    column(1, uiOutput("footer1")),
    column(1, uiOutput("footer2")),
    column(1, uiOutput("footer3")),
    column(1, uiOutput("footer4")),
    column(1, uiOutput("footer5"))
  ),
  theme = shinytheme("yeti"),
  tabPanel(
    "Data upload",
    fluidRow(
      column(
        3,
        h3("Upload a file with counts and metadata to start"),
        wellPanel(
          div(style = "float:right", actionLink("filebutton", "", icon = icon("info-circle", lib = "font-awesome"))),
          fileInput("countfile", "Select counts file"),
          radioButtons("countsep",
            label = "Column separator",
            inline = TRUE,
            choices = c(
              "<TAB>" = "\t",
              "comma" = ",",
              "semi-colon" = ";"
            ),
            selected = "\t"
          ),
          hr(),
          fileInput("metafile", "Select metadata file"),
          radioButtons("metasep",
            label = "Column separator",
            inline = TRUE,
            choices = c(
              "<TAB>" = "\t",
              "comma" = ",",
              "semi-colon" = ";"
            )
          )
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
        h4("Keep only loci for which:"),
        wellPanel(
          div(style = "float:right", actionLink("filterbutton", "", icon = icon("info-circle", lib = "font-awesome"))),
          numericInput("filtnum", value = NA, label = "Not more than this many samples:"),
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
  tabPanel(
    "EDA",
    fluidRow(
      column(
        3,
        h3("Exploratory data analysis"),
        h4("These plots will give you an idea of how your data is structured"),
        sliderInput("kmeans", "Number of gene clusters for heatmap", 50, 2000, 100, step = 50),
        actionButton("plotpca", "Plot now")
      ),
      column(
        9,
        shinycssloaders::withSpinner(
          plotOutput("EDA", width = "1000px", height = "900px")
        ),
        conditionalPanel(
          condition = "output.EDA", downloadButton("pcadown", "Download this plot")
        )
      )
    )
  ),
  tabPanel(
    "DGE",
    fluidRow(
      column(
        3,
        h3("Determine differentially expressed genes"),
        wellPanel(
          uiOutput("group1"), uiOutput("group2"),
          conditionalPanel(
            condition = "output.group2",
            actionButton("calc.dge", "Start DESeq2 calculation")
          )
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
        shinycssloaders::withSpinner(DT::dataTableOutput("samples"),
          image = sample(list.files("www"), 1),
          image.height = "200px"
        ),
        hr(),
        conditionalPanel(
          condition = "output.samples",
          shinycssloaders::withSpinner(plotOutput("MA", width = "100%")),
          downloadButton("madown", "Download this plot")
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
  tabPanel("About", fluidRow(column(5, includeMarkdown("about.md"))))
)

server <- function(input, output, server, session) {

  # increase maximum upload size
  options(shiny.maxRequestSize = 500 * 1024^2)

  # shut down R session when browser window is closed
  session$onSessionEnded(function() {
    stopApp()
  })

  observeEvent(input$filebutton, {
    showModal(modalDialog(
      title = "Format of input files",
      "The counts file should have 1 line per gene and 1 column per sample.
    The first column should contain the gene names. The metadata file
    should list the sample names in the first column, and the variable
    assignment in the second column. Both files should have column headers. 
    Change the column delimiter according to your data.  
    \n\nPlease note that the gene names in counts and metadata files must match!",
      size = "m",
      icon = icon("info-circle", lib = "font-awesome")
    ))
  })

  countstab <- reactive({
    req(input$countfile$datapath)
    read.table(input$countfile$datapath, header = TRUE, sep = input$countsep)
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
    read.table(input$metafile$datapath, header = TRUE, sep = input$metasep)
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

  observeEvent(metatab(), {
    filttab$meta <- metatab()
  })

  observeEvent(input$filterbutton, {
    showModal(modalDialog(
      title = "Filtering your data",
      HTML("It is often advisable to filter out genes with very low counts and/or
      those loci that only have mapped reads in a minority of your samples.
      Here you can specify a cutoff for filtering loci by the number of mapped
      reads per sample. The filtering will be stricter with lower sample
      numbers and higher read numbers.<br><br>For example, to keep only loci that
      have mapped reads in all of the samples, you would specify '0' and '0'
      for samples and reads, respectively. If you only wanted loci for which at
      least half of your 10 samples have mapped reads, you would specify '5' and '0'."),
      size = "m",
      icon = icon("info-circle", lib = "font-awesome")
    ))
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
      xlim(c(0, 200)) +
      ggtitle(paste("Before filtering:", nrow(countstab()), "genes")) +
      xlab("Observed counts") +
      theme_light(base_size = 18)

    pB <- filttab %>%
      mutate(total = rowSums(across(where(is.numeric)))) %>%
      ggplot(aes(x = total)) +
      geom_density() +
      xlim(c(0, 200)) +
      xlab("Observed counts") +
      ggtitle(paste("After filtering:", nrow(filttab), "genes")) +
      theme_light(base_size = 18)

    pA / pB
  })

  pcaplot <- eventReactive(input$plotpca, {
    req(filttab$meta)
    req(filttab$filt)
    rownames(filttab$filt) <- c()
    tab1 <- filttab$filt %>%
      column_to_rownames(var = colnames(filttab$filt)[1])
    tab2 <- filttab$meta %>%
      column_to_rownames(var = colnames(filttab$meta)[1])
    colnames(tab2) <- "sample.type"
    pca.res <- prcomp(tab1, scale = TRUE)

    p1 <- fviz_pca_biplot(pca.res,
      invisible = "ind",
      col.var = tab2$sample.type,
      geom.var = c("text", "point"),
      addEllipses = FALSE,
      labelsize = 4,
      pointsize = 6,
      pointshape = 19,
      repel = TRUE,
      title = "PCA",
      col.circle = NA
    ) +
      theme_light(base_size = 18) +
      theme(legend.title = element_blank())

    tab1.1 <- as.data.frame(cmdscale(dist(scale(t(tab1))))) %>%
      rownames_to_column(var = "sample")
    tab2.1 <- tab2 %>%
      rownames_to_column(var = "sample")

    p2 <- full_join(tab2.1, tab1.1) %>%
      ggplot(aes(x = V1, y = V2, col = sample.type, label = sample)) +
      geom_point(size = 6) +
      geom_text_repel(size = 4, force = 50) +
      theme_light(base_size = 18) +
      theme(
        axis.title = element_blank(),
        legend.title = element_blank()
      ) +
      ggtitle("MDS")

    hmcolors <- c(
      "#000004FF", "#1B0C42FF", "#4B0C6BFF", "#781C6DFF", "#A52C60FF",
      "#CF4446FF", "#ED6925FF", "#FB9A06FF", "#F7D03CFF", "#FCFFA4FF"
    )

    set.seed(3)
    kres <- kmeans(as.matrix(tab1), input$kmeans)
    gg_color_hue <- function(n) {
      hues <- seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    vec1 <- gg_color_hue(length(unique(tab2$sample.type)))
    names(vec1) <- sort(unique(tab2$sample.type))
    list1 <- list()
    list1$sample.type <- vec1
    plotlist <- ggheatmap(
      kres$centers,
      legendName = "",
      border = "white",
      scale = "row",
      color = hmcolors,
      cluster_rows = T,
      cluster_cols = T,
      annotation_cols = tab2,
      annotation_color = list1
    )

    plotlist$plotlist[[1]] <- plotlist$plotlist[[1]] +
      theme_light(base_size = 18) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()
      )

    plotlist$plotlist[[2]] <- plotlist$plotlist[[2]] +
      theme_void(base_size = 18) +
      theme(legend.title = element_blank())

    plotlist$plotlist[[4]] <- plotlist$plotlist[[4]] +
      theme_void(base_size = 18) + ggtitle("Heatmap\n")

    p3 <- ggplotify::as.ggplot(plotlist)

    layout <- "
      ##CCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      AACCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      BBCCC
      ##CCC
      "

    p1 + p2 + p3 + plot_layout(design = layout)
  })

  output$EDA <- renderPlot({
    pcaplot()
  })

  output$pcadown <- downloadHandler(
    filename = "pca.mds.pdf",
    content = function(file5) {
      pdf(file5, width = 1000 / 72, height = 800 / 72)
      print(pcaplot())
      dev.off()
    }
  )

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
    colnames(tab2) <- "sample.type"
    tab2$sample.type <- factor(tab2$sample.type)
    tab2$sample.type <- relevel(tab2$sample.type, ref = group1)
    samples_filt <- rownames(filter(tab2, sample.type %in% c(group1, group2)))
    dds <- DESeqDataSetFromMatrix(
      countData = round(select(tab1, samples_filt)),
      colData = filter(tab2, sample.type %in% c(group1, group2)),
      design = ~sample.type
    )
    DESeq(dds)
  })

  restab <- reactive({
    req(deg.table())
    results(deg.table(), format = "DataFrame", tidy = TRUE) %>%
      select(row, log2FoldChange, pvalue, padj)
  })

  output$samples <- DT::renderDataTable({
    DT::datatable(restab(),
      caption = htmltools::tags$caption(
        style = "font-size: 24px; font-weight: bold",
        paste(compare$group1, "vs.", compare$group2)
      ),
      selection = "none"
    ) %>%
      DT::formatStyle("padj",
        target = "row",
        backgroundColor = styleInterval(input$pvalue, c("#b9db92", "white"))
      )
  })

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

  maplot <- reactive({
    req(deg.table())
    p1 <- ggpubr::ggmaplot(results(deg.table()),
      size = 1,
      fdr = input$pvalue,
      top = 0
    ) +
      theme_light(base_size = 18) +
      ggtitle("MA plot",
        subtitle = paste(compare$group1, "vs.", compare$group2)
      )

    p2 <- restab() %>%
      ggplot(aes(x = padj)) +
      geom_histogram(bins = 100) +
      geom_vline(xintercept = input$pvalue, col = "red") +
      theme_light(base_size = 18) +
      ggtitle("Histogram of p-values",
        subtitle = paste(compare$group1, "vs.", compare$group2)
      )

    p1 | p2
  })

  output$MA <- renderPlot({
    maplot()
  })

  output$madown <- downloadHandler(
    filename = paste0(compare$group1, ".vs.", compare$group2, ".", input$pvalue, ".pdf"),
    content = function(file4) {
      pdf(file4, width = 1200 / 72, height = 600 / 72)
      print(maplot())
      dev.off()
    }
  )

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
        gene = genes[i], intgroup = "sample.type",
        returnData = TRUE
      )
      df1$geneID <- genes[i]
      genelist[[i]] <- df1
    }
    genelist %>%
      purrr::reduce(full_join) %>%
      ggplot(aes(x = sample.type, y = count, colour = sample.type)) +
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
      print(countP())
      dev.off()
    }
  )

  textcolors <- reactiveValues(
    data = "#D8D8D8",
    filter = "#D8D8D8",
    pca = "#D8D8D8",
    dge = "#D8D8D8",
    plots = "#D8D8D8"
  )

  observeEvent(
    {
      countstab()
      metatab()
    },
    {
      textcolors$data <- "#b9db92"
    }
  )

  observeEvent(input$filter, {
    textcolors$filter <- "#b9db92"
  })

  observeEvent(pcaplot(), {
    textcolors$pca <- "#b9db92"
  })

  observeEvent(restab(), {
    textcolors$dge <- "#b9db92"
  })

  observeEvent(input$countplot, {
    textcolors$plots <- "#b9db92"
  })

  output$message1 <- renderText({
    "DATA"
  })
  output$message2 <- renderText({
    "FILTER"
  })
  output$message3 <- renderText({
    "EDA "
  })
  output$message4 <- renderText({
    "DGE"
  })
  output$message5 <- renderText({
    "PLOTS"
  })

  output$footer1 <- renderUI({
    div(textOutput("message1"), style = paste("background-color:", textcolors$data, "; padding: 10px; border-radius: 20px ;text-align: center;"))
  })
  output$footer2 <- renderUI({
    div(textOutput("message2"), style = paste("background-color:", textcolors$filter, "; padding: 10px; border-radius: 20px ;text-align: center;"))
  })
  output$footer3 <- renderUI({
    div(textOutput("message3"), style = paste("background-color:", textcolors$pca, "; padding: 10px; border-radius: 20px ;text-align: center;"))
  })
  output$footer4 <- renderUI({
    div(textOutput("message4"), style = paste("background-color:", textcolors$dge, "; padding: 10px; border-radius: 20px ;text-align: center;"))
  })
  output$footer5 <- renderUI({
    div(textOutput("message5"), style = paste("background-color:", textcolors$plots, "; padding: 10px; border-radius: 20px ;text-align: center;"))
  })
}


# Run the application
shinyApp(ui = ui, server = server)
