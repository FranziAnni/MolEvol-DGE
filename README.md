# MolEvol-DGE
A simple DESeq2 wrapper for teaching purposes

To install, open `R` and install required packages:

```R
install.packages(c(
  "shiny",
  "DT",
  "factoextra",
  "tidyverse",
  "patchwork",
  "purrr",
  "shinycssloaders",
  "shinythemes",
  "ggrepel",
  "ggpubr",
  "pheatmap",
  "ggplotify"), dependencies = TRUE)
  
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

To run the app, type 

```R
shiny::runUrl("https://github.com/FranziAnni/MolEvol-DGE/archive/refs/heads/main.zip")
```

