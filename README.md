# MolEvol-DGE
A simple DESeq2 wrapper for teaching purposes

To install, open `R` and install required packages:

```
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
  "gridExtra"), dependencies = TRUE)
  
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

To run the app, type 

```
shiny::runGitHub("FranziAnni/MolEvol-DGE")`
```

