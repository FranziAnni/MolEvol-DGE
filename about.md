This is a `shiny` app for performing differential gene expression analyses (DGE). `shiny` is a framework for web apps using `R` code underneath. So by using this app you'll be running `R` code without typing any `R` commands. It is intended as a gentle introduction into differential gene expression analysis that lets you focus on the workflow and the results without having to worry about making your code work. 

**NOTE**

This app is a tool for teaching and learning. It is not intended as a replacement for a thorough RNAseq analysis.

## Manual

The app is organised into several modules that can be accessed by clicking the different tabs in the panel on the very top of this screen. The modules are sorted in the typical order of a DGE from left to right, and would typically be completed in this order as well. You can monitor your progress using the bar on the very bottom of the screen, which will highlight once a module is completed. In the following, a brief explanation of how to use the different modules is given, which should enable you to perform a simple DGE analysis on your own.

### Data upload

Use this module to upload your data. Two files are required:

* **Counts file** 
    
    This file contains the counts per transcript that you determined e.g., with `express`. There should be one transcript per line and one sample per column. The first column must contain the transcript names. All columns must have names, these will be used as sample names in the app. The name for the first column is not important. 
    
    The counts file should be a simple text file with the columns separated by either tabs, commas, or semicolons. You can choose the correct separator using the buttons (this also works after you have selected your file). To load your file, click `Browse...` navigate to your counts file and click `Open`. After upload is complete, you can preview the file you have just uploaded. Change the column separator if necessary. 
    
* **Metadata file**

    This file should summarise your experimental design. This app only supports experiments with a single factor, so this file is very simple: the first column should list all of your samples, and the second column should list the association to the factor. This might be e.g., "Tissue" or "Time point". The table should have headers which however doesn't have to follow any specific format. Again, you can choose the correct column separator using radio buttons. 
    
    After loading your file, make sure it looks as expected. You can now move on to the next module. 
    
### Filtering

It is often desirable to remove some of the loci before continuing with DGE analysis. Many of the transcripts will have very few counts, so they are likely not too important. Filtering also saves computational time in all following steps. In the module, you can filter your data by number of samples in which the locus has mapped reads, and by the number of reads that have mapped to the transcripts. Both can be regulated using the sliders. 

* **Samples**

    Specify here how many samples may be below the read threshold that you specify below. Transcripts for which more samples are below the read count will be removed. If you pick a lower number here, the filtering will be stricter.
    
* **Reads**

    Specify the read threshold here. If you pick a higher number here, the filtering will become stricter. 
    
Example: Let's assume you have 6 samples with 2 treatments (3 samples per treatment). To include only loci for which all of the samples have at least 5 mapped reads, you would specify "0" and "5" for samples and reads, respectively. To keep only the loci for which half of the samples have any mapped reads, you would do "3" and "0". 

The distribution of read counts will be shown for the dataset before and after filtering. in the main panel of this module. You can play around with the filter and see how that impacts your data. Once you are happy with the filter, don't forget to click `Filter now` before you continue. 

### EDA

Exploratory data analysis is used to get an overview of a dataset, especially one with many different variables whose relations are unclear or ones that are too big to allow visual skimming. The main aim here is to identify patterns which may have a biological background. It is important to remember that we are not using any statistical models here, we are simply exploring our data. 

Per default, three EDA plots are created in this module by clicking the button `Plot now`.

* **PCA** 

    **P**rincipal **C**omponent **A**nalysis is a method for reducing dimensions in the data while maintaining variability. Our dataset has many dimensions (as many as rows in the table we uploaded). PCA will take this multidimensional dataset and display it in only two (novel) dimensions that reflect a large part of the variability of the original dataset. The dimensions won't have any meaning, but will help us to understand relationships between the samples. 
    
* **MDS**

    **M**ulti**d**imensional **S**caling is similar to a PCA, however, distances between samples are calculated and used as input into a PCA. Like PCA, a nice approach to reduce dimensionality in the data, and to explore which samples are similar to each other. 
    
* **Heatmap** 

    You can think of a heatmap as our original counts spreadsheet, coloured by the cell content; higher numbers receive lighter colours, and lower numbers darker numbers. The heatmap we use here will however also cluster the data. Clustering means that transcripts with similar expression profiles will be combined into a single 'cell' of the heatmap. A dendrogram shows how these clusters relate to another – similar to a phylogenetic tree. You can choose how many clusters to calculate using the slider. More clusters means a higher resolution heatmap which however may be difficult to read. Low numbers of clusters are easy to read but may obscure the patterns in your data. 
    
You can download the plot using the button at the bottom of the screen. 

### DGE

This module runs the `DeSeq2` algorithm for determining differentially expressed genes using default settings. Using the dropdown menus, select the groups to compare. Not that the first group you select here will act as the 'baseline' in terms of determining if genes are over- or underrepresented in your expression data. After you press the button, calculation will commence. This will take up to a couple of minutes, so please be patient. 

When the calculation is done, you will see 

* **Table of results**

    This displays al of the transcripts, the estimated Log fold change, plus the p-values. Any significantly differentially expressed gene is highlighted in the table. You can adjust the p-value using the slider on the left hand side of the screen.
    
* **MA plot**

     Each point in this plot represents one transcript. On the x-axis, expression levels are displayed, and on the y-acis you can see the log fold change determined by DeSeq2. 
     
* **Histogram of p-values**

    This gives a nice overview of how often you found significantly differentially expressed genes, and how these change in relation to the p-value selected. 

Download buttons let you save the plots, and the list of differentially expressed genes. 

### Gene plots. 

There may be cases in which you want tp plot expression levels for individual loci. To do this, simply type (or copy & paste) in the names of the loci you want to display. Multiple genes can be typed in here, one per line. If the plot becomes to crowded, adjust the size using the slider. 