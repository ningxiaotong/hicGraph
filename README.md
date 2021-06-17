# hicGraph
hicGraph is an R package for visualizing Hi-C data and ChIP-Seq data.<br/>
So far, the visualizations that hicGraph can achieve include **heatmap** and **network** graph for Hi-C data, and **track** graph of ChIP-Seq data.<br/>
The realization of hicGraph function is based on the R language packages: **ggplot2**、**rtracklayer**、**visNetwork**、**shiny**.<br/>
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/) used to plot heatmap in hicGraph
- [rtrackalyer](http://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) used to plot track graph in hicGraph
- [visNetwork](https://cran.r-project.org/web/packages/visNetwork/) used to plot network in hicGraph
- [shiny](https://shiny.rstudio.com) used to intergrate heatmap、network and track hicgraph into a web page structure.



# Requirements
- [R](https://www.r-project.org) or [Rstudio](https://www.rstudio.com)


# Install
```
devtools::install_github('ningxiaotong/hicGraph')
```

Or you can click [here](./hicGraph_1.0.tar.gz) to download the .tar.gz to the local, and use the following method to install hicGraph.
```
install.packages("path-to/hicGraph_1.0.tar.gz", repos = NULL, type = "source")
```


# Usage
- Library packages
```
library(hicGraph)
library(visNetwork)
library(shiny)
```

- Run app
```
run_app() 
```

# Inputs
- Input the chr_num. Eg: chr18 (Please note that you must use this format for chr_num.)
- Input the data file:
  - Interaction matrix of 80k resolution. # format: `.csv`
  - Interaction matrix of 20k resolution. # format: `.csv`
  - A file containing information about the interaction between bins. # format: `.mcool`
  - A File containing information about the location of the compartments. # format: `.bed`
  - A File containing information about the location of the tads. # format: `.bed`
  - A File containing ChIP-Seq data. # format: `.bigwig`


# Test data
Please click [here](https://pan.baidu.com/s/1i5XEsVbNe-1YUh3K512IzQ) to download the test data.
The extraction code is `sf2g`.

# Demonstration of operation
You can click [here](https://pan.baidu.com/s/1ioQlEU04vK3cifSkiYARzQ) to view a demo of the operation.
The extraction code is `8ca9`.


# Note
When you Input your all input data, hicGraph needs about 4-5 hours to plot the graphs. Please wait paiently.
