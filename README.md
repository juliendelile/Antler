<div><img src="https://raw.githubusercontent.com/juliendelile/Antler/master/logo.png" width="300px" align="right"></div>

# ANTLER

**Antler** (**An**other **t**ranscriptome **l**ineage **e**xplore**r**) is an R package providing a set of methods for analysing single-cell RNA-seq experiments.

The idea behind **Antler** is to perform completely unbiased and data-driven analysis to enable the discovery of novel candidate genes and transcriptomic states.

Two vignettes demonstrate how **Antler** can be used to:

1. Identify, visualize and export the transcriptomic states of a dataset ([Link](https://juliendelile.github.io/Antler/articles/Transcriptomic-summary.html)).

2. Build the lineage tree of a cell population differentiating over time and infer the pseudotime dynamics of gene expression ([Link](https://juliendelile.github.io/Antler/articles/Carving-cell-trajectories.html)).

### Installation

**Antler** can be installed from Github via the [devtools](https://github.com/hadley/devtools) package:

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("juliendelile/Antler", dependencies = TRUE)
```
### Examples

Some of the plots produced in the two vignettes include:

* the transcriptomic summary of the cell states

<p align="center"><img src="https://raw.githubusercontent.com/juliendelile/Antler/master/vignettes/images/Transcriptome_summary_selection_I_hclust_selection_I_Normalized_logscaled.png" width="90%"></p>

* the cell state graph enabling the pseudotime ordering of the cells.

<p align="center"><img src="https://raw.githubusercontent.com/juliendelile/Antler/master/vignettes/images/PT_tree_cell_state_graph_landmarks.png" width="50%"></p>

* the lineage tree 

<p align="center"><img src="https://raw.githubusercontent.com/juliendelile/Antler/master/vignettes/images/PT_tree_communities.png" width="40%"></p>

* the reconstructed pseudotime dynamics of a gene

<p align="center"><img src="https://raw.githubusercontent.com/juliendelile/Antler/master/vignettes/images/PT_dynamics_lin09g01_Normalized.png" width="80%"></p>