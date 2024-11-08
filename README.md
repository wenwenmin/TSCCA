# TSCCA

### Introduction

In this study, we present a Tensor Sparse Canonical Correlation Analysis (TSCCA) method for identifying cancer-related miRNA-gene modules across multiple cancers. TSCCA is able to overcome the drawbacks of existing solutions and capture both the cancer-shared and specific miRNA-gene co-expressed modules with better biological interpretations. We comprehensively evaluate the performance of TSCCA using a set of simulated data and matched miRNA/gene expression data across 33 cancer types from the TCGA database.

This package is used to solve TSCCA model in our paper <a class="footnote-reference" href="#id2" id="id1">[1]</a>. 
More descriptions about these functions can be found in their annotation part. 


<p align="center"> 
<img src="https://github.com/wenwenmin/TSCCA/blob/master/Figures/Fig1_TSCCA.png">
</p>
Figure 1. Illustration of TSCCA to identify cancer-miRNA-gene functional modules.

### R code
A toy example explains how to use the TSCCA function. Before running the script, please first set the path for "script1_toy_example.R",
and then run the following R command in the Console. 

``` r
> source('script1_toy_example.R') 
```

Next, we applied TSCCA and SCCA to these simulated data in the second example. Please run the following R command in the Console.
``` r
> source('script2_simulation.R') 
```

### TCGA data 
In addition, we also comprehensively evaluated the performance of TSCCA using a set of the matched miRNA/gene expression data across 33 cancer types from the TCGA database (please see the folder "Results_tcga").

### References
<table class="docutils footnote" frame="void" id="id2" rules="none">
<colgroup><col class="label" /><col /></colgroup>
<tbody valign="top">
<tr><td class="label"><a class="fn-backref" href="#id2">[1]</a></td><td> 
Wenwen Min, Tsung-Hui Chang, Shihua Zhang* and Xiang Wan*. TSCCA: A tensor sparse CCA method for detecting microRNA-gene patterns from multiple cancers. PLOS Computational Biology 2021, DOI: 10.1371/journal.pcbi.1009044. 
</td></tr>
</tbody>
</table>
