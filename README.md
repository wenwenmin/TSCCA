# TSCCA

### Introduction

This package is used to solve Tensor Sparse Canonical Correlation Analysis (TSCCA) model in our paper <a class="footnote-reference" href="#id2" id="id1">[1]</a>. 
More descriptions about these functions can be found in their annotation part. 

#### TSCCA flowchart
The following flowchart illustrates how TSCCA identifies cancer-miRNA-gene functional modules. 
<p align="center"> 
<img src="https://github.com/wenwenmin/TSCCA/blob/master/Figures/Fig1_TSCCA.png">
</p>

### R code
Note that before running the codes, please first set the path for "script1_toy_example.R".
Please run the following R command in the Console and this simple example explains how to use the TSCCA function. 

``` r
> source('script1_toy_example.R') 
```

Next, we applied TSCCA and SCCA to these simulated data in the second example. Please run the following R command in the Console.
``` r
> source('script2_simulation.R') 
```

### TCGA data 
In addition, we also comprehensively evaluated the performance of TSCCA using a set of simulated data and matched miRNA/gene expression data across 33 cancer types from the TCGA database (please see the folder "Results_tcga").

### References
<table class="docutils footnote" frame="void" id="id2" rules="none">
<colgroup><col class="label" /><col /></colgroup>
<tbody valign="top">
<tr><td class="label"><a class="fn-backref" href="#id2">[1]</a></td><td> 
Wenwen Min, Tsung-Hui Chang, Shihua Zhang and Xiang Wan. TSCCA: A tensor sparse CCA method for detecting microRNA-gene patterns from multiple cancers. PLOS Computational Biology (Revision). 
</td></tr>
</tbody>
</table>
