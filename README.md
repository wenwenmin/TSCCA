# TSCCA

### Introduction

This package is used to solve Tensor Sparse Canonical Correlation Analysis (TSCCA) model in our paper <a class="footnote-reference" href="#id2" id="id1">[1]</a>. 

More descriptions about these functions can be found in their annotation part.

<p align="center"> 
<img src="https://github.com/wenwenmin/TSCCA/blob/master/Figures/TSCCA.png">
</p>

### R code
Note that before running the codes, please first set the path for "Fun2_simulation_study_main.R".
Please run the following R command in the Console. 

``` r
> source('Fun2_simulation_study_main.R') 
```

### TCGA data 
We comprehensively evaluate the performance of TSCCA using a set of simulated data and matched miRNA/gene expression data across 33 cancer types from the TCGA database.
<p align="center"> 
<img src="https://github.com/wenwenmin/TSCCA/blob/master/Fig_tcga_data_table.png">
</p>


### References
<table class="docutils footnote" frame="void" id="id2" rules="none">
<colgroup><col class="label" /><col /></colgroup>
<tbody valign="top">
<tr><td class="label"><a class="fn-backref" href="#id2">[1]</a></td><td> 
Wenwen Min, Tsung-Hui Chang, Shihua Zhang and Xiang Wan. TSCCA: A tensor sparse CCA method for detecting microRNA-gene patterns from multiple cancers. PLOS Computational Biology (Revision). 
</td></tr>
</tbody>
</table>
