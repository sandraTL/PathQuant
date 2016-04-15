# PathQuant (Pathway Quantify)

---------------

R package to *quantify* gene-metabolite associations (e.g. from mGWAS) to 
[KEGG](http://www.genome.jp/kegg/)'s pathways.

*Quantify*: calculate distance defined as the shortest path between a given gene
and a metabolite in pathway.

---------------

### Installation
1.
```r
--- code to run package ---
# If you don not have devtools installed
install.packages("devtools")

# Install PathQuant 
library(devtools)
devtools::install_github("SandraTL/PathQuant")
library("PathQuant")
```

### Information

License: [GNU General Public License (v3)](http://www.gnu.org/licenses/gpl-3.0.en.html)

### Instructions

* [Userguide Manuel](https://github.com/sandraTL/PathQuant/blob/master/manual.pdf)
* [Vignette](https://github.com/sandraTL/PathQuant/blob/master/vignettes/vignette.Rmd)

### TODO

* Add functionnalities for multiple KEGG map distance calculations
 


