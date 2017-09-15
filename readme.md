# Varying-Censoring Aware Matrix Factorization (VAMF)
VAMF is a probabilistic dimension reduction method for single cell RNA-Seq datasets. For details see the [biorxiv preprint](http://www.biorxiv.org/content/early/2017/07/21/166736). VAMF is also implemented as an [R package](https://github.com/willtownes/vamf).

## Repository Contents

* **algs**- implementations of the VAMF algorithm as well as wrapper functions for competing methods.
* **final_figs**- intermediate data and code for generating the figures used in the paper itself.
* **patel_2014**, **shalek_2013**, **trapnell_2014**- code for processing the data used in the three studies. For the former, see also the [data package](https://github.com/willtownes/patel2014gliohuman)
* **simulations**- code for the *noise only* and *latent clusters* scenarios. Start with `alg_tests.Rmd` then `sensitivity.Rmd`



## Contact

* Primary developer: [Will Townes](http://willtownes.github.io)
* Collaborators: [Martin Aryee](http://www.massgeneral.org/pathology/research/researchlab.aspx?id=1555), [Stephanie Hicks](http://www.stephaniehicks.com/), [Rafael Irizarry](http://rafalab.dfci.harvard.edu/)



