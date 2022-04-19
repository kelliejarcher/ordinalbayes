# ordinalbayes

## Installing ordinalbayes

1. Install JAGS https://mcmc-jags.sourceforge.io
2. Install R https://cran.r-project.org
3. Install runjags which is hosted on CRAN and 
http://runjags.sourceforge.net
4. Install Bioconductor packages 
<ol type="a">
  <li>See https://www.bioconductor.org/install/ to install base Bioconductor packages</li>
 <li>To run examples using finalSet and reducedSet the 
     DESeq2 package is needed. To install run</li>
     > BiocManager::install("DESeq2")</li>
</ol>
5. Install ordinalbayes using 
<ol type="a">
  <li>Within R use > install.packages("ordinalbayes")</li>
      or
 <li> > devtools::install_github("kelliejarcher/ordinalbayes", build_vignettes = 
TRUE)</li>
   This latter command takes longer because the vignette needs to be built.
  </ol>

## References
1. Archer KJ, Seffernick AE, Sun, S, Zhang, Y. ordinalbayes: Fitting ordinal Bayesian regression models to high-dimensional data using R. <i>STATS</i>, 5(2), 371-384, 2022. <a href>https://doi.org/10.3390/stats5020021</a>
2. Zhang Y, Archer KJ. Bayesian variable selection for high-dimensional data with an ordinal response: identifying genes associated with prognostic risk group in acute myeloid leukemia. <i>BMC Bioinformatics</i>, 22:539, 2021. <a href>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8565083/</a>
3. Zhang Y, Archer KJ. Bayesian penalized cumulative logit model for high-dimensional data with an ordinal response. <i>Statistics in Medicine</ii>, 40:1453-1481, 2021. <a href>https://pubmed.ncbi.nlm.nih.gov/33336826/</a>
