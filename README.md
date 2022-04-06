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
TRUE)
   This latter command takes longer because the vignette needs to be built.</li>
  </ol>
