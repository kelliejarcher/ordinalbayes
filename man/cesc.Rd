% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cesc}
\alias{cesc}
\title{Data Frame of Smaller Subset of The Cancer Genome Atlas Cervical Cancer HTSeq Data.}
\format{
This \code{data.frame} includes the following variables:
\describe{
\item{age_at_index}{a numeric vector storing age at cervical cancer diagnosis.}
\item{cigarettes_per_day}{a numeric vector storing the number of cigarettes smoked per day.}
\item{race}{a character vector representing race of the subject.}
\item{Stage}{an ordered factor representing stage of cervical cancer.}
\item{ENSG00000076344}{a numeric vector representing expression of ENSG00000076344.}
\item{ENSG00000077274}{a numeric vector representing expression of ENSG00000077274.}
\item{ENSG00000101888}{a numeric vector representing expression of ENSG00000101888.}
\item{ENSG00000115548}{a numeric vector representing expression of ENSG00000115548.}
\item{ENSG00000122884}{a numeric vector representing expression of ENSG00000122884.}
\item{ENSG00000125430}{a numeric vector representing expression of ENSG00000125430.}
\item{ENSG00000131370}{a numeric vector representing expression of ENSG00000131370.}
\item{ENSG00000135443}{a numeric vector representing expression of ENSG00000135443.}
\item{ENSG00000136457}{a numeric vector representing expression of ENSG00000136457.}
\item{ENSG00000138398}{a numeric vector representing expression of ENSG00000138398.}
\item{ENSG00000150636}{a numeric vector representing expression of ENSG00000150636.}
\item{ENSG00000161277}{a numeric vector representing expression of ENSG00000161277.}
\item{ENSG00000163510}{a numeric vector representing expression of ENSG00000163510.}
\item{ENSG00000164485}{a numeric vector representing expression of ENSG00000164485.}
\item{ENSG00000164651}{a numeric vector representing expression of ENSG00000164651.}
\item{ENSG00000166091}{a numeric vector representing expression of ENSG00000166091.}
\item{ENSG00000166342}{a numeric vector representing expression of ENSG00000166342.}
\item{ENSG00000171121}{a numeric vector representing expression of ENSG00000171121.}
\item{ENSG00000177173}{a numeric vector representing expression of ENSG00000177173.}
\item{ENSG00000180229}{a numeric vector representing expression of ENSG00000180229.}
\item{ENSG00000188817}{a numeric vector representing expression of ENSG00000188817.}
\item{ENSG00000197360}{a numeric vector representing expression of ENSG00000197360.}
\item{ENSG00000203601}{a numeric vector representing expression of ENSG00000203601.}
\item{ENSG00000225449}{a numeric vector representing expression of ENSG00000225449.}
\item{ENSG00000230201}{a numeric vector representing expression of ENSG00000230201.}
\item{ENSG00000233996}{a numeric vector representing expression of ENSG00000233996.}
\item{ENSG00000236138}{a numeric vector representing expression of ENSG00000236138.}
\item{ENSG00000236819}{a numeric vector representing expression of ENSG00000236819.}
\item{ENSG00000250602}{a numeric vector representing expression of ENSG00000250602.}
\item{ENSG00000253923}{a numeric vector representing expression of ENSG00000253923.}
\item{ENSG00000256980}{a numeric vector representing expression of ENSG00000256980.}
\item{ENSG00000259083}{a numeric vector representing expression of ENSG00000259083.}
\item{ENSG00000259134}{a numeric vector representing expression of ENSG00000259134.}
\item{ENSG00000260484}{a numeric vector representing expression of ENSG00000260484.}
\item{ENSG00000263612}{a numeric vector representing expression of ENSG00000263612.}
\item{ENSG00000264049}{a numeric vector representing expression of ENSG00000264049.}
\item{ENSG00000264954}{a numeric vector representing expression of ENSG00000264954.}
\item{ENSG00000265579}{a numeric vector representing expression of ENSG00000265579.}
\item{ENSG00000271711}{a numeric vector representing expression of ENSG00000271711.}
\item{ENSG00000272071}{a numeric vector representing expression of ENSG00000272071.}
\item{ENSG00000276517}{a numeric vector representing expression of ENSG00000276517.}
}
}
\source{
TCGAbiolinks TCGA-CESC
}
\usage{
cesc
}
\description{
This data.frame is a subset that includes 41 transcripts and four phenotypic variables for
242 subjects with squamous cell cervical cancer.
}
\examples{
data("cesc")
head(cesc)
}
\keyword{datasets}
