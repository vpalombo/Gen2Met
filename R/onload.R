.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: Gen2Met\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using Gen2Met v1.0!\n")
  #mymsg <- paste(mymsg,"For more information use: help(package = 'Gen2Met')\n")
  #mymsg <- paste(mymsg,"                          browseVignettes(package = 'Gen2Met')\n\n")
  mymsg <- paste(mymsg,"Citation:\n")
  mymsg <- paste(mymsg,"  Author: V. Palombo valentino.palombo@gmai.com\n")
  mymsg <- paste(mymsg,"  Title: Gen2Met: an R package in multi-omic Era\n")
  #mymsg <- paste(mymsg,"  Journal: \n")
  #mymsg <- paste(mymsg,"  \n")
  packageStartupMessage(mymsg)

  suppressWarnings(suppressPackageStartupMessages(library(bc3net)))
  suppressWarnings(suppressPackageStartupMessages(library(data.table)))
  suppressWarnings(suppressPackageStartupMessages(library(stringr)))
  suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
  suppressWarnings(suppressPackageStartupMessages(library(XML)))
  suppressWarnings(suppressPackageStartupMessages(library(xml2)))
  suppressWarnings(suppressPackageStartupMessages(library(visNetwork)))
  suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
  suppressWarnings(suppressPackageStartupMessages(library(BiocManager)))
  suppressWarnings(suppressPackageStartupMessages(library(KEGGREST)))
  suppressWarnings(suppressPackageStartupMessages(library(biomaRt)))

  require(bc3net)
  require(data.table)
  require(stringr)
  require(dplyr)
  require(XML)
  require(xml2)
  require(visNetwork)
  require(RColorBrewer)
  require(BiocManager)
  require(KEGGREST)
  require(biomaRt)

}

