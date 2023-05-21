#Script: Gen2Met.checkGenelist.R
#License: GPLv3 or later
#Modification date: 2023-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Check the genelist before run main PANEV script

Gen2Met.checkGenelist <- function(
  genelist=NULL
)
{
  #check the arguments in the genelist
  if (is.null(genelist)){
    cat("Gene list not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{
    cat("Gene list specified... ")
  }

  options(warn=-1)
  if ("external_gene_name" %in% colnames(genelist) & "entrezgene_id" %in% colnames(genelist)){
    cat("and correct! \n")
  }else{
    cat("but is incorrect. Please check! \n")
    cat("Remember your gene list must have two columns with 'entrez gene ID'  and 'gene symbol', named 'entrezgene_id' and 'external_gene_name' respectively \n")
    stop("Exit",call. = F)
  }
}
