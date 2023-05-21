#Script: Gen2Met.checkPathfile.R
#License: GPLv3 or later
#Modification date: 2023-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Check the genelist before run main PANEV script

Gen2Met.checkPathfile <- function(
  path.file=NULL
)
{
  #check the arguments in the genelist
  if (is.null(path.file)){
    cat("Path file not specified. All KEGG pathways will be considered in the final network! \n")
  }else{
    cat("Path file specified...")
    options(warn=-1)
    if ("external_gene_name" %in% colnames(genelist) & "entrezgene_id" %in% colnames(genelist)){
      cat("and correct! \n")
    }else{
      cat("but is incorrect. Please check! \n")
      cat("Remember to use the output file of Gen2Met.enrichment() function \n")
      stop("Exit",call. = F)
    }
  }
}
