#Script: Gen2Met.dataPreparation.R
#License: GPLv3 or later
#Modification date: 2023-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Prepare the dataset for Gen2Met core function. This function converts gene ID from ensembl to entrez (or vice versa) and to add the gene symbol.

Gen2Met.dataPreparation <- function(
  in.file,
  biomart.species=NULL
  )
{
  #check the arguments genelist, biomartspecies, gene_id
  if (file.exists(in.file)){
    genelist <- read.table(in.file, header=F, stringsAsFactors = F)
    cat("Input file imported! \n")
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }

  if (is.null(biomart.species)){
    cat("\n")
    cat("Biomart species not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{
    orgs <- biomaRt::listDatasets(useMart("ensembl"))
    if (biomart.species %in% orgs$dataset){
      cat("\n")
      cat("BiomaRt species correct! \n")
    }else{
      cat("\n")
      cat("BiomaRt species incorrect. Please check! \n")
      stop("Exit",call. = F)
    }
  }

  #generate output file name
  in.file <- gsub("\\..*", "", in.file)

  #convert to gene id in animal organism and export the result
  cat("\n")
  cat("Convertion from gene symbol to entrez ID ... \n")
  genelist2 <- biomaRt::getBM(attributes = c("entrezgene_id", "external_gene_name"),
                              filters = c("external_gene_name"),
                              values = genelist,
                              mart=(biomaRt::useMart("ensembl", biomart.species)), verbose = FALSE)
  colnames(genelist2) <- c("entrezgene_id", "external_gene_name")
  genelist2 <- genelist2[!is.na(genelist2$entrezgene_id),]
  cat("DONE\n")
  #cat("\n")
  #cat("Checking for corresponding genes in KEGG database ... \n")
  cat("n.", paste(sum(length(which(!is.na(genelist2$entrezgene_id))))), "out of", paste(sum(length(which(!is.na(genelist$V1))))), "genes have a corresponding entrez gene ID in BioMart database.", sep= " ")
  wd <- getwd()
  write.table(x = genelist2, file = paste(in.file,"_converted.txt",sep=""), row.names = F)
  cat("\n")
  cat("Gene list exported!")
  return(genelist2)
}
