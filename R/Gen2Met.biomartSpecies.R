#Script: Gen2Met.biomartSpecies.R
#License: GPLv3 or later
#Modification date: 2023-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Find the organism code among the list of those available on biomaRt.

Gen2Met.biomartSpecies <- function(
  string = NULL
  )
{
  #download the biomaRt ensembl dataset and show all the list of organisms available
  mart <- biomaRt::useMart("ensembl")
  orgs <- biomaRt::listDatasets(mart)
  colnames(orgs) <- c("organism_code", "description", "version")

  #check the argument
  if (is.null(string)){
    cat("The list of all available species for gene id convertion and data preparation was created! \n")
    cat("Remember to use the correct organism code for relative Gen2Met functions. \n")
    return(orgs)
  }else{
    #download the biomaRt ensembl dataset and show only the organism(s) matching the string of interest
    org <- grep(string, orgs$description, ignore.case = T, value=T)
    if (length(org)==0){
      cat("No match found for the requested string! \n")
    }else{
      orgs <- data.frame(orgs[orgs$description %in% org,1:3], row.names = NULL)
      colnames(orgs) <- c("organism_code", "description", "version")
      cat("The list of available species matched your string was created! \n")
      cat("Remember to use the correct organism code for relative Gen2Met functions. \n")
      return(orgs)
    }
  }
}
