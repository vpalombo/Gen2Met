#Script: Gen2Met.speciesCode
#License: GPLv3 or later
#Modification date: 2023-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Show the list of species available for Gen2Met

Gen2Met.speciesCode <- function(
  string = NULL
  )
{
  #download the list of KEGG available organisms
  orgs <- as.data.frame(KEGGREST::keggList("organism"))[2:3]
  #orgs <- as.data.frame(cbind(paste(orgs$organism), paste(orgs$species)))
  colnames(orgs) <- c("Gen2Met_code", "species" )
  #check for a specific species of interest
  if (is.null(string)){
    cat("The list of all species available for Gen2Met analysis was created! \n")
    cat("Remember to use the correct organism code for relative Gen2Met functions. \n")
    return(orgs)
  }else{
    #look for the specified string
    species <- grep(string, orgs$species, ignore.case = T, value=T)
    if (length(species) >= 1){
      Gen2Met_code <- paste(orgs[orgs$species %in% species,1])
      result <- as.data.frame(cbind(species,Gen2Met_code))
      cat("The list of available species, matching your string, was created! \n")
      cat("Remember to use the correct organism code for relative Gen2Met functions. \n")
      return(result)
    }else{
      cat("No match found for the requested string! \n")
    }
  }
}
