#Script: Gen2Met.enrichment.R
#License: GPLv3 or later
#Modification date: 2023-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: KEGG Pathway Enrichment analysis

Gen2Met.enrichment <- function(
  in.file=NULL,
  species=NULL,
  FDR=1,
  out.file="Gen2Met"
  )
{
  #cat("Arguments checking ... \n")
  #import genelist
  in.file.name <- paste(in.file)
  if (file.exists(in.file.name)){
    genelist <<- read.table(in.file.name, header=T, stringsAsFactors = F)
    cat("Input file is imported! \n")
    cat("\n")
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  #check genelist
  Gen2Met.checkGenelist(genelist = genelist)
  #check specie code
  Gen2Met.checkSpecies(species)

  # Function #
  cat("\nEnrichment analysis started ... ")
  #download genes list available with the relative root pathway(s)
  KEGGgenes <- data.frame(cbind(rn=names(KEGGREST::keggLink(paste(species),"pathway")),
                                V1=KEGGREST::keggLink(paste(species),"pathway")), stringsAsFactors = F, row.names = NULL)
  genelist$KEGGgenes <- paste(paste(species),":",genelist$entrezgene, sep="")
  genelist <- merge(x = genelist, y = KEGGgenes, by.x = "KEGGgenes", by.y = "V1", all=F)
  #download the available KEGG pathways list
  KEGGpath <- data.table::setDT(as.data.frame(cbind(KEGGREST::keggList("pathway", species))), keep.rownames = T)
  orgs <- as.data.frame(KEGGREST::keggList("organism"))[2:3]
  string <- paste(" -", paste(orgs[which(orgs$organism == species),2]))
  KEGGpath$V1 <- stringr::str_replace_all(KEGGpath$V1, stringr::fixed(string), "")
  #perform and enrichment analysis
  genes <- as.vector(unique(genelist$KEGGgenes))
  KEGGgenes2 <- as.data.frame(KEGGgenes)
  KEGGgenes2$V1 <- as.character(KEGGgenes2$V1)
  reference <- unique(paste(KEGGgenes2$V1))
  pathtofind <- unique(KEGGgenes2$rn)

  #remove general pathways
  path_to_remove <- c("map01100", "map01110", "map01120") #because too general pathways
  path_to_remove <- stringr::str_replace(path_to_remove, "map", paste("path:",species, sep = ""))
  pathtofind <- pathtofind[!(pathtofind %in% path_to_remove)]

  genesets <- list()
  KEGGpath$rn <- stringr::str_replace(KEGGpath$rn, paste(species), paste("path:",species, sep = ""))
  for (i in 1:length(pathtofind)) {
    genesets[[paste(pathtofind[i])]] <- subset(KEGGgenes2$V1, KEGGgenes2$rn==pathtofind[i])
  }
  enrich <- bc3net::enrichment(genes, reference, genesets, adj = "fdr", verbose = FALSE)
  enrich <- merge(x = enrich, y = KEGGpath, by.x = "TermID", by.y = "rn", all=F)
  colnames(enrich) <- c("pathway_ID", "n_genes", "all_genes", "pvalue", "FDR", "pathway_name")
  enrich <- enrich[order(enrich$FDR, enrich$pvalue),]
  enrich <- enrich[enrich$FDR <= FDR,]
  write.table(x = enrich, file = paste(out.file,"_enrichment.txt",sep=""), row.names = F)
  cat("\n")
  cat("and results exported! \n")
  #create and save the list of genes encompassed within the over-represented pathway(s)
  pathtofinterest <- enrich$pathway_ID
  genelist2 <- genelist[(genelist$rn %in% pathtofinterest),]
  genelist2 <- data.frame(genelist2$entrezgene_id, genelist2$external_gene_name)
  colnames(genelist2) <- c("entrezgene_id", "external_gene_name")
  genelist2 <- unique(genelist2)
  write.table(x = genelist2, file = paste(out.file,"_geneList.txt",sep=""), row.names = F)
  cat("\n")
  cat("The list of genes encompassed within enriched pathway(s) was created and exported! \n")
  #create and save the pathway(s) per gene(s) frequency table
}
