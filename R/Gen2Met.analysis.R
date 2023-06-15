#Script: Gen2Met.analysis.R
#License: GPLv3 or later
#Modification date: 2023-10-15
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Run Gen2Met core analysis.

Gen2Met.analysis <- function(
  in.file=NULL,
  out.file="Gen2Met_analysis",
  species=NULL
)
{
  #import genelist
  in.file.name <- paste(in.file)
  if (file.exists(in.file.name)){
    genelist <<- read.table(in.file.name, header=T, stringsAsFactors = F)
    cat("Input file imported! \n")
    cat("\n")
    #check the genelist arguments
    Gen2Met.checkGenelist(genelist = genelist)
    #delete gene(s) with no entrez gene id
    genelist <- genelist[!is.na(genelist$entrezgene_id),]
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  #check specie code
  Gen2Met.checkSpecies(species)
  cat("\n")
  cat("\nPrerequisite check passed!\n\nGen2Met is running ... \n")
  cat(" Please wait... It could take a while depending on the number of genes investigated! \n")

  genes <- paste(species, ":", genelist$entrezgene_id, sep="")
  cat("\n\nGen2Met is searching for modules and reactions in association with a total of", paste(length(genes), "genes  ... \n"))

  #lavoriamo sui geni
  path2genes <- NULL
  mod2genes <- NULL
  gene2Reaction <- NULL
  for (i in 1:length(genes)){
    out_nodes <- NULL
    print(i)
    print(Sys.time())
    ort <- names(KEGGREST::keggGet(genes[i])[[1]][["ORTHOLOGY"]]) #convertito in ortologo
    name <- KEGGREST::keggGet(ort)[[1]]$NAME
    #name <- keggGet(genes[i])[[1]]$SYMBOL
    if(is.null(name)){
      name <- genes[i]
    }
    dataDownloaded <- KEGGREST::keggGet(ort)[[1]]
    out_nodes <- data.frame(to = dataDownloaded$PATHWAY)
    if (dim(out_nodes)[1] != 0){
      out_nodes$id_to <- rownames(out_nodes)
      if(nrow(out_nodes)>0){
        rownames(out_nodes) <- NULL
        out_nodes$id_from <- genes[i]
        out_nodes$from <- name
        path2genes <- rbind(path2genes, out_nodes)
      }
    }
    out_nodes <- NULL
    out_nodes <- data.frame(to = dataDownloaded$MODULE)
    if (dim(out_nodes)[1] != 0){
      if(nrow(out_nodes)>0){
        out_nodes$id_to <- rownames(out_nodes)
        rownames(out_nodes) <- NULL
        out_nodes$id_from <- genes[i]
        out_nodes$from <- name
        mod2genes <- rbind(mod2genes, out_nodes)
      }
    }
    bridge <- data.frame(to = dataDownloaded$DBLINKS)
    if (dim(bridge)[1] != 0){
      subsetted <- bridge[grepl("RN: *", bridge$to),]
      if (length(subsetted) != 0){
        subsetted <- stringr::str_remove_all(subsetted, "RN: ")
        subsetted <- scan(text = subsetted, what = " ")
        if(length(subsetted)>0){
          if(length(subsetted)>1){
            exp_file <- data.frame(id_to = genes[i], id_from = do.call(rbind, strsplit(subsetted, "-", fixed=TRUE)))
            exp_file$to <- name
            for (k in 1:length(exp_file$id_to)){
              exp_file$from[k] <- KEGGREST::keggGet(exp_file$id_from[k])[[1]]$DEFINITION
            }
          }else{
            exp_file <- data.frame(id_to=genes[i], id_from=subsetted)
            exp_file$to <- name
            exp_file$from <- KEGGREST::keggGet(exp_file$id_from)[[1]]$DEFINITION
          }
          if(nrow(exp_file)>0){
            gene2Reaction <- rbind(gene2Reaction, exp_file)
          }
        }
      }
    }
  }
  path2genes <- unique(path2genes)
  mod2genes <- unique(mod2genes)
  gene2Reaction <- unique(gene2Reaction)
  #reorder columns based on common order used before
  gene2Reaction <- gene2Reaction[,colnames(path2genes)]
  path2mod <- NULL
  comp2mod <- NULL
  modules <- unique(mod2genes$id_to)
  cat("\nA total of", paste(length(modules), "modules were associated with your genelist based on information retrieved from KEGG!", "\n\nGen2Met is searching for compounds in association with these modules ... \n"))
  cat(" Please wait... It could take a while depending on the number of modules investigated! \n")
  #modules <- str_replace_all(modules, pattern = "vvi_", "")
  #query on KEGG
  if (!is.null(modules)){
    for (i in 1:length(modules)){
      out_nodes <- NULL
      print(i)
      print(Sys.time())
      name <- KEGGREST::keggGet(modules[i])[[1]]$NAME
      dataDownloaded <- KEGGREST::keggGet(modules[i])[[1]]
      out_nodes <- data.frame(to = dataDownloaded$PATHWAY)
      if (dim(out_nodes)[1] != 0){
        out_nodes$id_to <- rownames(out_nodes)
        if(nrow(out_nodes)>0){
          rownames(out_nodes) <- NULL
          out_nodes$id_from <- modules[i]
          out_nodes$from <- name
          path2mod <- rbind(path2mod, out_nodes)
        }
        out_nodes <- NULL
        out_nodes <- data.frame(to = dataDownloaded$COMPOUND)
        if(nrow(out_nodes)>0){
          out_nodes$id_to <- rownames(out_nodes)
          rownames(out_nodes) <- NULL
          out_nodes$id_from <- modules[i]
          out_nodes$from <- name
          comp2mod <- rbind(comp2mod, out_nodes)
        }
      }
    }
  }
  path2mod <- unique(path2mod)
  comp2mod <- unique(comp2mod)
  #lavoriamo sulle reazioni
  reaction <- unique(gene2Reaction$id_from)
  cat("\nA total of", paste(length(reaction), "reactions were associated with your genelist based on information retrieved from KEGG!", "\n\nGen2Met is searching for compounds in association with these reactions ... \n"))
  cat(" Please wait... It could take a while depending on the number of reactions investigated! \n")
  comp2Reaction <- NULL
  if (!is.null(reaction)){
    for (i in 1:length(reaction)){
      print(i)
      print(Sys.time())
      name <- KEGGREST::keggGet(reaction[i])[[1]]$DEFINITION
      dataDownloaded <- KEGGREST::keggGet(reaction[i])[[1]]
      bridge <- data.frame(to = dataDownloaded$EQUATION)
      if (dim(bridge)[1] != 0){
        subsetted <- bridge[grepl("C*", bridge)]
        subsetted <- stringr::str_remove_all(subsetted, c("[+]"))
        subsetted <- stringr::str_remove_all(subsetted, c("<=>"))
        subsetted <- stringr::str_remove_all(subsetted, c("[(n)]"))
        subsetted <- stringr::str_remove_all(subsetted, c("-1"))
        subsetted <- scan(text = subsetted, what = " ")
        subsetted <- subsetted[grepl("C", subsetted)]
        subsetted <- unique(subsetted)
        subsetted <- subsetted[!(subsetted %in% c("C013551", "C035411", "C007182", "C02128m", "C00039m"))] #compounds with no entry in KEGG
        if(length(subsetted)>0){
          if(length(subsetted)>1){
            exp_file <- data.frame(id_to = reaction[i], id_from = do.call(rbind, strsplit(subsetted, "-", fixed=TRUE)))
            exp_file$to <- name
            for (k in 1:length(exp_file$id_to)){
              exp_file$from[k] <- paste(KEGGREST::keggGet(exp_file$id_from[k])[[1]]$NAME, collapse = " ")
            }
          }else{
            exp_file <- data.frame(id_to=reaction[i], id_from=subsetted)
            exp_file$to <- name
            exp_file$from <- paste(KEGGREST::keggGet(exp_file$id_from)[[1]]$NAME, collapse = " ")
          }
          if(nrow(exp_file)>0){
            comp2Reaction <- rbind(comp2Reaction, exp_file)
          }
        }
      }
    }
  }
  comp2Reaction <- unique(comp2Reaction)
  #reorder columns based on common order used before
  comp2Reaction <- comp2Reaction[,colnames(path2genes)]

  mod2genes$interaction <- "gene2mod"
  comp2mod$interaction <- "mod2comp"
  gene2Reaction$interaction <- "reaction2gene"
  comp2Reaction$interaction <- "comp2reaction"

  data <- data.frame(rbind(mod2genes,
                           comp2mod,
                           gene2Reaction,
                           comp2Reaction
  ), stringsAsFactors = T)

  label <- data.frame(name=data$id_to, shared=data$to)
  label <- unique(label)
  label2 <- data.frame(name=data$id_from, shared=data$from)
  label2 <- unique(label2)
  label_total <- rbind(label, label2)
  label_total <- unique(label_total)

  library(data.table)
  label_total$obj <- "gene"
  label_total$obj[label_total$name %like% "M"] <- "module"
  label_total$obj[label_total$name %like% "map*"] <- "path"
  label_total$obj[label_total$name %like% "C"] <- "compound"
  label_total$obj[label_total$name %like% "G"] <- "glycan"
  label_total$obj[label_total$name %like% "R"] <- "reaction"

  compounds <- label_total[label_total$obj == "compound",]
  compounds$obj <- NULL
  colnames(compounds) <- c("KEGG_ID", "name")
  write.table(compounds, paste(out.file, "_compounds.txt", sep=""), col.names = T, row.names = F)

  reactions <- label_total[label_total$obj == "reaction",]
  reactions$obj <- NULL
  colnames(reactions) <- c("KEGG_ID", "name")
  write.table(reactions, paste(out.file, "_reactions.txt", sep=""), col.names = T, row.names = F)

  modules <- label_total[label_total$obj == "module",]
  modules$obj <- NULL
  colnames(modules) <- c("KEGG_ID", "name")
  write.table(modules, paste(out.file, "_modules.txt", sep=""), col.names = T, row.names = F)

  cat(paste(out.file, "_compounds.txt", sep=""),",", paste(out.file, "_reactions.txt", sep=""), "&",paste(out.file, "_modules.txt", sep=""), "files generated and exported in the working directory! \n")
}
