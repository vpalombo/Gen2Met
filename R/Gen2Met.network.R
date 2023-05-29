#Script: Gen2Met.network.R
#License: GPLv3 or later
#Modification date: 2023-10-15
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Run Gen2Met core analysis.

Gen2Met.network <- function(
  in.file=NULL,
  species=NULL,
  path.file=NULL,
  out.file="Gen2Met_network",
  physics = F
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

  path.file.name <- paste(path.file)
  if (file.exists(path.file.name)){
    path <<- read.table(path.file.name, header=T, stringsAsFactors = F)
    #cat("Path file imported! \n")
    cat("\n")
    #check the path file arguments
    Gen2Met.checkPathfile(path.file = path)
    }

  cat("\nPrerequisite check passed!\n\nGen2Met is running ... \n")
  cat(" Please wait... It could take a while depending on the number of genes investigated! \n")
  genes <- paste(species, ":", genelist$entrezgene_id, sep="")

  path2genes <- NULL
  mod2genes <- NULL
  gene2Reaction <- NULL
  for (i in 1:length(genes)){
    out_nodes <- NULL
    print(i)
    print(Sys.time())
    ort <- names(KEGGREST::keggGet(genes[i])[[1]][["ORTHOLOGY"]])
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


  #write.table(path2genes, paste(out_file, "_path2genes.txt", sep=""), col.names = T, row.names = F, sep="$")
  #write.table(mod2genes, paste(out_file, "_mod2genes.txt", sep=""), col.names = T, row.names = F, sep="$")
  #write.table(gene2Reaction, paste(out_file, "_gene2Reaction.txt", sep=""), col.names = T, row.names = F, sep="$")

  path2mod <- NULL
  comp2mod <- NULL
  modules <- unique(mod2genes$id_to)

  cat("\nA total of", paste(length(modules), "modules were associated with your genelist based on information retrieved from KEGG!", "\n\nGen2Met is searching for compounds in association with these modules ... \n"))
  cat(" Please wait... It could take a while depending on the number of modules investigated! \n")

  #modules <- str_replace_all(modules, pattern = "vvi_", "")
  #query on KEGG
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

  path2mod <- unique(path2mod)
  comp2mod <- unique(comp2mod)

  #eliminate pathway if path.file was specified
  if (file.exists(path.file.name)){
    path2mod <- path2mod[path2mod$id_to%in% stringr::str_replace_all(path$pathway_ID, pattern=paste("path:",species,sep = ""), "map"),]
  }

  #write.table(path2mod, paste(out_file, "_path2mod.txt", sep=""), col.names = T, row.names = F, sep="$")
  #write.table(comp2mod, paste(out_file, "_comp2mod.txt", sep=""), col.names = T, row.names = F, sep="$")

  #lavoriamo sulle reazioni
  reaction <- unique(gene2Reaction$id_from)

  cat("\nA total of", paste(length(reaction), "reactions were associated with your genelist based on information retrieved from KEGG!", "\n\nGen2Met is searching for compounds in association with these reactions ... \n"))
  cat(" Please wait... It could take a while depending on the number of reactions investigated! \n")

  comp2Reaction <- NULL
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
  comp2Reaction <- unique(comp2Reaction)
  #reorder columns based on common order used before
  comp2Reaction <- comp2Reaction[,colnames(path2genes)]
  #write.table(comp2Reaction, paste(out_file, "_comp2Reaction.txt", sep=""), sep="$", col.names = T, row.names = F)

  pathways <- unique(c(path2mod$id_to, path2genes$id_to))

  cat("\nGenes in your genelist are involved in a total of", paste(length(pathways),"pathways \n"))

  if (file.exists(path.file.name)){
    cat("\nBut you specified", paste(length(path$pathway_ID), "pathways !", "\n"))
    pathways <- pathways[pathways %in% stringr::str_replace_all(path$pathway_ID, pattern=paste("path:",species,sep = ""), "map")]
    cat("\n", paste(length(pathways)), "pathways are investigable!", "\n\nGen2Met is searching for interactions among these pathways ... \n")
  }else{
    cat("and you did not specified any pathway list of interest. All KEGG pathways will be considered in the analysis! \n")
  }
  cat(" Please wait... It could take a while depending on the number of pathways investigated! \n")

  #query on KEGG
  path2path <- NULL
  for (i in 1:length(pathways)){
    out_nodes <- NULL
    print(i)
    print(Sys.time())
    name <- KEGGREST::keggGet(pathways[i])[[1]]$PATHWAY_MAP
    dataDownloaded <- KEGGREST::keggGet(pathways[i])[[1]]
    out_nodes <- data.frame(to = KEGGREST::keggGet(pathways[i])[[1]]$REL_PATHWAY)
    if (dim(out_nodes)[1] != 0) {
      out_nodes$id_to <- rownames(out_nodes)
      out_nodes <- out_nodes[(out_nodes$id_to %in% pathways),]
      if(nrow(out_nodes)>0){
        rownames(out_nodes) <- NULL
        out_nodes$id_from <- pathways[i]
        out_nodes$from <- name
        path2path <- rbind(path2path, out_nodes)
      }
    }
  }
  path2path <- unique(path2path)


  #eliminate pathway if path.file was specified
  if (file.exists(path.file.name)){
    path2path <- path2path[path2path$id_to %in% stringr::str_replace_all(path$pathway_ID, pattern=paste("path:",species,sep = ""), "map"),]
    path2path <- path2path[path2path$id_from %in% stringr::str_replace_all(path$pathway_ID, pattern=paste("path:",species,sep = ""), "map"),]
  }

  #eliminate pathway if path.file was specified
  if (file.exists(path.file.name)){
    path2genes <- path2genes[path2genes$id_to %in% stringr::str_replace_all(path$pathway_ID, pattern=paste("path:",species,sep = ""), "map"),]
  }

  #eliminate pathway if path.file was specified
  if (file.exists(path.file.name)){
    path2mod <- path2mod[path2mod$id_to %in% stringr::str_replace_all(path$pathway_ID, pattern=paste("path:",species,sep = ""), "map"),]
  }


  cat(" Global Pathways will be excluded in the visual output! \n")
  path_to_remove <- c("map01100", "map01110", "map01120") #because too general pathways

  if (length(path2genes$to)>0){
    path2genes <- path2genes[!(path2genes$id_to %in% path_to_remove),]
    path2genes$interaction <- "gene2path"
  }

  if (length(path2genes$to)>0){
    path2genes <- path2genes[!(path2genes$id_to %in% path_to_remove),]
  }

  if (length(path2path$to)>0){
    path2path <- path2path[!(path2path$id_to %in% path_to_remove),]
    path2path <- path2path[!(path2path$id_from %in% path_to_remove),]
    path2path$interaction <- "path2path"
  }

  if (length(path2mod$to)>0){
    path2mod <-  path2mod[!(path2mod$id_to %in% path_to_remove),]
    path2mod$interaction <- "mod2path"
  }

  if (length(mod2genes$to)>0){
    mod2genes$interaction <- "gene2mod"
  }


  if (length(comp2mod$to)>0){
    comp2mod$interaction <- "mod2comp"
  }

  if (length(gene2Reaction$to)>0){
    gene2Reaction$interaction <- "reaction2gene"
    comp2Reaction$interaction <- "comp2reaction"
  }

  data <- data.frame(rbind(path2genes,
                           path2path,
                           path2mod,
                           mod2genes,
                           comp2mod,
                           gene2Reaction,
                           comp2Reaction
  ), stringsAsFactors = T)

  edges <- data.frame(source=data$id_from, target=data$id_to, interaction=data$interaction)
  #save the edges
  write.table(edges, paste(out.file, "_edges.txt", sep=""),  sep = "$", col.names = T, row.names = F)

  label <- data.frame(name=data$id_to, shared=data$to)
  label <- unique(label)
  label2 <- data.frame(name=data$id_from, shared=data$from)
  label2 <- unique(label2)
  label_total <- rbind(label, label2)
  label_total <- unique(label_total)

  #eliminate duplicate with similar name
  label_total2 <- NULL
  names <- unique(label_total$name)
  for (z in 1:length(names)){
    a <- label_total[label_total$name == names[z],]
    if (length(a$name)>1){
      a <- data.frame(a[1,])
    }
    label_total2 <- rbind(label_total2, a)
  }

  label_total <- unique(label_total2)
  label_total <- label_total[order(label_total$shared, decreasing = F),]

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
  #write.table(compounds, paste(out.file, "_compounds.txt", sep=""), col.names = T, row.names = F)
  write.table(label_total, paste(out.file, "_nodes.txt", sep=""), sep = "$", col.names = T, row.names = F)

  #create the attributes of diagram objects
  label_total$shape <- 1

  label_total$shape[label_total$obj == "path"] <- "square"
  label_total$shape[label_total$obj == "glycan"] <- "star"
  label_total$shape[label_total$obj == "gene"] <- "diamond"
  label_total$shape[label_total$obj == "compound"] <- "hexagon"
  label_total$shape[label_total$obj == "reaction"] <- "triangle"
  label_total$shape[label_total$obj == "module"] <- "triangleDown"

  label_total$type <- label_total$obj
  label_total$color<- "black"

  label_total$num <- c(1:length(label_total$name))

  dm <- edges[,1:2]
  colnames(dm) <- c("from", "to")

  #create diagram nodes
  ndf <- data.frame(
    id = label_total$num,
    label = paste(label_total$shared),
    shape = label_total$shape,
    group = label_total$type)
  #extract the interactions uncovered
  for (i in 1:length(label_total$num)){
    dm$from[paste(dm$from) == paste(label_total$name[i])] <- paste(label_total$num[i])
    dm$to[paste(dm$to) == paste(label_total$name[i])] <- paste(label_total$num[i])
  }
  #create the diagram edges
  edf <- data.frame(
    from = dm$from,
    to = dm$to)
  #legend color
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #create the diagram and export a html file
  graph <- visNetwork::visNetwork(ndf, edf, height = "1920px", width = "1080px", main = paste("Gen2Met visualization", out.file, sep=" "))
  if (physics == F){
    graph <- visNetwork::visIgraphLayout(graph)
  }
  graph <- visNetwork::visOptions(graph, highlightNearest = TRUE, nodesIdSelection = TRUE)
  visNetwork::visSave(graph, file = paste(out.file,".html", sep=""))
  cat("Well done! Diagram visualization was created and exported. \n")
}
