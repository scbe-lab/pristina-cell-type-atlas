id_module_kME <- function(modulecolors, datkme){
  genes <- names(modulecolors)
  ngenes <- length(genes)
  
  res <- data.frame(
    id = genes,
    module = modulecolors,
    kME = 0
  )
  
  colnames(datkme) <- sub("MM_", "", colnames(datkme))
  
  for (i in 1:ngenes) {
    genename <- genes[i]
    modulename <- modulecolors[i]
    kME <- datKME[
      rownames(datkme) == genename,
      colnames(datkme) == modulename
    ]
    res[i,3] <- kME
  }
  res$kME <-  as.numeric(res$kME)
  
  return(res)
}

ParseNetwork <- function(
  graph,
  list_attr = list(
    tflist, efflist, classANDcolor, td_hk, gfams, funcat, markerlist, ...
  )
) {
  # 1.1 asigna lo q tenga que asignar de datos para igraph (TF class, color..)
  df_attr <- data.frame(
    id = V(graph)$name,
    index = seq(
      1:length(V(graph))
    )
  )
  
  # merge data
  for (i in list_attr){
    df_attr <- merge(
      df_attr,
      i,
      by.x = 1,
      by.y = 1,
      all.x = T
    ) 
  }
  
  # cleanup data
  df_attr[is.na(df_attr)] <- ""
  df_attr <- df_attr[order(df_attr$index),] #IMPORTANT
  df_attr <- df_attr[!duplicated(df_attr$id),]
  rownames(df_attr) <- NULL
  
  # slap together with the network
  for (i in colnames(df_attr)[-1]){
    graph <- set_vertex_attr(graph, paste0(i), value = df_attr[[i]])
  }
  #' output of the function: a graph with additional 
  #' attributes to the different nodes and edges
  output <- list(
    graph,
    df_attr
  )
  return(output)
}

divide_into_components <- function(x,CCs,impl="create_from_scratch") {
  require(igraph)
  lis <- list()
  for( i in unique(CCs$member)) {
    print(i)
    vids_i <- which(
      V(x)$name %in% CCs$id[CCs$member == i]
    )
    nw_i <- induced_subgraph(
      graph = x,
      vids = vids_i,
      impl = impl
    )
    V(nw_i)$color <- V(nw_i)$genecolor
    lis[[i]] <- nw_i
  }
  return(lis)
}

centrality_by_network <- function(list_grns){
  l <- list()
  for (i in 1:length(list_grns)){
    g <- list_grns[[i]]
    v <- V(g)$name[V(g)$TFclass != ""] #add ifelses or cases here depending on user
    l[[i]] <- closeness(
      g,
      vids = v,
      mode = c("all"),
      weights = NULL,
      normalized = F
    )
    names(l)[i] <- names(rev(sort(table(V(g)$module))))[1]
  }
  return(l)
}

#' cross_connections, id_modules

#' There is a myriad ways to retrieve this information.
#' 
#' datKME can also be used for this purpose. datKME stores information on 
#' how well every gene connects with the eigengenes. The only problem of the 
#' datkME approach is that you are not technically counting how many genes 
#' your gene is connecting to. Thus it cannot retrieve probably the occassional
#' co-expression of two genes in the same cell type.
#' 
#' One could also try to parse the TOM matrix one by one and counting all
#' the prominent connection values.
#' 
#' For the time being i like this approach.

#' Function one:
#' generate a table with genes in rows and modules in columns.
#' Values correspond to number of genes on each module.
#' This can only be done if modules are factors. TAKE INTO ACCOUNT

connections_to_module_per_gene <- function(i, network, id_module){
  colnames(id_module) <- c("id","module")
  v <- names(
    unlist(
      ego(network, order=1, nodes = i, mode = "all", mindist = 0)
    )
  )
  e <- id_module$module[
    id_module$id %in% v
  ]
  r <- as.vector(table(e))
  return(r)
}

#' Function two:
#' Figures WHAT are the genes connecting.

cross_connected_genes <- function(cross_connections, id_modules){
  modules <- as.character(unique(id_modules[,2]))
  colnames(cross_connections) <- as.character(colnames(cross_connections))
  res <- list()
  t = 1
  for (i in modules) {
    
    col_module_i <- which(colnames(cross_connections) == i ) #grab module of the gene we measuring
    
    for (j in modules) {
      
      if (i == j ) next # skip if checking genes connecting to genes of its own module
      col_module_j <- which(colnames(cross_connections) == j ) # grab second module
      genes_i_j <- cross_connections$id[
        cross_connections[,col_module_i] != 0 & 
          cross_connections[,col_module_j] != 0
      ]
      if (length(genes_i_j) != 0) {
        if (
          paste0(j,"_",i) %in% names(res)
        ) { next } # skip if module i and j have been already compared
        else {
          res[[t]] <- genes_i_j
          names(res)[t] <- paste0(i,"_",j)
          t = t+1
        }
      }
    }
  }
  
  return(res)
  
  }
