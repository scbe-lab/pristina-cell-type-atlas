#' Params
#' @param genelist must be a object of type list() containing an undetermined number of character vectors. Each component of the list is a char vector of gene names. This function will perform gene ontology enrichment using the char vector as test sample and the totality of your organism's genes as population (Universe). 
#' @gene_universe is a vector of characters containing the genes you want to use as universe. This is normally the totality of genes in your species, or the totality of genes with recollected expression in your dataset.
#' @alg is the algorithm of choice, default is 'Classic', but 'elim' is advised
#' @cols is a char vector of the colors used in the barplots
#' 
require(topGO)
# getGOs(list(tfs,egs,active_tfs,target_genes,V(graph)$name),toy_universe)
getGOs <- function(
  genelist,
  gene_universe,
  gene2GO,
  alg = "Classic",
  stat = "fisher",
  cols=rainbow(length(genelist))
) {
  res <- list(
    GOtable = list(),
    GOplot = list()
  )
  geneID2GO <- gene2GO
  for (j in 1:length(genelist)){
    print(paste0("Starting analysis ",j," of ",length(genelist)))
    x <- genelist[[j]]
    numgenes <- length(x)
    allgenes <- data.frame(
      id = gene_universe)
    allgenes$diffreg <- 0
    allgenes$diffreg[allgenes$id %in% x] <- 1
    genelist_x <- allgenes$diffreg
    names(genelist_x) <- allgenes$id
    signif <- function (allScore) {return(allScore == 1)} # function for TopGO which retrieves the significant genes as those of interest
    # this function should be replaced for something much more efficient. Check other versions of these scripts run by other people
    
    #load data
    #' This is by far the longest and most computationally expensive step.
    #' Can this step be outside the loop and just keep inside the loop 
    #' the bit of which ones to keep?
    GOdata_x <- new( 
      "topGOdata",
      ontology = "BP",
      allGenes = genelist_x,
      geneSelectionFun = signif,
      annot = annFUN.gene2GO,
      gene2GO = geneID2GO
    )
    
    #test
    res_x <- runTest(GOdata_x, algorithm = alg, statistic = "fisher")
    
    #table
    allres_x <- GenTable(
      GOdata_x,
      classicFisher = res_x,
      orderBy = "classicFisher", #see if this can be replaced for a static thing that does not change with the test's name
      ranksOf = "classicFisher",
      topNodes = 30,
      numChar = 5000
    )
    allres_x <- allres_x[allres_x$Annotated > 5,]
    allres_x$classicFisher[allres_x$classicFisher=="< 1e-30"] <- 1e-30
    
    maxgenesplot <- ifelse(nrow(allres_x) > 10, 10, nrow(allres_x))
    res$GOtable[[j]] <- allres_x
    
    #plot
    res$GOplot[[j]] <- list(
      height = as.vector(
        rev(
          -log(
            as.numeric(allres_x$classicFisher[1:maxgenesplot]) * 
              length(x)
          )
        )
      ),
      horiz=T,
      border=F,
      names.arg=rev(allres_x$Term[1:maxgenesplot]),
      sub="-log10 p-value",
      main="GOs (elim,\n> 5 annot genes)",
      las=1,
      col=cols[j],
      space=c(0.5),
      cex.names=.7,
      plot=FALSE
    )
  }
  names(res$GOtable) <- names(genelist)
  return(res)
}
