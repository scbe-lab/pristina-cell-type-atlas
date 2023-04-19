## Gene Age assignment

require(tidyr)

#' Function to manipulate the gfam table, it loads a path to a file and returns
#' the gene--gfam association table, longformat

parse_gfam_table <- function(gfams_table,species_of_interest){
  #1 Load gene families table
  gfams_table <- gfams_table
  colnames(gfams_table)[1] <- "gfam"
  
  #2 subset for your species of interest
  # (one or more than one, parameter provided by user)
  species_filt <- which(colnames(gfams_table) %in% species_of_interest)
  gfams_table_filt <- gfams_table[,c(1,species_filt)]
  
  #3 Format the table of orthogroups to expand species <-> gfam
  gene_gfam <- data.frame(
    id = "0",
    gfam = "0"
  )
  for (species in species_of_interest) {
    species_filt <- which(colnames(gfams_table_filt) == species)
    df_sp <- separate_rows(
      gfams_table_filt[,c(1,species_filt)],
      species,
      sep=", "
    )[,2:1]
    colnames(df_sp) <- c("id","gfam")
    df_sp <- df_sp[df_sp$id != "",]
    gene_gfam <- rbind(gene_gfam,df_sp)
  }
  gene_gfam <- gene_gfam[gene_gfam$id != "0",]
  
  #4 return the result
  return(gene_gfam)

} # one might want to subdivide this table per species.


#' Once we have the table of gene <--> gfam, we do the same with the one
#' of fams and age. Input: dollo parsimony matrix. Output: gfam--age
#' 
#' IMPORTANT: this function only works if collumns corresponding to nodes are
#' properly sorted. It assumes the column with your species of interest is the
#' first one after the one with the gfam names. we might want to revise this!

parse_geneage_table <- function(dollo_table,ages_table){
  #1 Load the gfam age table
  dollo_table <- dollo_table
  colnames(dollo_table)[1] <- "gfam"
  
  ages_table$Name <- paste(
    formatC(
      seq(1:nrow(ages_table)),
      width = nchar(nrow(ages_table)),
      format = "d",
      flag = "0"
    ),
    ages_table$Name, sep = "_"
  )
  
  #2 Parse the table to relevant nodes
  dollo_filt <- which(colnames(dollo_table) %in% ages_table$Node)
  dollo_table_filt <- dollo_table[,c(1,dollo_filt)]
  colnames(dollo_table_filt) <- c("gfam",ages_table$Name)
  dollo_table_filt <- dollo_table_filt[dollo_table_filt[,2] != 0, ]
  dollo_table_filt[,2] <- 1
  
  #3 Create gene age table gfam -- age
  gfams_age <- data.frame(
    gfam = dollo_table_filt$gfam,
    age_num = rowSums(
      dollo_table_filt[,2:ncol(dollo_table_filt)]
    ),
    age = colnames(dollo_table_filt)[
      rowSums(dollo_table_filt[,2:ncol(dollo_table_filt)]) + 1
    ] #nice
  )
  
  #4 return the result
  return(gfams_age)
}


#' Merge both together, generating a table of gene--gfam--age

get_fams_and_age <- function(
  gfams_table,dollo_table,ages_table,species_of_interest
){
  
  if(
    length(species_of_interest) > 1
  ) {
    stop("You have picked more than one species of interest. Pick only one.")
  }
  
  # 1. Gene Families table
  gene_gfam <- parse_gfam_table(
    gfams_table = gfams_table, species_of_interest = species_of_interest
    )
  
  # 2. Gene Age Table
  gene_age <- parse_geneage_table(
    dollo_table = dollo_table, ages_table = ages_table
  )
  
  # 3. Combine the two
  gene_gfam_age <- merge(
    gene_gfam,
    gene_age,
    by.x = 2,
    by.y = 1,
    all = TRUE
  )[,c(2,1,3,4)]
  
  # 4. Return the table gene--gfam--age_num--age
  return(gene_gfam_age)
  
}

