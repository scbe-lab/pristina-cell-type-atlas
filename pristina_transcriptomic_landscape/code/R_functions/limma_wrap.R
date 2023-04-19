limma_wrap <- function(
  counts, samples.info, condition, MIN.PV = 0.05, MIN.FC = 1.5, adjust = "fdr",
  plot = FALSE){
  # counts = plei_counts[,1:2]
  # samples.info <- plei_samples_info[1:2,]
  # condition = "cellcycle_stage"
  # print(paste0("Condition chosen = ",condition))
  
  samples.info$condition <- samples.info[,colnames(samples.info) == condition ]
  
  print("Creating DGE objects")
  dge <- DGEList(
    counts=counts,
    samples=samples.info,
    group=samples.info$condition
  )
  dge
  keep <- filterByExpr(dge)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  # if(ncol(counts) <= 2) {dge <- t(dge)}
  
  print("Creating model matrix")
  design.condition <- model.matrix( ~ samples.info$condition, dge)
  colnames(design.condition) <- gsub("samples.info[$]condition", "", colnames(design.condition))
  colnames(design.condition) <- gsub("samples.info[$]sample", "", colnames(design.condition))
  names(attr(design.condition,"contrasts")) <- c("condition")
  rownames(design.condition) <- rownames(samples.info)
  
  print(paste0('The largest library has: ',max(dge$samples$lib.size),' counts'))
  print(paste0('The smallest library has: ',min(dge$samples$lib.size),' counts'))
  print(paste0("That's ", round(max(dge$samples$lib.size)/min(dge$samples$lib.size),1), " times bigger"))
  
  
  #ready to run LIMMA
  print("Ready to run limma.")
  
  print("Calculating factors")
  dge <- calcNormFactors(dge, method="TMM")
  
  print("voom normalisation")
  v <- voom(dge, design.condition, plot=plot) #it finds more DE genes using qnorm
  
  print("limma linear model fitting")
  fitv <- lmFit(v, design.condition)
  fitv <- eBayes(fitv)
  sumfit <- summary(decideTests(fitv))
  
  print("creating results")
  top.res.v <- topTable(fitv, coef=ncol(design.condition), 
                        adjust=adjust, sort.by="B",
                        number=Inf)
  
  top.res.v$rank <- ifelse(is.na(top.res.v$P.Value) | is.na(top.res.v$logFC), NA,
                           -log10(top.res.v$P.Value)*top.res.v$logFC *
                             ifelse(top.res.v$P.Value <= MIN.PV & abs(top.res.v$logFC) >= MIN.FC,
                                    ifelse(top.res.v$adj.P.Val <= MIN.PV, 1e4, 1e2), 1));
  
  print("adjusting pvalues")
  top.res.v$DGE.pval <- as.factor(
    ifelse(top.res.v$P.Value <= MIN.PV & abs(top.res.v$logFC) >= MIN.FC,
           ifelse(top.res.v$logFC >= MIN.FC, "up", "down"), "no-sig") );
  
  top.res.v$DGE.padj <- as.factor(
    ifelse(top.res.v$adj.P.Val <= MIN.PV & abs(top.res.v$logFC) >= MIN.FC,
           ifelse(top.res.v$logFC >= MIN.FC, "UP", "DOWN"), "no-sig") );
  
  UPDN <- c(nrow(top.res.v[ top.res.v$DGE.padj == "UP", ]),
            nrow(top.res.v[ top.res.v$DGE.padj == "DOWN", ]))
  updn <- c(nrow(top.res.v[ top.res.v$DGE.pval == "up", ]),
            nrow(top.res.v[ top.res.v$DGE.pval == "down", ]))
  
  
  print("creating plot")
  volcanoplot.limma.v <- EnhancedVolcano(
    data.frame(top.res.v), x = 'logFC', y = 'adj.P.Val',
    lab = rownames(top.res.v),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    ylim = c(0, 6),
    pCutoff = MIN.PV, FCcutoff = MIN.FC, pointSize = 1.0, labSize = 2.0,
    title = paste0("Piwi vs non-Piwi"),
    subtitle = 'limma-voom analysis',
    caption = paste0('log2 FC cutoff: ', MIN.FC, '; p-value cutoff: ',
                     MIN.PV, '\nTotal = ', nrow(top.res.v),
                     ' markers  [ ',UPDN[1],'UP, ',UPDN[2],'DOWN ]'),
    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0,
    drawConnectors = TRUE, widthConnectors = 0.25,
    colConnectors = 'grey30',
    gridlines.major = FALSE,
    gridlines.minor = TRUE,
    col = c("gray","#a2d9a6","darkgray","#4eaf55")
    # max.overlaps = 10,
  ) + coord_flip()
  
  if(plot == TRUE){
    print(volcanoplot.limma.v)
  }
  
  res <- list(
    input_information = samples.info,
    design = design.condition,
    summary_fit = sumfit,
    results = top.res.v,
    volcanoplot = volcanoplot.limma.v
  )
  
  print("DONE.")
  
  
  return(res)
  
}