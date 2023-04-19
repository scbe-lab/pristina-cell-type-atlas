#' from https://stackoverflow.com/questions/44646488/stacked-bar-plot-with-hierarchical-clustering-dendrogram
#', and https://wilkelab.org/cowplot/articles/shared_legends.html
library(tidyverse)
library(ggdendro)
library(vegan)
library(colorspace)
library(cowplot)

barplot_cluster <- function(x){ #add clumethod,color pal, ...
  hc=hclust(dist(x),method="ward.D2")
  hc=reorder(
    hc,
    wts=-as.matrix(x)%*%seq(ncol(x))^2
    ) # vegan::reorder.hclust
  tree=ggdendro::dendro_data(as.dendrogram(hc),type="rectangle")
  
  p1=ggplot(ggdendro::segment(tree))+
    geom_segment(aes(x=y,y=x,xend=yend,yend=xend),lineend="round",size=.4)+
    scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don'x crop half of line between top-level nodes
    scale_y_continuous(limits=.5+c(0,nrow(x)),expand=c(0,0))+
    theme(
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
      axis.title=element_blank(),
      panel.background=element_rect(fill="white"),
      panel.grid=element_blank(),
      plot.margin=margin(5,5,5,0)
    )
  
  x=x[hc$labels[hc$order],]
  x2=data.frame(V1=rownames(x)[row(x)],V2=colnames(x)[col(x)],V3=unname(c(x)))
  lab=round(100*x2$V3)
  lab[lab==0]=""
  
  p2=ggplot(x2,aes(x=factor(V1,level=rownames(x)),y=V3,fill=V2))+
    geom_bar(stat="identity",width=1,position=position_fill(reverse=T))+
    # geom_text(aes(label=lab),position=position_stack(vjust=.5,reverse=T),size=3.5)+
    coord_flip() +
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    scale_fill_manual(values=colorspace::hex(HSV(head(seq(0,360,length.out=ncol(x)+1),-1),.5,1)))+
    theme(
      axis.text=element_text(color="black",size=11),
      axis.text.x=element_blank(),
      axis.ticks=element_blank(),
      axis.title=element_blank(),
      plot.margin=margin(5,0,5,5)
    )
  
  prow <- cowplot::plot_grid(
    p2 + theme(legend.position="none"),
    p1,
    labels = c("A", "B"), rel_widths=c(1,.4))
  # extract the legend from one of the plots
  p2_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    p2 + theme(legend.box.margin = margin(0, 0, 0, 6))
  )
  
  cowplot::plot_grid(prow,p2_legend,rel_widths = c(1, .4))
  
}
