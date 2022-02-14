#This script describes the workflow of data analysis used in Figure 1C
###by Haikuo Li @ Humphreys Lab


library(Matrix)
library(monocle)
library(stringr)
library(plyr)

bulk_count<-as.matrix(bulk_count)
#bulk_count is a bulk integrated count matrix file (sampleXgene). See Methods of our manuscript for more information.

gene <- read.table("gene.txt", quote="\"", comment.char="")
colnames(gene)<-'gene_short_name'
sample <- read.csv("sample.txt", header=FALSE)
#sample is a meta file
#Please define the 3rd column according to your meta file (e.g., healthy or UUO D2)

colnames(sample)<-c('group','index','condition')


pd <- new("AnnotatedDataFrame", data = sample)
fd <- new("AnnotatedDataFrame", data = gene)
cds_kidney <- newCellDataSet(as(bulk_count, "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.1, expressionFamily=negbinomial.size())
cds_kidney <- estimateSizeFactors(cds_kidney)
cds_kidney <- estimateDispersions(cds_kidney)

cds_kidney <- detectGenes(cds_kidney, min_expr = 1)
#print(head(fData(cds_kidney)))
expressed_genes <- row.names(subset(fData(cds_kidney),
                                    num_cells_expressed >= 2))

fData(cds_kidney)$use_for_ordering <-
  fData(cds_kidney)$num_cells_expressed > 0.1 * ncol(cds_kidney)
#print(head(fData(cds_kidney)))

plot_pc_variance_explained(cds_kidney, return_all = F)

#print(head(pData(cds_kidney)))
clustering_DEG_genes <-
  differentialGeneTest(cds_kidney[expressed_genes,],
                       fullModelFormulaStr = '~condition',
                       cores = 8) 
##


cds_kidney_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:5000]
cds_kidney <-setOrderingFilter(cds_kidney, ordering_genes = cds_kidney_ordering_genes)

cds_kidney <-reduceDimension(cds_kidney,max_components =15, method = 'DDRTree')
cds_kidney <- orderCells(cds_kidney)

Health_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 2){
    T0_counts <- table(pData(cds)$State, pData(cds)$condition)[,"Health"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

cds_kidney <-orderCells(cds_kidney, root_state = Health_state(cds_kidney))

plot_cell_trajectory(cds_kidney, color_by = "condition",cell_size = 7,cell_link_size = 1,
                     show_branch_points=FALSE) + theme_grey()+
  theme(legend.position = "right",legend.text=element_text(size=15),legend.title=element_blank(),
        legend.key = element_rect(fill = "white"))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))+
  scale_color_manual(breaks = c('Health','IRI 6hrs','IRI Day2','IRI Day7','IRI Day14','IRI Day28',
                               'UUO Day2','UUO Day4','UUO Day6','UUO Day10','UUO Day14'),
                     values=c('#98df8a','#eb3636','#eb6c36','#eb9d36','#ebcd36','#caeb36',
                              '#36d6eb','#36afeb','#369deb','#3675eb','#3636eb'))

plot_cell_trajectory(cds_kidney, color_by = "State",
                     cell_size = 7,cell_link_size = 1.5,show_branch_points=FALSE)+theme_grey()+
  theme(legend.position = "right",legend.text=element_text(size=15),legend.title = element_text( size = 20))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))                                 

plot_cell_trajectory(cds_kidney, color_by = "Pseudotime",
                     cell_size = 7,cell_link_size = 1.5,show_branch_points=FALSE)+theme_grey()+
  theme(legend.position = "right",legend.text=element_text(size=15),legend.title = element_text( size = 20))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))                                 

