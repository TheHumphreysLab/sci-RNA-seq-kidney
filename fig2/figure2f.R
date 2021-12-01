#This script describes the workflow of data analysis used in Fig. 2f and following BEAM analysis used in supplemental figures
###by Humphreys Lab


#Fig. 2f (left panel)

library(SeuratData)
library(SeuratDisk)
sciseq <- LoadH5Seurat("read_PT_h5seurat.h5seurat")
library(monocle)
data <- as(as.matrix(sciseq@assays$RNA@data), 'sparseMatrix')

pData <- data.frame(cell_name = colnames(data), PT_0503=sciseq@meta.data$PT_0503,
                    cell_color=sciseq@meta.data$cell_color,sample_color=sciseq@meta.data$sample_color,
                    sample=sciseq@meta.data$sample2_used,row.names = colnames(data))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new("AnnotatedDataFrame", data = pData)
fd <- new("AnnotatedDataFrame", data = fData)

PT_cds <- newCellDataSet(data,phenoData = pd,featureData = fd)
PT_cds <- estimateSizeFactors(PT_cds)
PT_cds <- estimateDispersions(PT_cds)
PT_cds <- detectGenes(PT_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(PT_cds),num_cells_expressed >= 30))

fData(PT_cds)$use_for_ordering <-fData(PT_cds)$num_cells_expressed > 0.05 * ncol(PT_cds)


Vcam1_id <-row.names(subset(fData(PT_cds), gene_short_name == "Kcnip4"))
Fgf1_id <-row.names(subset(fData(PT_cds), gene_short_name == "Fgf1"))
Havcr1_id <-row.names(subset(fData(PT_cds), gene_short_name == "Nrg1"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth,"FR",classify_func = function(x) { x[Vcam1_id,] >= 1 })
cth <- addCellType(cth,"Health",classify_func = function(x) { x[Fgf1_id,] >= 1 })
cth <- addCellType(cth,"Injury",classify_func = function(x) { x[Havcr1_id,] >= 1 })

PT_cds <- classifyCells(PT_cds, cth)

marker_diff <- markerDiffTable(PT_cds[expressed_genes,],cth,cores = 10)

semisup_clustering_genes <-
  row.names(marker_diff)[order(marker_diff$qval)][1:800]

PT_cds <- setOrderingFilter(PT_cds, semisup_clustering_genes)
plot_ordering_genes(PT_cds)
PT_cds <- reduceDimension(PT_cds, max_components = 3,method = 'DDRTree')

Injury_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 2){
    T0_counts <- table(pData(cds)$State, pData(cds)$PT_0503)[,"Type1 Injured PT"]
    return(as.numeric(names(T0_counts)[which (T0_counts == max(T0_counts))]))
  } else {return (1)}}

PT_cds <- orderCells(PT_cds)

plot_cell_trajectory(PT_cds, color_by = "PT_0503",cell_size = 1.5,cell_link_size = 1,
                     show_branch_points=FALSE)+
  theme(legend.position = "right",legend.text=element_text(size=12),legend.title=element_blank(),
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(override.aes = list(size = 7)))+
  scale_color_manual(breaks = c('Healthy PT','Type1 Injured PT','Acute Injury','Repairing','Failed repair'),
                     values=c('#0e8239','#ed5579','#d62728','#ff7f0e','#945015'))




###BEAM analysis
p.heat<-plot_genes_branched_heatmap(PT_cds[row.names(subset(BEAM_res2,
                                                          qval < 1e-30)),],
                                    branch_point = 2,
                                    num_clusters = 2,
                                    cores = 10,
                                    use_gene_short_name = T,
                                    show_rownames = T,
                                    return_heatmap = T,
                                    hmcols = colorRampPalette(c("#FDE725FF","white","#440154FF"))(60))

library(reshape2)
my_row<-p.heat$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)
genes.c1<-my_row[my_row$cluster==1,]$gene
genes.c2<-my_row[my_row$cluster==2,]$gene
heatmap.matrix<-p.heat$heatmap_matrix
hm.failed<-heatmap.matrix[,rev(1:100)]
colnames(hm.failed)<-1:100

hm.repair<-heatmap.matrix[,101:200]
colnames(hm.repair)<-1:100

hm.repair.c1<-hm.repair[genes.c1,]
hm.repair.c1<-as.data.frame(hm.repair.c1)
hm.repair.c1$genes<-rownames(hm.repair.c1)
hm.repair.c1.melted<-melt(hm.repair.c1)
hm.repair.c1.melted$group<-'Successful Repair'

hm.failed.c1<-hm.failed[genes.c1,]
hm.failed.c1<-as.data.frame(hm.failed.c1)
hm.failed.c1$genes<-rownames(hm.failed.c1)
hm.failed.c1.melted<-melt(hm.failed.c1)
hm.failed.c1.melted$group<-'Failed Repair'

hm.c1.melted<-rbind(hm.repair.c1.melted,hm.failed.c1.melted)
hm.c1.melted$genes2<-paste(hm.c1.melted$genes, hm.c1.melted$group, sep = '_')

library(ggplot2)
library(cowplot)
ggplot(hm.c1.melted, aes(x=variable, y=value, color=group))+
  geom_line(aes(group=genes2),alpha=0.08)+
  geom_smooth(aes(group=group),se=F, size=2,method = 'loess')+
  scale_color_manual(values=c('#1267de','#4db030'))+
  theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = 'top',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        legend.title=element_blank(),legend.text=element_text(size=18),axis.title = element_text(size = 18)
  )+
  panel_border(colour='black',size = 0.5)+
  labs(y="Expression", x = "Pseudotime")


hm.repair.c2<-hm.repair[genes.c2,]
hm.repair.c2<-as.data.frame(hm.repair.c2)
hm.repair.c2$genes<-rownames(hm.repair.c2)
hm.repair.c2.melted<-melt(hm.repair.c2)
hm.repair.c2.melted$group<-'Successful Repair'

hm.failed.c2<-hm.failed[genes.c2,]
hm.failed.c2<-as.data.frame(hm.failed.c2)
hm.failed.c2$genes<-rownames(hm.failed.c2)
hm.failed.c2.melted<-melt(hm.failed.c2)
hm.failed.c2.melted$group<-'Failed Repair'

hm.c2.melted<-rbind(hm.repair.c2.melted,hm.failed.c2.melted)
hm.c2.melted$genes2<-paste(hm.c2.melted$genes, hm.c2.melted$group, sep = '_')

ggplot(hm.c2.melted, aes(x=variable, y=value, color=group))+
  geom_line(aes(group=genes2),alpha=0.08)+
  geom_smooth(aes(group=group),se=F, size=2,method = 'loess')+
  scale_color_manual(values=c('#1267de','#4db030'))+
  theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = 'top',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        legend.title=element_blank(),legend.text=element_text(size=18),axis.title = element_text(size = 18)
  )+
  panel_border(colour='black',size = 0.5)+
  labs(y="Expression", x = "Pseudotime")

save.image("~/Analysis2/seurat/monocle2/IRI_output/IRI_beam_e-30.RData")





#plot branch
plot_genes_branched_pseudotime.2 <- function (cds, 
                                              branch_states = NULL, 
                                              branch_point=1,
                                              branch_labels = NULL,
                                              method = "fitting", 
                                              min_expr = NULL, 
                                              cell_size = 0.75,
                                              cell_alpha= 1,
                                              line_size = 1,
                                              nrow = NULL, 
                                              ncol = 1, 
                                              panel_order = NULL, 
                                              color_by = "State",
                                              expression_curve_linetype_by = "Branch", 
                                              trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch", 
                                              reducedModelFormulaStr = NULL, 
                                              label_by_short_name = TRUE,
                                              relative_expr = TRUE,
                                              #gene_pairs = NULL,
                                              ...)
{
  Branch <- NA  
  if (is.null(reducedModelFormulaStr) == FALSE) {
    pval_df <- branchTest(cds, 
                          branch_states=branch_states,
                          branch_point=branch_point,
                          fullModelFormulaStr = trend_formula,
                          reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", 
                          ...)
    fData(cds)[, "pval"] <- pval_df[row.names(cds), 'pval']
  }
  if("Branch" %in% all.vars(terms(as.formula(trend_formula)))) { #only when Branch is in the model formula we will duplicate the "progenitor" cells
    cds_subset <- buildBranchCellDataSet(cds = cds, 
                                         branch_states = branch_states, 
                                         branch_point=branch_point,
                                         branch_labels = branch_labels, 
                                         progenitor_method = 'duplicate',
                                         ...)
  }
  else {
    cds_subset <- cds
    pData(cds_subset)$Branch <- pData(cds_subset)$State
  }
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
  }
  if (integer_expression) {
    CM <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(CM)))
  }
  else {
    cds_exprs <- reshape2::melt(exprs(cds_subset))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)
  
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  if (integer_expression) {
    cds_exprs$adjusted_expression <- round(cds_exprs$expression)
  }
  else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  cds_exprs$Branch <- as.factor(cds_exprs$Branch) 
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
  
  full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula, 
                                            relative_expr = T, new_data = new_data)
  colnames(full_model_expectation) <- colnames(cds_subset)
  
  cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
  if(!is.null(reducedModelFormulaStr)){
    reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
                                                 relative_expr = T, new_data = new_data)
    colnames(reduced_model_expectation) <- colnames(cds_subset)
    cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
  }
  
  # FIXME: If you want to show the bifurcation time for each gene, this function
  # should just compute it. Passing it in as a dataframe is just too complicated
  # and will be hard on the user. 
  # if(!is.null(bifurcation_time)){
  #     cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
  # }
  if (method == "loess")
    cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
  cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
  
  if(!is.null(reducedModelFormulaStr)){
    cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
    cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
  }
  
  cds_exprs$State <- as.factor(cds_exprs$State)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  # if (!is.null(bifurcation_time)) {
  #   q <- q + geom_vline(aes(xintercept = bifurcation_time),
  #                       color = "black", linetype = "longdash")
  # }
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(fill = color_by), size = I(cell_size),shape=21,stroke=0,alpha=cell_alpha)
  }
  if (is.null(reducedModelFormulaStr) == FALSE)
    q <- q + scale_y_log10() + facet_wrap(~feature_label +
                                            pval, nrow = nrow, ncol = ncol, scales = "free_y")
  else q <- q + scale_y_log10() + facet_wrap(~feature_label,
                                             nrow = nrow, ncol = ncol, scales = "free_y")
  if (method == "loess")
    q <- q + stat_smooth(aes(fill = Branch, color = Branch),
                         method = "loess")
  else if (method == "fitting") {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
                                  col = "Branch"),size=I(line_size), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Pre-branch", "AT1", "AT2", "AT1", "AT2")
  }
  
  if(!is.null(reducedModelFormulaStr)) {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
                       color = 'black', linetype = 2, data =  cds_exprs)   
  }
  
  q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
  
  q <- q + theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
  q + expand_limits(y = min_expr)
}

#fr_genes <- row.names(subset(fData(PT_cds),gene_short_name %in% c("Vcam1", "Kcnip4",'Sema5a','Dcdc2a')))
fr_genes <- row.names(subset(fData(PT_cds),gene_short_name %in% c("Vcam1",'Dcdc2a')))
plot_genes_branched_pseudotime.2(PT_cds[fr_genes,], branch_point = 2,color_by = "PT_0503",panel_order=c('Vcam1','Dcdc2a'),
                               cell_size=1.2,ncol = 2)+
  scale_color_manual(breaks = c('Y_35','Y_37'),values=c('tomato','royalblue'))+ 
  theme(legend.direction = "horizontal",legend.position = "none",legend.box = "vertical",axis.title = element_text(size = 15)
        )+
  scale_fill_manual(values=c('#d62728','#945015','#0e8239','#ff7f0e','#ed5579'))+
  labs(y="Expression", x = "Pseudotime")

#r_genes <- row.names(subset(fData(PT_cds),gene_short_name %in% c("Fgf1", "Lrp2",'Ghr','Cubn')))
r_genes <- row.names(subset(fData(PT_cds),gene_short_name %in% c("Fgf1", "Lrp2")))
plot_genes_branched_pseudotime.2(PT_cds[r_genes,],
                                 branch_point = 2,
                                 color_by = "PT_0503",
                                 ncol = 2)

inj_genes <- row.names(subset(fData(PT_cds),gene_short_name %in% c("Havcr1", "Nrg1")))
plot_genes_branched_pseudotime.2(PT_cds[inj_genes,],
                                 branch_point = 2,
                                 color_by = "PT_0503",
                                 ncol = 2)

all_genes <- row.names(subset(fData(PT_cds),gene_short_name %in% c("Vcam1",'Dcdc2a',"Fgf1", "Lrp2","Havcr1", "Nrg1")))
plot_genes_branched_pseudotime.2(PT_cds[all_genes,], branch_point = 2,color_by = "PT_0503",
                                 panel_order=c('Vcam1','Dcdc2a',"Lrp2","Fgf1", "Havcr1", "Nrg1"),
                                 ncol = 2)+
  scale_color_manual(breaks = c('Y_35','Y_37'),values=c('#4db030','#1267de'))+ 
  theme(legend.direction = "horizontal",legend.position = "none",legend.box = "vertical",axis.title = element_text(size = 15)
  )+
  scale_fill_manual(values=c('#d62728','#945015','#0e8239','#ff7f0e','#ed5579'))+
  labs(y="Expression", x = "Pseudotime")



#Fig. 2f (right panel)

sciseq <- LoadH5Seurat("read_PT_h5seurat_2.h5seurat")

data <- as(as.matrix(sciseq@assays$RNA@data), 'sparseMatrix')

pData <- data.frame(cell_name = colnames(data), PT_0503=sciseq@meta.data$PT_0503,
                    cell_color=sciseq@meta.data$cell_color,sample_color=sciseq@meta.data$sample_color,
                    sample=sciseq@meta.data$sample2_used,row.names = colnames(data))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new("AnnotatedDataFrame", data = pData)
fd <- new("AnnotatedDataFrame", data = fData)

PT_cds <- newCellDataSet(data,phenoData = pd,featureData = fd)
PT_cds <- estimateSizeFactors(PT_cds)
PT_cds <- estimateDispersions(PT_cds)
PT_cds <- detectGenes(PT_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(PT_cds),num_cells_expressed >= 30))

fData(PT_cds)$use_for_ordering <-fData(PT_cds)$num_cells_expressed > 0.05 * ncol(PT_cds)


Kcnip4_id <-row.names(subset(fData(PT_cds), gene_short_name == "Kcnip4"))
Fgf1_id <-row.names(subset(fData(PT_cds), gene_short_name == "Fgf1"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth,"FR",classify_func = function(x) { x[Kcnip4_id,] >= 1 })
cth <- addCellType(cth,"Health",classify_func = function(x) { x[Fgf1_id,] >= 1 })


PT_cds <- classifyCells(PT_cds, cth)

marker_diff <- markerDiffTable(PT_cds[expressed_genes,],cth,cores = 10)

semisup_clustering_genes <-
  row.names(marker_diff)[order(marker_diff$qval)][1:200]

PT_cds <- setOrderingFilter(PT_cds, semisup_clustering_genes)
plot_ordering_genes(PT_cds)

PT_cds <- reduceDimension(PT_cds, max_components = 10,method = 'DDRTree')
PT_cds <- orderCells(PT_cds)

#
Health_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 2){
    T0_counts <- table(pData(cds)$State, pData(cds)$PT_0503)[,"Healthy pt"]
    return(as.numeric(names(T0_counts)[which (T0_counts == max(T0_counts))]))
  } else {return (1)}}

plot_cell_trajectory(PT_cds, color_by = "PT_0503",cell_size = 1.5,cell_link_size = 1,
                     show_branch_points=FALSE)+
  theme(legend.position = "right",legend.text=element_text(size=12),legend.title=element_blank(),
        legend.key = element_rect(fill = "white"))+
  guides(color = guide_legend(override.aes = list(size = 7)))+
  scale_color_manual(breaks = c('Healthy pt','Type2 Injured pt','Acute Injury','Repairing','Failed repair'),
                     values=c('#0e8239','#4e97c2','#d62728','#ff7f0e','#945015'))

PT_cds <- orderCells(PT_cds, root_state = Health_state(PT_cds))

plot_cell_trajectory(PT_cds, color_by = "Pseudotime",
                     cell_size = 1.5,cell_link_size = 1.5,show_branch_points=FALSE)+theme_grey()+
  theme(legend.position = "right",legend.text=element_text(size=15),legend.title = element_text( size = 20))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))  

plot_cell_trajectory(PT_cds, color_by = "State",
                     cell_size = 1.5,cell_link_size = 1.5,show_branch_points=FALSE)+theme_grey()+
  theme(legend.position = "right",legend.text=element_text(size=15),legend.title = element_text( size = 20))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22)) 

to_be_tested <- row.names(subset(fData(PT_cds),
                                 gene_short_name %in% c('Vcam1','Dcdc2a',"Lrp2","Fgf1", "Havcr1", "Nrg1")))
cds_subset <- PT_cds[to_be_tested,]
plot_genes_in_pseudotime(cds_subset, color_by = "PT_0503",ncol=2,
                         panel_order=c('Vcam1','Dcdc2a',"Lrp2","Fgf1", "Havcr1", "Nrg1"))+
  scale_color_manual(breaks = c('Healthy pt','Type2 Injured pt','Acute Injury','Repairing','Failed repair'),
                     values=c('#0e8239','#4e97c2','#d62728','#ff7f0e','#945015'))+
  theme(legend.position = "none",axis.title = element_text(size = 15))+
  labs(y="Expression", x = "Pseudotime")
