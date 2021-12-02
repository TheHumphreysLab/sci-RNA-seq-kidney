##Commands used in CellphoneDB analysis

cellphonedb method statistical_analysis sciseq_meta.txt sciseq.h5ad --counts-data=gene_name --output-path sciseq_cpdb --iterations=2000 --result-precision=6 --threads=20

cellphonedb plot heatmap_plot --pvalues-path ./pvalues.txt --output-path heatmap_plot --count-name heatmap_count.png --log-name heatmap_log_count.png ../sciseq_meta.txt

#then we ran similar procedures on demultiplexed datasets to obtain CCI of each disease condition
#See Methods for more information