# Comprehensive single-cell transcriptional profiling defines epithelial injury responses during kidney fibrogenesis
This repository documents the scripts to generate data for our manuscript studying mouse kidney fibrogenesis with sci-RNA-seq3. <link> <br>

The work is published in this [manuscript](https://doi.org/10.1016/j.cmet.2022.09.026) at Cell Metabolism. For citation:
```
Li, H., Dixon, E. E., Wu, H., & Humphreys, B. D. (2022). Comprehensive single-cell transcriptional profiling defines shared and unique epithelial injury responses during kidney fibrosis. Cell metabolism, 34(12), 1977-1998.
DOI: https://doi.org/10.1016/j.cmet.2022.09.026 (PMID: 36265491)
```

Our detailed step-by-step protocol of sci-RNA-seq library generation can be found [here](https://doi.org/10.1016/j.xpro.2022.101904) at STAR Protocols.
```
Li, H., & Humphreys, B. D. (2022). Mouse kidney nuclear isolation and library preparation for single-cell combinatorial indexing RNA sequencing. STAR protocols, 3(4), 101904.
DOI: https://doi.org/10.1016/j.xpro.2022.101904 (PMID: 36595916)
```
<br>

Project workflow:<br>
<img src="https://haikuoli.github.io/files/sciseq-scheme.png" alt="Project workflow"><br><br>

The raw (.fastq) and pre-processed (count matrix) data and metadata have been deposited in NCBI’s Gene Expression Omnibus and are available through GEO Series accession number GSE190887.  <br> <br>
A searchable database, including gene expression in all kidney cell types and PT cell types is available at our Kidney Interactive Transcriptomics (K.I.T.) website: http://humphreyslab.com/SingleCell/.<link> <br><br>
Pre-processing of raw fastq files was performed as previously described (Cao et al. Nature 2019; Cao et al. Science 2020): https://github.com/JunyueC/sci-RNA-seq3_pipeline<br>

### Descriptions

#### 1. Scripts for Figure 1<br>
Quality control and single-cell clustering<br>
Pseudobulk trajectory ordering (Fig. 1c)<br>
Single-cell visualization (Fig. 1d)<br>
Gene expression visualization (Fig. 1e)<br>


#### 2. Scripts for Figure 2<br>
Single-cell visualization (Fig. 2a)<br>
Gene expression visualization (Fig. 2b)<br>
Single-cell TF activity analysis (Fig. 2c)<br>
Single-cell pathway activity analysis (Fig. 2d)<br>
Single-cell trajectory ordering (Fig. 2f)<br>
Fate prediction analysis with time-series scRNAseq data (Fig. 2g)<br>

#### 3. Scripts for Figure 6<br> 
Quality control and single-cell subclustering<br>
Single-cell visualization (Fig. 6a/b)<br>
Gene expression visualization (Fig. 6c/d/f)<br>

#### 4. Scripts for Figure 7<br>
Quality control and single-cell subclustering<br>
Single-cell visualization (Fig. 7a/b)<br>
Gene module scoring (Fig. 7e)<br>
Commands used for CellPhoneDB analysis<br>
Statistical analysis (Fig. 7f/g)<br>

#### 5. Meta data<br>
RT barcode - Sample Lookup Table<br>
Meta info of all cells used for clustering<br>

<br>
**************<br>


Find us on Twitter: 
<br/>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @HumphreysLab</a>
<br/><br/>
