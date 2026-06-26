# Temporal Uncoupling of Radial Glia Lineage Progression in Cortical Organoid

Raw figures/tables as well as scripts for Stouffer, Miranda, Pauler, Pipicelli et al. Nature 2026

---

## Quick guide

### scRNA-Seq analysis followed this order of analysis
- cellranger
- initial_analysis_organoids
- analysis_embryo_comp

### MADM-CloneSeq
- MADM-CloneSeq data

### statistics and figures for clonal data
- clonal_analysis

### statistics for histology data
- statistics

---

## More detailed description of folder content

### cellranger
contains config files for cellranger analysis as well as meta data linking barcodes used for demultiplexing

### initial_analysis_organoids
initial processing of organoid data

### analysis_embryo_comp
main part of the analysis and all figures

Quick guide to script - figure link (all in analysis_embryo_comp):

- Fig. 1g, h, j: step07.R
- Fig. 1i: step09.R
- Ext. Data Fig 2, 3: step05.R, step06
- Fig. 4b-h: step14.R
- Fig. 4i: step15.R
- Fig. 4j-l: step17.R
- Ext.Data Fig 4a-d: step08.R
- Ext. Data Fig 4e: step11
- Ext. Data Fig 4f: step10
- Ext. Data Fig 5, 6: step04.R
- Ext. Data Fig 7a,b: step07.R
- Ext. Data Fig 7c: step20.R
- Ext. Data Fig 16abc: clonal_analysis
- Ext. Data Fig 16def: clonal_analysis
- Ext. Data Fig 21: clonal_analysis
- Ext. Data Fig 22a: step16.R
- Ext. Data Fig 22b,d:step14.R
- Ext. Data Fig 23a,b,c,d: step14.R
- Ext. Data Fig 23e: step21
- Ext. Data Fig 24b,c,d,e,f,g: step19.R
- Ext. Data Fig 25: step17.R

### MADM-CloneSeq: 
contains info about MADM-CloneSeq used or created in above analysis as well as some stand alone analyses

Note that PDFs were compressed to save space

