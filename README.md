# CM_16S_cross-sectionalF245
Data processing and analysis of cross-sectional microbiome

Monzó et al., Dietary restriction mitigates the age-associated decline in mouse B cell receptor repertoire diversity, Cell Reports (2023), https://doi.org/10.1016/j.celrep.2023.112722

### Running cutadapt file:
```
slurm cmd_fastq_trimm.sh
```
### Running general dada2 pipeline  
It calls to R script "dada2_F245.Rmd".  
This script runs the dada2 pipeline until taxa annotation, generating intermediate files for later analysis, QC and processing comparisons.
```
Rscript dada2_F245.Rmd
```
### Running general QC
```
jupyter notebook depth_plotting_F2.ipynb
``` 
### Generating alpha and beta diversity tables for analysis  
Script "dada2_F1_abphylo.R" calculates basic tables for analysis of alpha and beta diversity
```
Rscript dada2_F2_abphylo.R --path ~/workspace/16S_final/CM_16S_cross-sectionalF245/ --nochim ../analysis/seqtab_merge3/CLEAN_merged_seqtabNochim_20210215.rds --taxa ../analysis/seqtab_merge3/CLEAN_taxonomy_merged_20210215.rds --metadata ../metadata/metadata_ready2.csv --count_tab ../analysis/seqtab_merge3/mergedQC/CLEAN_ASVs_counts_merged_20210215.tsv
``` 
### Plotting and general analysis of alpha and beta diversity  
Alpha and beta diversity plotting and basic per-timepoint stats are run as ipythons, as it makes it easier to keep track of the plots as they are being generated.  
```
jupyter notebook alpha_diversity_F2.ipynb  
jupyter notebook beta_diversity_F2.ipynb
```
### Running DESeq2 in each timepoint for comparisons
```
Rscript deseq_timepoints.Rmd
```
### Comparing numbers of differentially abundant ASVs
```
jupyter notebook Num_DiffAb_ASVs.ipynb
```
