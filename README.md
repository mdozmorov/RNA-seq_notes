# A continually expanding collection of RNA-seq tools

These notes are not intended to be comprehensive. They include notes about methods, packages and tools I would like to explore. For a comprehensive overview of the subject, consider [other bioinformatics resources](https://github.com/mdozmorov/blogs/tree/master/Bioinformatics). Issues with suggestions and pull requests are welcome!

# Table of content

* [Quality control](#quality-control)
* [Clustering](#clustering)
* [Differential expression](#differential-expression)
* [Functional enrichment](#functional-enrichment)
* [Alternative splicing](#alternative-splicing)
  * [circRNAs](#circrnas)
* [Gene fusion](#gene-fusion)
* [Structural variations](#structural-variations)
* [Networks](#networks)
* [Motif enrichment](#motif-enrtichment)
* [Classification](#classification)
* [Visualization](#visualization)
* [Data](#data)
* [Misc](#misc)


## Quality control

- `MultiQC` - Summarization and visualization QC results for multiple samples in one report. Recognizes multiple QC tools. http://multiqc.info/

- `ngsReports` - An R Package for managing FastQC reports and other NGS related log files. https://github.com/UofABioinformaticsHub/ngsReports, http://biorxiv.org/content/early/2018/05/02/313148.abstract

- `sickle` - A windowed adaptive trimming tool for FASTQ files using quality. Post-adapter trimming step. https://github.com/najoshi/sickle

## Clustering

- `clust` - Python package for identification of consistently coexpressed clusters of genes, within and across datasets. Consensus clustering principle. ~50% of genes do not cluster well and thus shouldn't be considered. Compared with seven tools (cross-clustering, k-means, SOMs, MCL, HC, Click, WGCNA) using seven different cluster validation metrics. Outperforms all, produces more focused and significant functional enrichment results. https://github.com/BaselAbujamous/clust
    - Abu-Jamous, Basel, and Steven Kelly. “Clust: Automatic Extraction of Optimal Co-Expressed Gene Clusters from Gene Expression Data.” Genome Biology 19, no. 1 (December 2018). https://doi.org/10.1186/s13059-018-1536-8.


## Differential expression

- `ideal` - Interactive Differential Expression AnaLysis. http://bioconductor.org/packages/release/bioc/html/ideal.html

- Publication-ready volcano plots with enhanced colouring and labeling. https://github.com/kevinblighe/EnhancedVolcano

## Functional enrichment

- `data_analysis_portals.xlsx` - 25 data analysis portals, from the Metascape paper. [Source](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09234-6/MediaObjects/41467_2019_9234_MOESM3_ESM.xlsx)

- `gProfileR` - enrichment of gene lists, from GO to KEGG and others, organism-specific, ortholog conversion

- `Metascape` - multi-gene lists functional analysis, auto gene ID recognition, >40 databases. http://metascape.org/gp/index.html#/main/step1
    - Zhou, Yingyao, Bin Zhou, Lars Pache, Max Chang, Alireza Hadj Khodabakhshi, Olga Tanaseichuk, Christopher Benner, and Sumit K. Chanda. “Metascape Provides a Biologist-Oriented Resource for the Analysis of Systems-Level Datasets.” Nature Communications 10, no. 1 (December 2019): 1523. https://doi.org/10.1038/s41467-019-09234-6.


## Alternative splicing

- `RegTools` - integration of somatic variants from DNA-seq and splice junctions from RNA-seq data to identify variants causing aberrant splicing in cancer. https://regtools.readthedocs.io/en/latest/
    - Feng, Yang-Yang, Avinash Ramu, Kelsy C Cotto, Zachary L Skidmore, Jason Kunisaki, Donald F Conrad, Yiing Lin, et al. “RegTools: Integrated Analysis of Genomic and Transcriptomic Data for Discovery of Splicing Variants in Cancer,” November 25, 2018. https://doi.org/10.1101/436634.

- `vast-tools` - A toolset for profiling alternative splicing events in RNA-Seq data. https://github.com/vastgroup/vast-tools

### circRNAs

- `CIRCpedia` database of cornRNAs from human, mouse, and some model organisms. Ribo-, poly(A)-, RNAse R methods for enriching for circRNAs. http://www.picb.ac.cn/rnomics/circpedia/. `CIRCexplorer2` for the analysis of such experiments, https://circexplorer2.readthedocs.io/en/latest/
    - Zhang et al., “Diverse Alternative Back-Splicing and Alternative Splicing Landscape of Circular RNAs.”

## Gene fusion

- `Arriba` - Fast and accurate gene fusion detection from RNA-Seq data. https://github.com/suhrig/arriba

- `GeneFuse` - Gene fusion detection and visualization. https://github.com/OpenGene/GeneFuse

- `EricScript` is a computational framework for the discovery of gene fusions in paired end RNA-seq data, https://sites.google.com/site/bioericscript/


## Structural variations

- `SQUID` - transcriptomic structural variation caller. Genome segment graph, then rearrange segments so that as many read alignments as possible are concordant with the rearranged sequence. Compared with MUMmer3, DELLY2, LUMPY in simulated settings, and with SOAPfuse, deFuse, FusionCatcher, JAFFA, INTEGRATE tools using real data. https://github.com/Kingsford-Group/squid
    - Ma, Cong, Mingfu Shao, and Carl Kingsford. “SQUID: Transcriptomic Structural Variation Detection from RNA-Seq.” Genome Biology 19, no. 1 (12 2018): 52. https://doi.org/10.1186/s13059-018-1421-5.

- `transindel` - Indel caller for DNA-seq or RNA-seq, https://github.com/cauyrd/transIndel


## Networks

- `GENIE3` - random forest regression detection of gene modules. Input - expression matrix, output - gene x gene square co-regulation matrix. https://github.com/aertslab/GENIE3


## Motif enrichment

- `RcisTarget` - finding enriched motifs in cis-regulatory regions in a gene list. https://github.com/aertslab/RcisTarget


## Classification

- `MLSeq` - Machine learning interface for RNA-Seq data. https://www.bioconductor.org/packages/release/bioc/html/MLSeq.html, Zararsız, Gökmen, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Gozde Erturk Zararsiz, Izzet Parug Duru, and Ahmet Ozturk. “A Comprehensive Simulation Study on Classification of RNA-Seq Data.” Edited by Christian Schönbach. PLOS ONE 12, no. 8 (August 23, 2017): e0182507. https://doi.org/10.1371/journal.pone.0182507.


## Visualization

- `chromoMap` - Interactive Visualization and Mapping of Human Chromosomes, https://cran.r-project.org/web/packages/chromoMap/index.html

- genomation: a toolkit for annotation and visualization of genomic data. https://bioconductor.org/packages/devel/bioc/vignettes/genomation/inst/doc/GenomationManual.html

- `karyoploteR` - An R/Bioconductor package to plot arbitrary data along the genome. https://github.com/bernatgel/karyoploteR

- `pcaExplorer` - Interactive Visualization of RNA-seq Data Using a Principal Components Approach. https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html

- `WIlsON` - Web-based Interactive Omics VisualizatioN, accepts, text files, SummarizedExperiment datasets. 


## Data

- `cBioPortalData` - cBioPortal data as MultiAssayExperiment objects. https://github.com/waldronlab/cBioPortalData

- `curatedTCGAData` - Curated Data From The Cancer Genome Atlas (TCGA) as MultiAssayExperiment objects. https://github.com/waldronlab/curatedTCGAData

- `ARCHS4` - Massive Mining of Publicly Available RNA-seq Data from Human and Mouse. http://amp.pharm.mssm.edu/archs4/download.html

- `PINS` - A novel method for data integration and disease subtyping. http://www.cs.wayne.edu/tinnguyen/PINS/PINS.html

- `recount2` - https://bioconductor.org/help/workflows/recountWorkflow/

- `GREIN` - re-analysis of RNA-seq datasets from GEO. Download processed data, visualization, power analysis, differential expression, functional enrichment analysis, connectivity analysis with LINCS L1000 data. Web-app, http://www.ilincs.org/apps/grein/, source code, https://github.com/uc-bd2k/grein, Docker image, https://hub.docker.com/r/ucbd2k/grein/
    - Al Mahi, Naim, Mehdi Fazel Najafabadi, Marcin Pilarczyk, Michal Kouril, and Mario Medvedovic. “GREIN: An Interactive Web Platform for Re-Analyzing GEO RNA-Seq Data,” October 27, 2018. https://doi.org/10.1101/326223.


## Misc

- `BioJupies` - analysis of GEO/GTEx data or your own gene expression table/FASTQ in autogenerated Jupyter notebook. Rich set of tools for EDA (PCA, Clustergrammer, Library size analysis), Differential expression analysis (Volcano, MA plots), Enrichment analysis (Enrichr, GO, Pathway, TF, Kinase, miRNA enrichments), L1000 signatures. Best suited for two-group analysis. Includes Methods for the selected tools. https://amp.pharm.mssm.edu/biojupies/

- `HGNChelper` - Handy Functions for Working with HGNC Gene Symbols and Affymetrix Probeset Identifiers. https://cran.r-project.org/web/packages/HGNChelper/index.html

- `tximport` - importing transcript abundance datasets from Salmon, Sailfish, kallisto, RSEM, and differential analysis

- `rpkmforgenes` - a Python script for calculating gene expression for RNA-Seq data. http://sandberg.cmb.ki.se/media/data/rnaseq/instructions-rpkmforgenes.html

- Python interface to access reference genome features (such as genes, transcripts, and exons) from Ensembl, https://github.com/openvax/pyensembl
