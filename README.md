# A continually expanding collection of RNA-seq tools

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

RNA-seq related tools and genomics data analysis resources. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Pipelines](#pipelines)
  - [Preprocessing](#preprocessing)
  - [Aligners](#aligners)
    - [Long-read](#long-read)
  - [Analysis](#analysis)
- [Quality control](#quality-control)
- [Imputation](#imputation)
- [Deconvolution](#deconvolution)
- [Batch effect](#batch-effect)
- [Clustering](#clustering)
- [Timecourse](#timecourse)
- [Allele-specific expression](#allele-specific-expression)
- [Differential expression](#differential-expression)
- [Pathways, Functional enrichment](#pathways-functional-enrichment)
  - [Transcription regulators](#transcription-regulators)
- [Non-canonical RNAs](#non-canonical-rnas)
  - [Alternative splicing](#alternative-splicing)
  - [miRNAs](#mirnas)
  - [lncRNAs](#lncrnas)
  - [circRNAs](#circrnas)
  - [Gene fusion](#gene-fusion)
  - [Isoforms](#isoforms)
- [CNVs and Structural variations](#cnvs-and-structural-variations)
- [Networks](#networks)
  - [Transcription regulators](#transcription-regulators)
- [Integrative](#integrative)
- [Classification](#classification)
- [Visualization](#visualization)
- [Data](#data)
  - [Genes](#genes)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Pipelines

- [The ENCODE pipelines](https://www.encodeproject.org/pages/pipelines/) (RNA-seq, microRNA-seq, long-read RNA, ChIP-seq, ATAC-seq, DNAse-seq, methylation, WGBS, 3D), including quality control steps. Figures - pipeline schemas. Implemented in [WDL](https://openwdl.org/), available on DockerHub. Supplementary tools - [CAPER](https://github.com/ENCODE-DCC/caper) - Cromwell/WDL wrapper for Python, provides backends for HPC job submission systems, [CROO](https://github.com/ENCODE-DCC/croo) - Cromwell output organizer, [accession](https://github.com/ENCODE-DCC/accession). <details>
    <summary>Paper</summary>
    Hitz, Benjamin C., Jin-Wook Lee, Otto Jolanki, Meenakshi S. Kagda, Keenan Graham, Paul Sud, Idan Gabdank, et al. “The ENCODE Uniform Analysis Pipelines.” Preprint. Bioinformatics, April 6, 2023. https://doi.org/10.1101/2023.04.04.535623.
</details>

- [RNA-seq pipeline developed by NASA GeneLab](https://github.com/nasa/GeneLab_Data_Processing). FastQC/MultiQC, TrimGalore, STAR (two-pass mode), RSEM (to quantify isiforms), DESeq2. Can normalize using ERCCs. produces unnormalized and normalized counts, more. PCA, heatmaps, other visualization. [Step-by-step instructions and commands](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/GL-DPPD-7101-C.md)
    - Galazka, Jonathan. “[NASA GeneLab RNA-Seq Consensus Pipeline: Standardized Processing of Short-Read RNA-Seq Data](https://doi.org/10.1101/2020.11.06.371724),” bioRxiv, November 10, 2020

- [GEO2RNAseq](https://anaconda.org/xentrics/r-geo2rnaseq) - an R-based RNA-seq processing pipeline, from FASTQ files or directly from GEO. Reads meta-data, oriented for two-group differential analysis. Competitors: Galaxy, the Total RNA-seq Analysis Package for R (TRAP), EasyRNASeq, READemption. Installable through Conda https://anaconda.org/xentrics/r-geo2rnaseq
    - Seelbinder, Bastian, Thomas Wolf, Steffen Priebe, Sylvie McNamara, Silvia Gerber, Reinhard Guthke, and Joerg Linde. “[GEO2RNAseq: An Easy-to-Use R Pipeline for Complete Pre-Processing of RNA-Seq Data](https://doi.org/10.1101/771063).” Preprint. Bioinformatics, September 16, 2019

- [RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow) - A repository for setting up a RNAseq workflow. Detailed instructions and code for each analysis and visualization step.

### Preprocessing

- [Check strandedness of RNA-Seq fastq files](https://github.com/betsig/how_are_we_stranded_here)

- [Illumina Instrument Type from fastq](https://www.biostars.org/p/198143/)

- [adapters](https://github.com/stephenturner/adapters) - Adapter sequences for trimming, by Stephen Turner

- [bamcov](https://github.com/fbreitwieser/bamcov) - Quickly calculate and visualize sequence coverage in alignment files in command line

- [bamtocov](https://github.com/telatin/bamtocov) - coverage extraction from BAM/CRAM files to wig format

- [bamcount](https://github.com/BenLangmead/bamcount) - BigWig and BAM utilities, coverage, by Ben Langmead

- [bigwig-nim](https://github.com/brentp/bigwig-nim) - single static binary + liberal license of tool to convert bed to bigwig (and back) and get fast coverage stats from bigwig, by Brent Pedersen, [Twitter](https://twitter.com/brent_p/status/1162765386548305923?s=03)

- [cgpBigWig](https://github.com/cancerit/cgpBigWig) - Package of C scripts for generation of BigWig coverage files. 

- [covviz](https://github.com/brwnj/covviz) - calculate and view coverage based variation

- [indexcov](https://github.com/brentp/goleft/tree/master/indexcov) - Quickly estimate coverage from a whole-genome bam or cram index

- [fastq-pair](https://github.com/linsalrob/fastq-pair) - Match up paired end fastq files quickly and efficiently

- [FastUniq](https://sourceforge.net/projects/fastuniq/) - an ultrafast de novo duplicates removal tool for paired short DNA sequences.

- [fastqwiper](https://github.com/mazzalab/fastqwiper) - An ensemble method to recover corrupted FASTQ files, drop or fix pesky lines, remove unpaired reads, and settle reads interleaving.

- [faster](https://github.com/angelovangel/faster) - A (very) fast program for getting statistics about a fastq file, written in Rust. Get the read length, GC content, mean Phred scores, trim frong and tail, regex search. Compiled binaries are available

- [fgbio tools](https://fulcrumgenomics.github.io/fgbio/tools/latest/) - general purpose tools for working with sequencing data. Examples of tools for manipulating FASTA/FASTQ files (FastqToBam, TrimFastq), RNA-seq (EstimateRnaSeqInsertSize), SAM/BAM (FilterBam, FindTechnicalReads, SplitBam, TrimPrimers), UMI (CollectDuplexSeqMetrics), VCF/BCF (AssessPhasing, FilterSomaticVcf, FixVcfPhaseSet). [GitHub](https://github.com/fulcrumgenomics/fgbio)

- [rasusa](https://github.com/mbhall88/rasusa) - Randomly subsample sequencing reads to a specified coverage, single- and paired end reads

- [RSEM tutorial](https://github.com/bli25/RSEM_tutorial) - A short tutorial on how to use RSEM

### Aligners

- [Chromap](https://github.com/haowenz/chromap) - ultra-fast aligner (>10X faster) for ChIP-seq, Hi-C, scATAC-seq. Based on the minimizer sketch. Memory depends only on genome index size, \~20Gb for human. <details>
    <summary>Paper</summary>
    Zhang, Haowen, Li Song, Xiaotao Wang, Haoyu Cheng, Chenfei Wang, Clifford A. Meyer, Tao Liu, et al. “Fast Alignment and Preprocessing of Chromatin Profiles with Chromap.” Nature Communications, 12 November 2021, https://doi.org/10.1038/s41467-021-26865-w
</details>

- [LJA](https://github.com/AntonBankevich/LJA) (La Jolla Assembler) - long-read genome assembler using multiplex de Bruijn graphs. Utilizes three modules/ideas: jumboDBG (constructing large de Bruijn graphs), mowerDBG (error-correcting reads), and multiplexDBG (utilizing the entire read-length for resolving repeats). Includes the LJApolish module that expands the collapsed homopolymer runs in the resulting assembly. Benchmarked against hifiasm and HiCanu, significantly reduced the number of misassemblies. <details>
    <summary>Paper</summary>
    Bankevich, Anton, Andrey Bzikadze, Mikhail Kolmogorov, Dmitry Antipov, and Pavel A. Pevzner. “LJA: Assembling Long and Accurate Reads Using Multiplex de Bruijn Graphs.” Preprint. Bioinformatics, December 11, 2020. https://doi.org/10.1101/2020.12.10.420448.
</details>

- [SNAP](https://github.com/amplab/snap) - paired-read short-read (70-300bp) aligner based on fussy set intersection. 2-5x faster than BWA-mem2, Bowtie2.  When used with Haplotype Caller from the Genome Analysis Toolkit, SNAP produces better concordance with known-truth sets than other aligners for most of the genome-in-a-bottle and Illumina Platinum genomes. Additonal features: accepts SAM and BAM, outputs sorted, duplicate marked and indexed file. Binaries for Windows, Mac, Linux. [Tweet](https://twitter.com/razoralign/status/1465664608711036932?s=20). <details>
    <summary>Paper</summary>
    Bolosky, William J., Arun Subramaniyan, Matei Zaharia, Ravi Pandya, Taylor Sittler, and David Patterson. “Fuzzy Set Intersection Based Paired-End Short-Read Alignment.” Preprint. Bioinformatics, November 23, 2021. https://doi.org/10.1101/2021.11.23.469039.
</details>

#### Long-read

- [Verkko](https://github.com/marbl/verkko) - long-read diploid genome assembly, for PacBio HiFi and Oxford Nanopore reads. Telomere-to-telomere phased, diploid genome assembly. <details>
    <summary>Paper</summary>
    Rautiainen, Mikko, Sergey Nurk, Brian P. Walenz, Glennis A. Logsdon, David Porubsky, Arang Rhie, Evan E. Eichler, Adam M. Phillippy, and Sergey Koren. “Telomere-to-Telomere Assembly of Diploid Chromosomes with Verkko.” Nature Biotechnology, February 16, 2023. https://doi.org/10.1038/s41587-023-01662-6.
</details>

- [Minimap2](https://github.com/lh3/minimap2) - aligner for long- (SMRT, ONT technologies, over 1kb) and short- (over 100bp, paired-end supported) reads. Spli-read alignment, gap cost for long insertions and deletions, reduces spurious alignment. 3-4 tiimes faster than short-read aligners (C and Python implementation), over 30 times faster than long-read aligners (BLASR, BWA-MEM, GraphMap, minialign, NGMLR). Presets of parameters. <details>
    <summary>Paper</summary>
    Li, Heng. “Minimap2: Pairwise Alignment for Nucleotide Sequences.” Edited by Inanc Birol. Bioinformatics 34, no. 18 (September 15, 2018): 3094–3100. https://doi.org/10.1093/bioinformatics/bty191.
</details>

- [LongRead_tutorials](https://github.com/timkahlke/LongRead_tutorials) - Workflows and tutorials for LongRead analysis with specific focus on Oxford Nanopore data. [Website](https://timkahlke.github.io/LongRead_tutorials/)

- [lorax](https://github.com/tobiasrausch/lorax) - A long-read analysis toolbox for cancer genomics applications. Requires matched tumor-normal data sequenced using long-reads.

- [NGMLR](https://github.com/philres/ngmlr) - long-read mapper designed to align PacBio or Oxford Nanopore (standard and ultra-long) to a reference genome with a focus on reads that span structural variations

- [Sniffles](https://github.com/fritzsedlazeck/Sniffles) - structural variation caller using third generation sequencing (PacBio or Oxford Nanopore).


### Analysis

- [GeneTonic](https://bioconductor.org/packages/GeneTonic) - R package for visualization and interpretation of differential expression and functional enrichment analyses results. Heatmaps, volcano plots, PCA, KEGG pathways with gene expression overlay. Shiny app, tutorials, HTML report. Input: differentially expressed genes (DESeq2 preferred), enrichment results (clusterProfiler, Enrichr, others), the `shaker` function processes many data formats. [Alternatives](https://docs.google.com/spreadsheets/d/167XV0w18P0FSld1dt6owN4C2Esxl5FU2QTo4D-wclz0/edit#gid=0). [Code to reproduce the paper](https://github.com/federicomarini/GeneTonic_supplement). [GitHub](https://github.com/federicomarini/GeneTonic). [Demo webserver](http://shiny.imbei.uni-mainz.de:3838/GeneTonic/) <details>
    <summary>Paper</summary>
    Marini, Federico, Annekathrin Ludt, Jan Linke, and Konstantin Strauch. “GeneTonic: An R/Bioconductor Package for Streamlining the Interpretation of RNA-Seq Data.” BMC Bioinformatics 22, no. 1 (December 2021): 610. https://doi.org/10.1186/s12859-021-04461-5.
</details>

- [Omics Playground](https://github.com/bigomics/omicsplayground) - self-service bioinformatics platform to analyze, visualize, and integrate omics data. Descriptive statistics, differential gene expression, clustering, GSEA/KEGG, signature, cell type prediction. Input - text-based files. R/Shiny implementation, Docker image. [GitHub](https://github.com/bigomics/omicsplayground), [Documentation](https://omicsplayground.readthedocs.io/)
    - Akhmedov, Murodzhon, Axel Martinelli, Roger Geiger, and Ivo Kwee. “[Omics Playground: A Comprehensive Self-Service Platform for Visualization, Analytics and Exploration of Big Omics Data](https://doi.org/10.1093/nargab/lqz019).” NAR Genomics and Bioinformatics 2, no. 1 (March 1, 2020)

- [Shiny-Seq](https://github.com/schultzelab/Shiny-Seq) - Shiny app, and Docker image. Input - a count table, or kallisto output, and annotations. Normalization (DESeq2), batch effect removal (limma or SVA),  differential expression analysis (DeSeq2), co-expression network analysis (WGCNA), functional enrichment analysis (clusterprofiler), TFBS motif overrepresentation (pcaGopromoter). Visualization (heatmaps, volcano plots). [GitHub](https://github.com/schultzelab/Shiny-Seq). [RNA-seq blog post](https://www.rna-seqblog.com/shiny-seq-advanced-guided-transcriptome-analysis/)
    - Sundararajan, Zenitha, Rainer Knoll, Peter Hombach, Matthias Becker, Joachim L. Schultze, and Thomas Ulas. “[Shiny-Seq: Advanced Guided Transcriptome Analysis](https://doi.org/10.1186/s13104-019-4471-1).” BMC Research Notes 12, no. 1 (December 2019). 

- [ARMOR](https://github.com/csoneson/ARMOR) - modular RNA-seq analysis, from raw files, QC, alignment, quantification, differential expression and transcript usage analyses, enrichment analysis, shiny app data exploration, comparible with iSEE (Figure 1). Implemented using Snakemake and Conda for reproducibility. Input - FASTQ files and a metadata file. config.yaml - workflow settings. [Data example](https://github.com/csoneson/ARMOR/tree/chiron_realdataworkflow/E-MTAB-7029). <details>
    <summary>Paper</summary>
    Orjuela, Stephany, Ruizhu Huang, Katharina M. Hembach, Mark D. Robinson, and Charlotte Soneson. “ARMOR: An Automated Reproducible MOdular Workflow for Preprocessing and Differential Analysis of RNA-Seq Data.” G3&amp;#58; Genes|Genomes|Genetics, May 14, 2019, g3.400185.2019. https://doi.org/10.1534/g3.119.400185.
</details>

- [3D_RNA-seq](https://github.com/wyguo/ThreeDRNAseq) - R package and Shiny app, and Docker image, for differential expression, differential alternative splicing, and differential Transcript Usage. Two- and mutliple group analysis, time course. Input - Salmon/Kallisto transcript quantification files, or .csv. Diagnostic plots, PCA, batch removal using RUVseq, limma-voom for differential expression, iso-kTSP and TSIS for isoform switching between groups and in time course, respectively. Heatmaps, barplots, volcano plots, Venn diagrams. 3D stands for **D**ifferential Expression, **D**ifferential Alternative Splicing, and **D**ifferential Transcript Usage analyses. [Web](https://ics.hutton.ac.uk/3drnaseq/),  [Manual](https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md),  [GitHub](https://github.com/wyguo/ThreeDRNAseq)
    - Guo, Wenbin, Nikoleta Tzioutziou, Gordon Stephen, Iain Milne, Cristiane Calixto, Robbie Waugh, John WS Brown, and Runxuan Zhang. “[3D RNA-Seq - a Powerful and Flexible Tool for Rapid and Accurate Differential Expression and Alternative Splicing Analysis of RNA-Seq Data for Biologists](https://doi.org/10.1101/656686).” Preprint. Bioinformatics, May 31, 2019.

- [DrEdGE](https://github.com/ptgolden/dredge) - Differential Expression Gene Explorer, Takes in a table of transcript abundance counts, experimental design, other input can be generated in R. Plot, table, heatmap visualization options. Web: http://dredge.bio.unc.edu/ , GitHub: https://github.com/ptgolden/dredge
    - Tintori, Sophia C, Patrick Golden, and Bob Goldstein. “[Differential Expression Gene Explorer (DrEdGE): A Tool for Generating Interactive Online Data Visualizations for Exploration of Quantitative Transcript Abundance Datasets](https://doi.org/10.1101/618439).” Preprint. Genomics, April 25, 2019. 

- [Phantasus](https://github.com/ctlab/phantasus) - interactive exploratory analyses of genomic data, from clustering, PCA, to enrichment, network, and pathway analyses. Works on user data, ARCHS4, TCGA. https://github.com/ctlab/phantasus. [Bioconductor R package](http://www.bioconductor.org/packages/release/bioc/html/phantasus.html) and [Docker image](https://hub.docker.com/r/dzenkova/phantasus)

## Quality control

- [Qualimap](http://qualimap.conesalab.org/) - Qualimap 2 is a platform-independent application written in Java and R that provides both a Graphical User Inteface (GUI) and a command-line interface to facilitate the quality control of alignment sequencing data and its derivatives like feature counts. Supported types of experiments include: Whole-genome sequencing, Whole-exome sequencing, RNA-seq (speical mode available), ChIP-seq

- [fastp](https://github.com/OpenGene/fastp) - fast C++ parallelized tool for FASTQ quality control, adapter trimming, quality filtering, pruning, polyX (polyG) trimming, works with single- and paired-end data. [GitHub](https://github.com/OpenGene/fastp)
    - Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. “[Fastp: An Ultra-Fast All-in-One FASTQ Preprocessor](https://doi.org/10.1093/bioinformatics/bty560).” Bioinformatics 34, no. 17 (September 1, 2018)

- [MultiQC](http://multiqc.info/) - Summarization and visualization QC results for multiple samples in one report. Recognizes multiple QC tools

- [ngsReports](https://github.com/UofABioinformaticsHub/ngsReports) - An R Package for managing FastQC reports and other NGS related log files. [biorXiv](http://biorxiv.org/content/early/2018/05/02/313148.abstract)

- [sickle](https://github.com/najoshi/sickle) - A windowed adaptive trimming tool for FASTQ files using quality. Post-adapter trimming step

- [fastqc](https://github.com/Malarkey73/fastqc) - an R package for quality control (QC) of short read fastq files, analog of the [original FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- [FastQt](https://github.com/labsquare/FastQt) - FastQC port to Qt5: A quality control tool for high throughput sequence data

- [fastqcheck](https://github.com/VertebrateResequencing/fastqcheck) - Generate statistics on and validate fastq files

## Imputation

- [fancyimpute](https://github.com/iskandr/fancyimpute) - Multivariate imputation and matrix completion algorithms implemented in Python. Algorithms: SimpleFill, KNN, SoftImpute, IterativeImputer, IterativeSVD, MatrixFactorization, NuclearNormMinimization, BiScaler.

- [softImpute](https://CRAN.R-project.org/package=softImpute) - R package for Matrix Completion via Iterative Soft-Thresholded SVD, by Trevor Hastie and Rahul Mazumder.

## Deconvolution

See also [Cancer_notes/Deconvolution](https://github.com/mdozmorov/Cancer_notes#deconvolution)

- [AutoGeneS](https://github.com/theislab/AutoGeneS) - Python code for RNA-seq deconvolution. Novel feature selection method (multi-object optimization) eliminating collinear genes and improving deconvolution accuracy. Signature matrix generation from annotated scRNA-seq data, as well as using sorted/purified cell signatures. Nu-Support Vector Regression (Nu-SVR) better controlling for outliers. Compatible with the scanpy pipeline. [Documentation](https://autogenes.readthedocs.io/en/latest/). <details>
    <summary>Paper</summary>
    Aliee, Hananeh, and Fabian J. Theis. "AutoGeneS: Automatic gene selection using multi-objective optimization for RNA-seq deconvolution." Cell Systems 12, no. 7 (2021): 706-715. https://doi.org/10.1016/j.cels.2021.05.006
</details>

- [DeCompress](https://github.com/bhattacharya-a-bt/DeCompress) - R package for a semi-reference-free deconvolution method for targeted panels. Input: a target gene expression panel (RNA-seq or microarray, raw scale), and a reference dataset from similar tissue (not compartment-specific profiles) to expand the feature space of targeted panels using compressed sensing. Output: cell type proportions and gene expression. Ensemble reference-free deconvolution is performed on this artificially expanded dataset to estimate cell-type proportions and gene signatures. Three steps: 1) selection of compartment-specific genes from the reference, 2) compressed sensing to expand the targeted panel, and 3) ensemble deconvolution of the expanded dataset. Benchmarked against deconf, CellDistinguisher, Linseed, DeconICA, TOAST, on in-silico mixtures (GTeX), published mixing experiments, the Carolina Breast Cancer Study. Survival, eQTL analyses demonstrating increased biological relevance. [Scripts to reproduce the results](https://github.com/bhattacharya-a-bt/DeCompress_supplement). <details>
    <summary>Paper</summary>
    Bhattacharya, Arjun, Alina M. Hamilton, Melissa A. Troester, and Michael I. Love. “DeCompress: Tissue Compartment Deconvolution of Targeted MRNA Expression Panels Using Compressed Sensing.” Preprint. Bioinformatics, August 14, 2020. https://doi.org/10.1101/2020.08.14.250902.
</details>

- [GEDIT](https://github.com/BNadel/GEDIT) - cell type deconvolution from gene expression data. Supports microarray, RNA-seq data, blood and stromal tissue types, human and mouse species. Includes reference data from 8 sources (Human Skin Signatures, Human Body Atlas, Human Primary Cell Atlas, Blueprint, Encode, LM22, 10X single cell dataset, Immunostates, Tabula Muris, Mouse Body Atlas, ImmGen). Input - gene expression matrix and reference data. Joint quantile normalization, signature score calculation, expression range normalization, non-negative linear regression with no regularization. Compared with CIBERSORT, dtangle, DeconRNASeq, performs better in most cases. <details>
    <summary>Paper</summary>
    Nadel, Brian B, David Lopez, Dennis J Montoya, Feiyang Ma, Hannah Waddel, Misha M Khan, Serghei Mangul, and Matteo Pellegrini. “The Gene Expression Deconvolution Interactive Tool (GEDIT): Accurate Cell Type Quantification from Gene Expression Data.” GigaScience 10, no. 2 (January 29, 2021): giab002. https://doi.org/10.1093/gigascience/giab002.
</details>

- [MuSiC](https://github.com/xuranw/MuSiC) - cell type deconvolution method to use scRNA-seq data (pre-determined cell types) to deconvolve bulk RNA-seq data. Gene weighting to prioritize stable and reliably expressed genes to build cell signatures. Compared with CIBERSORT, Nonnegative least squares (NNLS), BSEQ-sc. <details>
    <summary>Paper</summary>
    Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao Li. “Bulk Tissue Cell Type Deconvolution with Multi-Subject Single-Cell Expression Reference.” Nature Communications 10, no. 1 (December 2019). https://doi.org/10.1038/s41467-018-08023-x.
</details>

- [SCDC](https://meichendong.github.io/SCDC/) - Bulk Gene Expression Deconvolution by Integrating Multiple Single-Cell RNA Sequencing References. Contrasted with other methods (Bseq-SC, DWLS, MuSiC, Bisque, CIBERSORTx) that utilize one scRNA-seq dataset. Outperforms existing methods even if using one single-cell dataset. R package. <details>
    <summary>Paper</summary>
    Dong, Meichen, Aatish Thennavan, Eugene Urrutia, Yun Li, Charles M. Perou, Fei Zou, and Yuchao Jiang. “SCDC: Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA Sequencing References.” Preprint. Bioinformatics, August 22, 2019. https://doi.org/10.1101/743591.
</details>

## Batch effect

- [ComBat-seq](https://github.com/zhangyuqing/sva-devel) - batch effect correction for RNA-seq data using negative binomial regression. Maintains count nature of RNA-seq data. Tested on simulated data (polyester package), and experimental data. Achieves the highest true positive rate. [Code to reproduce paper](https://github.com/zhangyuqing/ComBat-seq)
    - Zhang, Yuqing, Giovanni Parmigiani, and W. Evan Johnson. “[ComBat-Seq: Batch Effect Adjustment for RNA-Seq Count Data](https://doi.org/10.1101/2020.01.13.904730).” Preprint. Bioinformatics, January 14, 2020. 

- [combat.py](https://github.com/brentp/combat.py) - python / numpy / pandas / patsy version of ComBat for removing batch effects


## Clustering

- [clust](https://github.com/BaselAbujamous/clust) - Python package for identification of consistently coexpressed clusters of genes, within and across datasets. Consensus clustering principle. Aproximately 50% of genes do not cluster well and thus shouldn't be considered. Compared with seven tools (cross-clustering, k-means, SOMs, MCL, HC, Click, WGCNA) using seven different cluster validation metrics. Outperforms all, produces more focused and significant functional enrichment results. 
    - Abu-Jamous, Basel, and Steven Kelly. “[Clust: Automatic Extraction of Optimal Co-Expressed Gene Clusters from Gene Expression Data](https://doi.org/10.1186/s13059-018-1536-8).” Genome Biology 19, no. 1 (December 2018)

## Timecourse

- Benchmarking of time course differential analysis tools for RNA-seq. Classical pairwise comparison outperforms specially designed methods in terms of overall performance and robustness to noise. Tested on stimulated data (generated using parameters estimated from ReCount data) and experimental data ([GSE69822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69822)). AdaptiveGP, DyNB, EbSeqHMM, edgeR/DESeq2 (used as reference gold standard), FunPat, ImpulseDE2, lmms, Next maSigPro, splineTimeR, TimeSeq. Brief description of each method. ImpulseDE2 - overall best performing, comparable with pairwise comparison. Followed by SplineTC on long time series, lmms and others (Discussion). Despite differences, functional enrichment results are similar. [GitHub](https://github.com/daniel-spies/rna-seq_tcComp) with data and analysis scripts. <details>
    <summary>Paper</summary>
    Spies, Daniel, Peter F Renz, Tobias A Beyer, and Constance Ciaudo. “Comparative Analysis of Differential Gene Expression Tools for RNA Sequencing Time Course Data.” Briefings in Bioinformatics 20, no. 1 (January 18, 2019): 288–98. https://doi.org/10.1093/bib/bbx115.
</details>

- [LPWC](https://gitter-lab.github.io/LPWC/) - Lag Penalized Weighted Correlation, a similarity measure to group pairs of time series that are not perfectly synchronized. Review of previous approaches (hierarchical clustering, partition-based, Bayesian models). Correlation-based, with the lag penalty for shift, two options to select it. Best for 5 or more time points, for shorter time course - either use high penalty or another tool STEM. Tested on simulated data (ImpulsDE data) using the adjusted Rand index
    - Chandereng, Thevaa, and Anthony Gitter. “[Lag Penalized Weighted Correlation for Time Series Clustering](https://doi.org/10.1186/s12859-019-3324-1).” BMC Bioinformatics 21, no. 1 (December 2020)

- Time course gene expression analysis review. Biological scenarios requiring a time-course, analytical approaches, Table 1 - software for time course analysis (EDGE, BETR, clustering tools, network analysis).
    - Bar-Joseph, Ziv, Anthony Gitter, and Itamar Simon. “[Studying and Modelling Dynamic Biological Processes Using Time-Series Gene Expression Data](https://doi.org/10.1038/nrg3244).” Nature Reviews. Genetics 13, no. 8 (July 18, 2012)

- [DREM 2.0](http://sb.cs.cmu.edu/drem/) - time course analysis of gene expression data. Detects patterns of gene expression changes and the corresponding transcription factors driving them, motif discovery using protein-DNA (ChIP-seq, ChIP-chip, computational) data, differential motif analysis (DECOD method). Hidden Markov Model-based algorithm. Java tool, GUI and command line interface. 
    - Schulz, Marcel H, William E Devanny, Anthony Gitter, Shan Zhong, Jason Ernst, and Ziv Bar-Joseph. “[DREM 2.0: Improved Reconstruction of Dynamic Regulatory Networks from Time-Series Expression Data](https://doi.org/10.1186/1752-0509-6-104).” BMC Systems Biology 6, no. 1 (2012)


## Differential expression

- [lmerSeq](https://github.com/stop-pre16/lmerSeq) - linear mixed models for RNA-seq data (VST, variance stabilization transformed). Wraps lme4, lmerTest, nlme packages. Tested on synthetic and experimental data vs. DREAM, rmRNA-seq, almost uniformly more powerful, has the best control of FDR. [GitHub](https://github.com/stop-pre16/lmerSeq) with tutorials. <details>
    <summary>Paper</summary>
    Vestal, Brian E., Elizabeth Wynn, and Camille M. Moore. “LmerSeq: An R Package for Analyzing Transformed RNA-Seq Data with Linear Mixed Effects Models.” BMC Bioinformatics 23, no. 1 (November 16, 2022): 489. https://doi.org/10.1186/s12859-022-05019-9.
</details>

- [Glimma 2.0](https://github.com/hasaru-k/GlimmaV2) - R package for interactive visualization of DESeq2, edgeR, Limma objects. MDS plot, MA plot, Volcano plot. D3, htmlwidgets, plotly, dygraphs. Plots are embeddable in RMarkdown. Export data as CSV, PNG, SVG. [Bioconductor](https://bioconductor.org/packages/Glimma/). <details>
    <summary>Paper</summary>
    Kariyawasam, Hasaru, Shian Su, Oliver Voogd, Matthew E Ritchie, and Charity W Law. “Dashboard-Style Interactive Plots for RNA-Seq Analysis Are R Markdown Ready with Glimma 2.0.” NAR Genomics and Bioinformatics 3, no. 4 (October 4, 2021): lqab116. https://doi.org/10.1093/nargab/lqab116.
</details>

- [Topconfect](https://bioconductor.org/packages/release/bioc/html/topconfects.html) - differential gene expression using confidence intervals of log fold changes. Confect units are interpretable as effect sizes. Uses TREAT (limma) functionality, Methods, also wraps edgeR and DESeq2 functionality. Tested using simulated and tumor-normal data. R package, [GitHub](https://github.com/pfh/topconfects), [Bioconductor](https://bioconductor.org/packages/release/bioc/html/topconfects.html)
    - Harrison, Paul F., Andrew D. Pattison, David R. Powell, and Traude H. Beilharz. “[Topconfects: A Package for Confident Effect Sizes in Differential Expression Analysis Provides a More Biologically Useful Ranked Gene List](https://doi.org/10.1186/s13059-019-1674-7).” Genome Biology 20, no. 1 (December 2019)

- [Degust](http://degust.erc.monash.edu/) -  interactive RNA-seq analysis, RNA-seq exploration, analysis and visualisation

- [diffexpr](https://github.com/wckdouglas/diffexpr) - A python package using rpy2 to port DESeq2 into Python.

- [ideal](http://bioconductor.org/packages/release/bioc/html/ideal.html) - Interactive Differential Expression AnaLysis. 

- [EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano) - Publication-ready volcano plots with enhanced colouring and labeling, by Kevin Blighe

- [ALDEx](https://bioconductor.org/packages/release/bioc/html/ALDEx2.html) - differential abundance analysis of compositional data. ANOVA-like approach, partitioning within-condition to between-condition variation. Uses a Dirichlet-multinomial model to infer abundance from counts. Methods description for transforming proportional data into independent components
    - Fernandes, Andrew D., Jean M. Macklaim, Thomas G. Linn, Gregor Reid, and Gregory B. Gloor. “[ANOVA-like Differential Expression (ALDEx) Analysis for Mixed Population RNA-Seq](https://doi.org/10.1371/journal.pone.0067019).” PloS One 8, no. 7 (2013)

## Allele-specific expression

- [ASEP](https://github.com/Jiaxin-Fan/ASEP) - allele-specific (differential) gene expression in a population from RNA-seq. A generalized linear mixed-effects model, the subject-specific random effect accounts for correlation of multiple SNPs within the same gene, and sample heterogeneity is modeled by a two-component mixture distribution. A pseudo-phasing approach (MBASED), 'majority voting'. Type I error is approx. 1%. Input - gene-specific SNP alleles and read counts. <details>
    <summary>Paper</summary>
    Fan, Jiaxin, Jian Hu, Chenyi Xue, Hanrui Zhang, Katalin Susztak, Muredach P. Reilly, Rui Xiao, and Mingyao Li. “ASEP: Gene-Based Detection of Allele-Specific Expression across Individuals in a Population by RNA Sequencing.” Edited by Xiaofeng Zhu. PLOS Genetics 16, no. 5 (May 11, 2020): e1008786. https://doi.org/10.1371/journal.pgen.1008786.
</details>

- [ASEQ](https://bitbucket.org/CibioBCG/aseq/src/master/) - allele-specific gene/transcript expression analysis from paired genomic and transcriptomic sequencing. Does not require paternal and maternal genome data. Concordant with AlleleSeq and MBASED tools. Fast, C implementation. <details>
    <summary>Paper</summary>
    Romanel, Alessandro, Sara Lago, Davide Prandi, Andrea Sboner, and Francesca Demichelis. “ASEQ: Fast Allele-Specific Studies from next-Generation Sequencing Data.” BMC Medical Genomics 8, no. 1 (December 2015): 9. https://doi.org/10.1186/s12920-015-0084-2.
</details>

- [SNPsplit](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/) - allele-specific splitting of SAM/BAM alignments using known SNP genotypes. Based on alignment to SNP-masked genomes (N). Perl, three scripts (SNPsplit, SNPsplit_genome_preparation, tag2sort). Input - SAM/BAM files, an annotation file containing the positions of all SNPs in the genome. Performs 1) read tagging and 2) read sorting. Classifies into Allele 1-specific, Allele 2-specific, Unassigned, Conflicting. Works with DNA-seq, RNA-seq, Hi-C, Bisulfite-seq data, single- and paired end, aligned by various aligners. [GitHub](https://github.com/FelixKrueger/SNPsplit). <details>
    <summary>Paper</summary>
    Krueger, Felix, and Simon R. Andrews. “SNPsplit: Allele-Specific Splitting of Alignments between Genomes with Known SNP Genotypes.” F1000Research 5 (July 2016): 1479. https://doi.org/10.12688/f1000research.9037.2.
</details>

- [phASER](https://github.com/secastel/phaser) - phasing and Allele Specific Expression from RNA-seq. Performs haplotype phasing using read alignments in BAM format from both DNA and RNA based assays, and provides measures of haplotypic expression for RNA based assays.

## Pathways, Functional enrichment

- [pathDIP](http://ophid.utoronto.ca/pathDIP/) - **path**way **D**ata **I**ntegration **P**ortal integrating 24 major databases (5366 for human), adds protein-protein interactions to pathways. Enrichment analysis. API for Java, Python, R. <details>
    <summary>Paper</summary>
    Rahmati, Sara, Mark Abovsky, Chiara Pastrello, et al. “pathDIP 4: An Extended Pathway Annotations and Enrichment Analysis Resource for Human, Model Organisms and Domesticated Species.” Nucleic Acids Research, November 16, 2019, gkz989. https://doi.org/10.1093/nar/gkz989.
    Pastrello C, Kotlyar M, Abovsky M, Lu R, Jurisica I. PathDIP 5: improving coverage and making enrichment analysis more biologically meaningful. Nucleic Acids Res. 52(D1):D663-D671, 2024. 
</details>

- [rGREAT](https://bioconductor.org/packages/rGREAT/) - GREAT method for TSS-centric enrichment of genomic regions. Interface with online GREAT, local implementation. Integrates GO, MSigDB, supports more than 600 organisms via [BioMartGOGeneSets](https://github.com/jokergoo/BioMartGOGeneSets), custom organisms/annotations support. Results viewable via Shiny app. Different TSS annotations, despite differences, produce similar results. [GitHub](https://github.com/jokergoo/rGREAT). <details>
    <summary>Paper</summary>
    Gu, Zuguang, and Daniel Hübschmann. “RGREAT: An R/Bioconductor Package for Functional Enrichment on Genomic Regions.” Edited by Tobias Marschall. Bioinformatics, November 17, 2022, btac745. https://doi.org/10.1093/bioinformatics/btac745.
</details>

- [eVITTA](https://tau.cmmt.ubc.ca/eVITTA/) - web tool for interactive gene expression and functional enrichment (GSEA, overrepresentation) analyses. Three modules: 1) easyGEO, retrieval and analysis of GEO datasets; 2) easyGSEA, GSEA or overrepresentation enrichment analyses; 3) easyVizR, comparison among experimental groups (overlap of gene lists and enrichment results). Input: DESeq2/edgeR output, (ranked) gene lists. Output: interactive barplots, heatmaps, volcano plots (plotly), rank-rank hypergeometric overlap (RRHO) plot, networks, enrichment tables. Figure 1, Table 1 - overview of analyses, inputs, outputs. [GitHub](https://github.com/easygsea/eVITTA). <details>
    <summary>Paper</summary>
    Cheng, Xuanjin, Junran Yan, Yongxing Liu, Jiahe Wang, and Stefan Taubert. “EVITTA: A Web-Based Visualization and Inference Toolbox for Transcriptome Analysis.” Nucleic Acids Research 49, no. W1 (July 2, 2021): W207–15. https://doi.org/10.1093/nar/gkab366.
</details>

- [GSEABenchmarking](https://github.com/waldronlab/GSEABenchmarking) - an R package for systematic testing of gene set enrichment analyses. 10 major enrichment analyses tested on runtime, % significant sets, type I error rate, relevance to phenotype. Microarray and TCGA RNA-seq data. Best performing - overrepresentation analysis, aka Fisher's, hypergeometric test. [Tweet by Levi Waldron](https://twitter.com/LeviWaldron1/status/1142092301403115521?s=03)
    - Geistlinger, Ludwig, Gergely Csaba, Mara Santarelli, Marcel Ramos, Lucas Schiffer, Charity W Law, Nitesh Turaga, et al. “[Towards a Gold Standard for Benchmarking Gene Set Enrichment Analysis](https://doi.org/10.1101/674267).” Preprint. Bioinformatics, June 19, 2019. 

- [PaintOmics 3](http://www.paintomics.org/) - web tool for KEGG pathway enrichment analysis and visualization of gene expression (also, metabolite, protein, region-based data) over pathway diagrams. Competitors: MapMan, KaPPA-View, Pathview Web. Auto-detection of IDs. Analyzes fold change, time course
    - Hernández-de-Diego, Rafael, Sonia Tarazona, Carlos Martínez-Mira, Leandro Balzano-Nogueira, Pedro Furió-Tarí, Georgios J. Pappas, and Ana Conesa. “[PaintOmics 3: A Web Resource for the Pathway Analysis and Visualization of Multi-Omics Data](https://doi.org/10.1093/nar/gky466).” Nucleic Acids Research 46, no. W1 (July 2, 2018)

- [singscore](https://bioconductor.org/packages/singscore/) - a rank-based single-sample enrichment of molecular signatures (R/Bioconductor). Does not depend on background samples, stable for the analysis of small number of samples (<25). Compared woth GSVA, Z-score, PLAGE, ssGSEA, had overall better power, fast. [PySingscore](https://github.com/DavisLaboratory/PySingscore) - Python implementation. [GitHub](https://github.com/DavisLaboratory/singscore). <details>
    <summary>Paper</summary>
    Foroutan, Momeneh, Dharmesh D. Bhuva, Ruqian Lyu, Kristy Horan, Joseph Cursons, and Melissa J. Davis. “Single Sample Scoring of Molecular Phenotypes.” BMC Bioinformatics 19, no. 1 (December 2018): 404. https://doi.org/10.1186/s12859-018-2435-4.
</details>

- [CEMiTool](https://bioconductor.org/packages/release/bioc/html/CEMiTool.html) - gene co-expression analysis, reimplements WGCNA, includes selection of a soft-thresholding power using Cauchi distribution, gene enrichment analysis and, optionally, PPI network. Good overview of WGCNA algorithm. 
    - Russo, Pedro S. T., Gustavo R. Ferreira, Lucas E. Cardozo, Matheus C. Bürger, Raul Arias-Carrasco, Sandra R. Maruyama, Thiago D. C. Hirata, et al. “[CEMiTool: A Bioconductor Package for Performing Comprehensive Modular Co-Expression Analyses](https://doi.org/10.1186/s12859-018-2053-1).” BMC Bioinformatics 19, no. 1 (20 2018)

- [EnrichmentBrowser](https://bioconductor.org/packages/release/bioc/html/EnrichmentBrowser.html) - R package for microarray/RNA-seq normalization, ID mapping, differential analysis, functional enrichment (many methods) and network analyses and visualization
     - Geistlinger, Ludwig, Gergely Csaba, and Ralf Zimmer. “[Bioconductor’s EnrichmentBrowser: Seamless Navigation through Combined Results of Set- & Network-Based Enrichment Analysis](https://doi.org/10.1186/s12859-016-0884-1).” BMC Bioinformatics 17 (January 20, 2016)

- [Metascape](http://metascape.org/gp/index.html#/main/step1) - multi-gene lists functional analysis, auto gene ID recognition, >40 databases. [data_analysis_portals.xlsx]((https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09234-6/MediaObjects/41467_2019_9234_MOESM3_ESM.xlsx)) - 25 data analysis portals
    - Zhou, Yingyao, Bin Zhou, Lars Pache, Max Chang, Alireza Hadj Khodabakhshi, Olga Tanaseichuk, Christopher Benner, and Sumit K. Chanda. “[Metascape Provides a Biologist-Oriented Resource for the Analysis of Systems-Level Datasets](https://doi.org/10.1038/s41467-019-09234-6).” Nature Communications 10, no. 1 (December 2019)

- [Positional Gene Enrichment analysis](http://silico.biotoul.fr/pge/) of gene sets for high resolution identification of overrepresented chromosomal regions.  
    - [De Preter, K., Barriot, R., Speleman, F., Vandesompele, J., Moreau, Y., 2008, Nucleic Acids Res.](https://academic.oup.com/nar/article/36/7/e43/2410597)

- [blitzgsea](https://github.com/MaayanLab/blitzgsea) - Fast Gene Set Enrichment Analysis (GSEA) implementation of the prerank algorithm. Use Loess interpolation of bimodal ES distribution for accurate p-value estimation. Python.

- [gProfileR](https://CRAN.R-project.org/package=gProfileR) - enrichment of gene lists, from GO to KEGG and others, organism-specific, ortholog conversion. [Web interface](https://biit.cs.ut.ee/gprofiler/gost)

- [msigdb4mouse](https://github.com/markziemann/msigdb4mouse) - MSigDB gene sets for mouse (GO, KEGG, Reactome).

- [Rummagene](https://rummagene.com/about) - Finding gene sets from supplementary materials matching query. [Youtube video](https://www.youtube.com/watch?v=RftBcfWsSOs), 47m.

### Transcription regulators

- [ChEA3](https://amp.pharm.mssm.edu/chea3/) - predicting regulatory TFs for sets of user-provided genes. Improved backend reference gene set data (six datasets), ranking of the most significantly enriched TFs. Benchmarking against several other TF prioritization tools (overviewed in intro). [Docker](https://hub.docker.com/r/maayanlab/chea3), [API web-interface and downloadable data](https://amp.pharm.mssm.edu/chea3/)
    - Keenan, Alexandra B., Denis Torre, Alexander Lachmann, Ariel K. Leong, Megan L. Wojciechowicz, Vivian Utti, Kathleen M. Jagodnik, Eryk Kropiwnicki, Zichen Wang, and Avi Ma’ayan. “[ChEA3: Transcription Factor Enrichment Analysis by Orthogonal Omics Integration](https://doi.org/10.1093/nar/gkz446).” Nucleic Acids Research, May 22, 2019. 

- [RABIT](http://rabit.dfci.harvard.edu/) - find TFs regulating a list of genes. Integrated ChIP-seq and gene expression data, regression framework. Tested in experimental KO data, tumor-profiling cohorts. 
    - Jiang, Peng, Matthew L. Freedman, Jun S. Liu, and Xiaole Shirley Liu. “[Inference of Transcriptional Regulation in Cancers](https://doi.org/10.1073/pnas.1424272112).” Proceedings of the National Academy of Sciences 112, no. 25 (June 23, 2015)

- [RcisTarget](https://github.com/aertslab/RcisTarget) - finding enriched motifs in cis-regulatory regions in a gene list

## Non-canonical RNAs

- [ITAS](https://github.com/EpiEpiMSU/ITAS) - database of transcript annotation for small RNAs. Filtered, corrected, and integrated transcript annotations for several types of small RNAs (miRNA, piRNA, tRNA, rRNA, tsRNA) and several species (human, mouse, rat, fly, worm). Data from miRBase, piRNAdb, GtRNAdb, UCSC, tRFdb, MINTbase. Compared with [SPORTS](https://github.com/junchaoshi/sports1.1) annotation pipeline, detects more transcripts, differentially expressed. [Processing scripts](https://github.com/EpiEpiMSU/ITAS_scripts). <details>
    <summary>Paper</summary>
    Stupnikov, Alexey, Vitaly Bezuglov, Ivan Skakov, Victoria Shtratnikova, J. Richard Pilsner, Alexander Suvorov, and Oleg Sergeyev. “ITAS: Integrated Transcript Annotation for Small RNA.” Non-Coding RNA 8, no. 3 (May 2, 2022): 30. https://doi.org/10.3390/ncrna8030030.
</details>

- [oddgenes](https://github.com/cmdcolin/oddgenes) - A list of weird gene annotations or things that break bioinformatics assumptions

### Alternative splicing

- Benchmarking review, two types of alternative splicing analysis: differential splicing and differential isoform detection. DESeq2, DEXSeq, Limma and NOISeq perform well overall. [GitHub](https://github.com/gamerino/benchmarkingDiffExprAndSpl)
    - Merino, Gabriela A, Ana Conesa, and Elmer A Fernández. “[A Benchmarking of Workflows for Detecting Differential Splicing and Differential Expression at Isoform Level in Human RNA-Seq Studies](https://doi.org/10.1093/bib/bbx122).” Briefings in Bioinformatics 20, no. 2 (March 25, 2019)

- Differential splicing tool benchmarking. he three different methodological categories: exon-based (DEXSeq, edgeR, JunctionSeq, limma), isoform-based (cuffdiff2, DiffSplice) and event-based methods (dSpliceType, MAJIQ, rMATS, SUPPA). Exon-based methods perform well (MAJIQ, rMATS, DEXSeq). <details>
    <summary>Paper</summary>
    Mehmood, Arfa, Asta Laiho, Mikko S Venäläinen, Aidan J McGlinchey, Ning Wang, and Laura L Elo. “Systematic Evaluation of Differential Splicing Tools for RNA-Seq Studies.” Briefings in Bioinformatics 21, no. 6 (December 1, 2020): 2052–65. https://doi.org/10.1093/bib/bbz126.
</details>

- [ASpli](https://bioconductor.org/packages/ASpli/) - integrating several independent measures of alternative splicing. Bin-level analysis of genes/splice junctions (. edgeR to test for differences. Estimates PSI, PIR, PJU, novel junction. Input: BAM files and a genome annotation file. Bioconductor R package
    - Mancini, Estefania, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, and Ariel Chernomoretz. “[ASpli: Integrative Analysis of Splicing Landscapes through RNA-Seq Assays](https://doi.org/10.1093/bioinformatics/btab141),” Bioinformatics, 02 March 2021

- [SUPPA2](https://github.com/comprna/SUPPA) - differential splicing from RNA-seq analysis (differences in PSI). TPM-quantified transcripts, handles replicates, accounts for sequencing depth, covariates. Clusters alternatively spliced events (DBSCAN, OPTICS). Compared with rMATS, DEXseq, MAJIQ on simulated data, very fast, accurate. [Tutorial based on Salmon quantification](https://github.com/comprna/SUPPA/wiki/SUPPA2-tutorial)
    - Trincado, Juan L., Juan C. Entizne, Gerald Hysenaj, Babita Singh, Miha Skalic, David J. Elliott, and Eduardo Eyras. “[SUPPA2: Fast, Accurate, and Uncertainty-Aware Differential Splicing Analysis across Multiple Conditions](https://doi.org/10.1186/s13059-018-1417-1).” Genome Biology, (December 2018)

- [LeafCutter](https://github.com/davidaknowles/leafcutter) - intron splicing analysis. Identifies variable splicing events from RNA-seq and finds events of high complexity. Does not measure alternative transcription start sites and alternative polyadenylation directly, hence does not require read assembly or inference of isoforms supported by ambigious reads. GTEx reanalysis identified 31.5% unannotated alternatively excised introns. Outperforms Cufflinks2, rMATS, MAJIQ, also by speed and memory requirements. [Visualization Shiny app](https://leafcutter.shinyapps.io/leafviz/). <details>
    <summary>Paper</summary>
    Li, Yang I., David A. Knowles, Jack Humphrey, Alvaro N. Barbeira, Scott P. Dickinson, Hae Kyung Im, and Jonathan K. Pritchard. “Annotation-Free Quantification of RNA Splicing Using LeafCutter.” Nature Genetics 50, no. 1 (January 2018): 151–58. https://doi.org/10.1038/s41588-017-0004-9.
</details>

- [MAJIQ](https://majiq.biociphers.org/) - local splicing variation analysis. Detects canonical and alternative splicing events. Quantifies as Percent Selected In (PSI). Differential splicing as delta PSI. Visualization using VOILA package. Python 3.
    - Vaquero-Garcia, Jorge, Alejandro Barrera, Matthew R. Gazzara, Juan González-Vallinas, Nicholas F. Lahens, John B. Hogenesch, Kristen W. Lynch, and Yoseph Barash. “[A New View of Transcriptome Complexity and Regulation through the Lens of Local Splicing Variations](https://doi.org/10.7554/eLife.11752).” ELife 5 (February 1, 2016)

- [MISO](https://miso.readthedocs.io/en/fastmiso/) (Mixture-of-Isoforms) is a probabilistic framework that quantitates the expression level of alternatively spliced genes from RNA-Seq data, and identifies differentially regulated isoforms or exons across samples. - By modeling the generative process by which reads are produced from isoforms in RNA-Seq, the MISO model uses Bayesian inference to compute the probability that a read originated from a particular isoform.- MISO treats the expression level of a set of isoforms as a random variable and estimates a distribution over the values of this variable. - The estimation algorithm is based on sampling, and falls in the family of techniques known as Markov Chain Monte Carlo (“MCMC”). 
    - Katz, Yarden, Eric T. Wang, Edoardo M. Airoldi, and Christopher B. Burge. “[Analysis and Design of RNA Sequencing Experiments for Identifying Isoform Regulation](https://doi.org/10.1038/nmeth.1528).” Nature Methods 7, no. 12 (December 2010)

- [RegTools](https://regtools.readthedocs.io/en/latest/) - integration of somatic variants from DNA-seq and splice junctions from RNA-seq data to identify variants causing aberrant splicing in cancer. 
    - Feng, Yang-Yang, Avinash Ramu, Kelsy C Cotto, Zachary L Skidmore, Jason Kunisaki, Donald F Conrad, Yiing Lin, et al. “[RegTools: Integrated Analysis of Genomic and Transcriptomic Data for Discovery of Splicing Variants in Cancer](https://doi.org/10.1101/436634),” November 25, 2018. 

- [rMATS](http://rnaseq-mats.sourceforge.net/) alternative splicing detection tool. Using paired samples.RNA-seq depth and alternative splicing power - 200M reads minimum. [rMATS-turbo](https://github.com/Xinglab/rmats-turbo/) - 100X faster implementation, [Tweet](https://twitter.com/YiXing77/status/1267834667715235843?s=20)
    - Shen, Shihao, Juw Won Park, Zhi-xiang Lu, Lan Lin, Michael D. Henry, Ying Nian Wu, Qing Zhou, and Yi Xing. “[RMATS: Robust and Flexible Detection of Differential Alternative Splicing from Replicate RNA-Seq Data](https://doi.org/10.1073/pnas.1419161111).” Proceedings of the National Academy of Sciences of the United States of America 111, no. 51 (December 23, 2014)

- [MMSEQ](https://github.com/eturro/mmseq) - Haplotype, isoform and gene level expression analysis using multi-mapping RNA-seq reads. Including differential analysis. 

- [tappAS](http://tappas.org/) - functional impact of alternative splicing. Input - transcript-level count matrix

- [vast-tools](https://github.com/vastgroup/vast-tools) - A toolset for profiling alternative splicing events in RNA-Seq data
    - Irimia, Manuel, Robert J. Weatheritt, Jonathan D. Ellis, Neelroop N. Parikshak, Thomas Gonatopoulos-Pournatzis, Mariana Babor, Mathieu Quesnel-Vallières, et al. “[A Highly Conserved Program of Neuronal Microexons Is Misregulated in Autistic Brains](https://doi.org/10.1016/j.cell.2014.11.035).” Cell 159, no. 7 (December 18, 2014)

### miRNAs

- miRNA expression atlas across 196 cell type. Reanalysis of 175 studies using miRBase and miRge3.0. UCSCGenome Browser [tracks](https://genome.ucsc.edu/cgi-bin/hgHubConnect). [microRNAome](https://bioconductor.org/packages/microRNAome/) R package, SummarizedExperiment for the microRNAome project. [GitHub](https://github.com/mhalushka/miROme) with code, tutorial, figures, data. <details>
    <summary>Paper</summary>
    

- [CancerMIRNome](http://bioinfo.jialab-ucr.org/CancerMIRNome/) - web server for exploratory miRNA analysis in TCGA cancers and circulating microRNA studies. Query individual miRNAs, cancers. Differential expression, ROC for predicting tumor-normal distinction, survival plots, miRNA-target correlation, functional enrichment of targets. [GitHub](https://github.com/rli012/CancerMIRNome)
    - Li, Ruidong, Han Qu, Shibo Wang, Xuesong Wang, Yanru Cui, Lei Yu, John M. Chater, et al. “[CancerMIRNome: A Web Server for Interactive Analysis and Visualization of Cancer MiRNome Data](https://doi.org/10.1101/2020.10.04.325670).” Preprint. Bioinformatics, October 5, 2020

- [PharmacomiR](http://www.pharmaco-mir.org/) - miRNA-drug associations analysis

- [microRNAome](https://bioconductor.org/packages/microRNAome/) - read counts for microRNAs across tissues, cell-types, and cancer cell-lines, SummarizedExperiment R package

- [miRNAmeConverter](https://bioconductor.org/packages/miRNAmeConverter/) - Convert miRNA Names to Different miRBase Versions

- [MIENTURNET](http://userver.bio.uniroma1.it/apps/mienturnet/) - web tool for miRNA-target enrichment analysis, prioritization, network visualization, functional enrichment for microRNA target genes. 
    - Licursi, Valerio, Federica Conte, Giulia Fiscon, and Paola Paci. “[MIENTURNET: An Interactive Web Tool for MicroRNA-Target Enrichment and Network-Based Analysis](https://doi.org/10.1186/s12859-019-3105-x).” BMC Bioinformatics 20, no. 1 (December 2019)

- [miRDB](http://mirdb.org) - database for miRNA target prediction and functional annotations. The targets were predicted by MirTarget from RNA-seq and CLIP-seq data. Five species: human, mouse, rat, dog and chicken. Custom target prediction. Cell line-specific. Integrative analysis of target prediction and Gene Ontology data. 
    - Chen, Yuhao, and Xiaowei Wang. “[MiRDB: An Online Database for Prediction of Functional MicroRNA Targets](https://doi.org/10.1093/nar/gkz757).” Nucleic Acids Research, August 31, 2019

- [TAM2](http://www.lirmed.com/tam2/) - miRNA enrichment analysis. Manually curated and established miRNA sets. Single list analysis, up vs downregulated. Complementary tools - [miSEA](http://www.baskent.edu.tr/~hogul/misea/), [miEAA](https://ccb-compute2.cs.uni-saarland.de/mieaa2)
    - Li, Jianwei, Xiaofen Han, Yanping Wan, Shan Zhang, Yingshu Zhao, Rui Fan, Qinghua Cui, and Yuan Zhou. “[TAM 2.0: Tool for MicroRNA Set Analysis](https://doi.org/10.1093/nar/gky509).” Nucleic Acids Research 46, no. W1 (July 2, 2018)

- [miRsponge](http://bio-bigdata.hrbmu.edu.cn/miRSponge/) - identification and analysis of miRNA sponge interaction networks and modules. Seven methods for miRNA sponge interaction detection (miRHomology, pc, sppc, hermes, ppc, muTaME, and cernia), and integrative method, description of each method. Four module detection methods (FN, MCL, LINKCOMM, MCODE), description of each. Enrichment analyses - disease (DO, DisGeNet, Network of Cancer Genes), functions (GO, KEGG, REACTOME). Survival analysis.
    - Zhang, Junpeng, Lin Liu, Taosheng Xu, Yong Xie, Chunwen Zhao, Jiuyong Li, and Thuc Duy Le. “[MiRsponge: An R/Bioconductor Package for the Identification and Analysis of MiRNA Sponge Interaction Networks and Modules](https://doi.org/10.1101/507749).” BioRxiv, January 1, 2018

- [MirGeneDB](http://mirgenedb.org/) - standardized microRNA database, 1288 microRNA families across 45 species. Downloadable, FASTA, GFF, BED files. Nomenclature refs 19, 20. 
    - Fromm, Bastian, Diana Domanska, Eirik Hoye, Vladimir Ovchinnikov, Wenjing Kang, Ernesto Aparicio-Puerta, Morten Johansen, et al. “[MirGeneDB 2.0: The Metazoan MicroRNA Complement](https://doi.org/10.1101/258749).” BioRxiv, August 13, 2019. 

- [miRPathDB](https://mpd.bioinf.uni-sb.de/about.html) - miRNA-pathway association database, human, mouse
    - Backes, Christina, Tim Kehl, Daniel Stöckel, Tobias Fehlmann, Lara Schneider, Eckart Meese, Hans-Peter Lenhof, and Andreas Keller. “[MiRPathDB: A New Dictionary on MicroRNAs and Target Pathways](https://doi.org/10.1093/nar/gkw926).” Nucleic Acids Research 45, no. D1 (January 4, 2017)

### lncRNAs

- [ncFANs v2.0](http://ncfans.gene.ac/) - functional annotation of non-coding RNAs. Three modules: ncFANs-NET, ncFANs-eLnc and ncFANs-CHIP, for annotations using pre-built coexpression/comethylation/lncRNA-gene regulatory networks (GTeX, TCGA), enhancer-derived lncRNAs (data from GRO-seq, de novo transcript assembly), microarray-based analysis, random forest-predicted networks. Input: lists of ncRNA Ensembl IDs or gene symbols. <details>
    <summary>Paper</summary>
    Zhang, Yuwei, Dechao Bu, Peipei Huo, Zhihao Wang, Hao Rong, Yanguo Li, Jingjia Liu et al. "ncFANs v2. 0: an integrative platform for functional annotation of non-coding RNAs." Nucleic Acids Research 49, no. W1 (29 May 2021): W459-W468. https://doi.org/10.1093/nar/gkab435
</details>

- [LncSEA](http://bio.liclab.net/LncSEA/index.php) - long non-coding RNA database and enrichment analysis. Covers over 50K lncRNAs, contains reference sets in 18 categories (Accessible chromatin, enhancer, super enhancer, transcription factor, survival, Drug, Disease, Cancer hallmark, subsellular location etc., Supplementary Table 2 and 3 - data sources) and 66 subcategories (based on specific attributes, overlap/proximal/closest, cancer subtypes, etc.), Table 1. Hypergeometric enrichment, Jaccard, Simpson overlaps, Correction for multiple testing (BH, Bonferroni). ID conversion. Supplementary Material 2 - details of categories. Previous tools: Co-LncRNA, Lnc-GFP, FARNA, [LnCompare](http://www.rnanut.net/lncompare/) (Supplementary Table 1 - Comparison of LncSEA with other databases and tools). [Supplementary material](https://academic.oup.com/nar/article/49/D1/D969/5921297#supplementary-data). [GitHub](https://github.com/lxy-boy/LncSEA-Code). <details>
    <summary>Paper</summary>
    Chen, Jiaxin, Jian Zhang, Yu Gao, Yanyu Li, Chenchen Feng, Chao Song, Ziyu Ning, et al. “LncSEA: A Platform for Long Non-Coding RNA Related Sets and Enrichment Analysis,” Nucleic Acids Research, 8 January 2021. https://doi.org/10.1093/nar/gkaa806
</details>

- [lncRNAKB](http://psychiatry.som.jhmi.edu/lncrnakb/) - database of long noncoding RNAs. lncRNAs are typically less conserved, expressed low on average and highly tissue-specific. Combines six resources (CHESS, LNCipedia, NONCODE, FANTOM, MiTranscriptome, BIGTranscriptome). Information about tissue-specific expression, eQTL, WGCNA co-expression to predict functions in a tissue-specific manner, random forest prediction of protein-coding score. Data: GTF gene annotation, tissue-specific expression (TPM, counts, eQTL). [RNA-seq blog post](https://www.rna-seqblog.com/lncrnakb-a-comprehensive-knowledgebase-of-long-non-coding-rnas/)
    - Seifuddin, Fayaz, Komudi Singh, Abhilash Suresh, Yun-Ching Chen, Vijender Chaitankar, Ilker Tunc, Xiangbo Ruan, et al. “[LncRNAKB: A Comprehensive Knowledgebase of Long Non-Coding RNAs](https://doi.org/10.1101/669994).” Preprint. Bioinformatics, June 13, 2019. 

- [UClncR](http://bioinformaticstools.mayo.edu/research/uclncr-pipeline/) - detecting and quantifying expression of unknown and known lncRNAs. Works for unstranded and stranded RNA-seq. Incorporates StringTie, Sebnif for novel lncRNA detection, iSeeRNA for assessing noncoding potential. Annotates lncRNAs by the nearby protein-coding genes. Tested on real data using Gencode annotations with parts of lncRNA annotations removed. 
    - Sun, Zhifu, Asha Nair, Xianfeng Chen, Naresh Prodduturi, Junwen Wang, and Jean-Pierre Kocher. “[UClncR: Ultrafast and Comprehensive Long Non-Coding RNA Detection from RNA-Seq](https://doi.org/10.1038/s41598-017-14595-3).” Scientific Reports 7, no. 1 (December 2017)


### circRNAs

- [DCC](https://github.com/dieterich-lab/DCC) Python scripts and [CircTest](https://github.com/dieterich-lab/CircTest) R visualization package - circular RNA detection. DCC uses STAR output (chimeric.out.junction) and detects back-splice junctions, filters, integrated replicate data. A much higher precision than competitors (CIRI, KNIFE), similar sensitivity. Tests for host gene-independence of circRNA expression across different experimental conditions. <details>
    <summary>Paper</summary>
    Cheng, Jun, Franziska Metge, and Christoph Dieterich. "Specific identification and quantification of circular RNAs from sequencing data." Bioinformatics 32, no. 7 (2016): 1094-1096. https://doi.org/10.1093/bioinformatics/btv656
</details>

- [CIRCpedia](http://www.picb.ac.cn/rnomics/circpedia/) database of cornRNAs from human, mouse, and some model organisms. Ribo-, poly(A)-, RNAse R methods for enriching for circRNAs. [CIRCexplorer2](https://circexplorer2.readthedocs.io/en/latest/) for the analysis of such experiments
    - Zhang et al., “Diverse Alternative Back-Splicing and Alternative Splicing Landscape of Circular RNAs.”

### Gene fusion

- [MINTIE](https://github.com/Oshlack/MINTIE) - identifying novel, rare transcriptional variants in cancer RNA-seq data. Detects fusions, transcribed structural variants (>=7bp), novel splice variants (flanked by >=20bp), complex variants (Figure 2). Filters, annotates, and prioritizes variants. Case(s) vs. control(s) analysis (single case vs. N controls). Four steps: transcriptome assembly of the case sample (SOAPdenovo-Trans), pseudo-alignment of cases and controls to an index composed of the assembled and reference transcripts (CHESS, Salmon), differential expression to identify upregulated novel features, and annotation of novel transcripts. Outperforms eight other variang detection methods on [simulated](https://github.com/Oshlack/MINTIE/tree/master/simu/run_simu.py) and experimental datasets. <details>
    <summary>
    Cmero, Marek, Breon Schmidt, Ian J. Majewski, Paul G. Ekert, Alicia Oshlack, and Nadia M. Davidson. “MINTIE: Identifying Novel Structural and Splice Variants in Transcriptomes Using RNA-Seq Data.” Genome Biology 22, no. 1 (December 2021): 296. https://doi.org/10.1186/s13059-021-02507-8.
</summary>

- [MetaFusion](https://github.com/ccmbioinfo/MetaFusion) - gene fusion caller by filtering and aggregating calls from multiple (7 by default) fusion callers (included in Docker/Singularity images, orchestrated by GenPipes). Results are summarized into new Common Fusion Format. Includes FusionAnnotator tool. [Documentation](https://github.com/ccmbioinfo/MetaFusion/wiki). <details>
    <summary>Paper</summary>
    Apostolides, Michael, Yue Jiang, Mia Husić, Robert Siddaway, Cynthia Hawkins, Andrei L Turinsky, Michael Brudno, and Arun K Ramani. “MetaFusion: A High-Confidence Metacaller for Filtering and Prioritizing RNA-Seq Gene Fusion Candidates.” Edited by Janet Kelso. Bioinformatics 37, no. 19 (October 11, 2021): 3144–51. https://doi.org/10.1093/bioinformatics/btab249.
</details>

- [CICERO](https://platform.stjude.cloud/workflows/rapid_rna-seq) - gene fusion detection, uses longer (>75bp) reads, a local assembly-based. Prioritizes candidates. Outperforms ChimeraScan, deFuse, FusionCatcher, Arriba on TCGA brain tumor data. [FusionEditor](https://proteinpaint.stjude.org/FusionEditor/) imports CICERO's output for visualization. Imports paired-end FASTQs or aligned BAMs. Supports hg19 only. [Web](https://platform.stjude.cloud/workflows/rapid_rna-seq), [GitHub](https://github.com/stjude/Cicero)
    - Tian, Liqing, Yongjin Li, Michael N. Edmonson, Xin Zhou, Scott Newman, Clay McLeod, Andrew Thrasher, et al. “[CICERO: A Versatile Method for Detecting Complex and Diverse Driver Fusions Using Cancer RNA Sequencing Data](https://doi.org/10.1186/s13059-020-02043-x).” Genome Biology 21, no. 1 (December 2020)

- [annoFuse](https://github.com/d3b-center/annoFuse/) - an R package for standartization, filtering and annotation of fusion calls detected by STAR-Fusion and Arriba, two best methods for fusion detection. Visualization options. Applied to OpenPBTA data. 
    - Gaonkar, Krutika S., Komal S. Rathi, Payal Jain, Yuankun Zhu, Miguel A. Brown, Bo Zhang, Pichai Raman, et al. “[AnnoFuse: An R Package to Annotate and Prioritize Putative Oncogenic RNA Fusions](https://doi.org/10.1101/839738).” Preprint. Bioinformatics, November 12, 2019. 

- [ChimerDB](http://203.255.191.229:8080/chimerdbv31/mindex.cdb) is a comprehensive database of fusion genes encompassing analysis of deep sequencing data and manual curations. In this update, the database coverage was enhanced considerably by adding two new modules of TCGA RNA-Seq analysis and PubMed abstract mining

- [TUMOR FUSION GENE DATA PORTAL](https://www.tumorfusions.org/) - Landscape of cancer-associated fusions using the Pipeline for RNA sequencing Data Analysis. 

- [FusionScan](http://fusionscan.ewha.ac.kr/) – prediction of fusion genes from RNA-Seq data. [RNA-seq blog post](https://www.rna-seqblog.com/fusionscan-prediction-of-fusion-genes-from-rna-seq-data/), [GitHub](https://github.com/iamlife/FusionScan)

- [Arriba](https://github.com/suhrig/arriba) - Fast and accurate gene fusion detection from RNA-Seq data 

- [FuSeq](https://github.com/nghiavtr/FuSeq) - fast fusion detection. Compared with FusionMap, TRUP, TopHat-Fusion, JAFFA, SOAPfuse
    - Vu, Trung Nghia, Wenjiang Deng, Quang Thinh Trac, Stefano Calza, Woochang Hwang, and Yudi Pawitan. “[A Fast Detection of Fusion Genes from Paired-End RNA-Seq Data](https://doi.org/10.1186/s12864-018-5156-1).” BMC Genomics 19, no. 1 (December 2018)

- [GeneFuse](https://github.com/OpenGene/GeneFuse) - Gene fusion detection and visualization

- [EricScript](https://sites.google.com/site/bioericscript/) is a computational framework for the discovery of gene fusions in paired end RNA-seq data

### Isoforms

- [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html) - An R package to Identify, Annotate and Visualize Alternative Splicing and Isoform Switches with Functional Consequences (from RNA-seq data)
    - Vitting-Seerup, Kristoffer, and Albin Sandelin. “[IsoformSwitchAnalyzeR: Analysis of Changes in Genome-Wide Patterns of Alternative Splicing and Its Functional Consequences](https://doi.org/10.1093/bioinformatics/btz247).” Bioinformatics, April 15, 2019
    - Vitting-Seerup, Kristoffer, and Albin Sandelin. “[The Landscape of Isoform Switches in Human Cancers](https://doi.org/10.1158/1541-7786.MCR-16-0459).” Molecular Cancer Research 15, no. 9 (September 2017) - Isoform switching analysis of TCGA data, tumor vs. normal. Consequences, survival prediction. [Supplementary data](https://mcr.aacrjournals.org/content/15/9/1206.figures-only) has isoform switching analysis results for all TCGA cancers

## CNVs and Structural variations

- [SuperFreq](https://github.com/ChristofferFlensburg/SuperFreq) - CNV analysis from exome data adapted for RNA-seq data. Based on log fold-change variance estimation with the neighbour correction. R package, input - BAM files (reference normal needed), variant calls from samtools or other tools, output - visualization of CNAs, other variant-related plots
    - Flensburg, Christoffer, Alicia Oshlack, and Ian J Majewski. “[Detecting Copy Number Alterations in RNA-Seq Using SuperFreq](https://www.biorxiv.org/content/10.1101/2020.05.31.126888v1)” June 01, 2020

- [CaSpER](https://github.com/akdess/CaSpER) - identification of CNVs from  RNA-seq data, bulk and single-cell (full-transcript only, like SMART-seq). Utilized multi-scale smoothed global gene expression profile and B-allele frequency (BAF) signal profile, detects concordant shifts in signal using a 5-state HMM (homozygous deletion, heterozygous deletion, neutral, one-copy-amplification, high-copy-amplification). Reconstructs subclonal CNV architecture for scRNA-seq data. Tested on GBM scRNA-seq, TCGA, other. Compared with HoneyBADGER. R code and tutorials
    - Serin Harmanci, Akdes, Arif O. Harmanci, and Xiaobo Zhou. “[CaSpER Identifies and Visualizes CNV Events by Integrative Analysis of Single-Cell or Bulk RNA-Sequencing Data.](https://doi.org/10.1038/s41467-019-13779-x)” Nature Communications 11, no. 1 (December 2020)

- [CNAPE](https://github.com/WangLabHKUST/CNAPE) - CNV detection from RNA-seq data. Regularized logistic regression (Lasso), trained on TCGA samples. Prediction accuracy >80%. R implementation
    - Mu, Quanhua, and Jiguang Wang. “[CNAPE: A Machine Learning Method for Copy Number Alteration Prediction from Gene Expression](https://doi.org/10.1109/TCBB.2019.2944827).” IEEE/ACM Transactions on Computational Biology and Bioinformatics, 2019, 1–1. 

- [CNVkit-RNA](https://github.com/etal/cnvkit) - CNV estimation from RNA-seq data. Improved moving average approach, corrects for GC content, gene expression level, gene length, correlation of gene expression and CNV (estimated from TCGA). [Docs](https://cnvkit.readthedocs.io/en/stable/), [Video tutorial](https://youtu.be/JGOIXYoLG6w)
    - Talevich, Eric, and A. Hunter Shain. “[CNVkit-RNA: Copy Number Inference from RNA-Sequencing Data.](https://doi.org/10.1101/408534)” Preprint. Bioinformatics, September 4, 2018

- [InferCNV](https://www.bioconductor.org/packages/release/bioc/html/infercnv.html) - Inferring copy number alterations from tumor single cell RNA-Seq data. R package. [GitHub wiki](https://github.com/broadinstitute/inferCNV/wiki). Part of [Trinity Cancer Transcriptome Analysis Toolkit](https://github.com/NCIP/Trinity_CTAT/wiki)

- [SQUID](https://github.com/Kingsford-Group/squid) - transcriptomic structural variation caller. Genome segment graph, then rearrange segments so that as many read alignments as possible are concordant with the rearranged sequence. Compared with MUMmer3, DELLY2, LUMPY in simulated settings, and with SOAPfuse, deFuse, FusionCatcher, JAFFA, INTEGRATE tools using real data
    - Ma, Cong, Mingfu Shao, and Carl Kingsford. “[SQUID: Transcriptomic Structural Variation Detection from RNA-Seq](https://doi.org/10.1186/s13059-018-1421-5).” Genome Biology 19, no. 1 (12 2018)

- [transindel](https://github.com/cauyrd/transIndel) - Indel caller for DNA-seq or RNA-seq


## Networks

- [ANANSE](https://github.com/vanheeringen-lab/ANANSE) (ANalysis Algorithm for Networks Specified by Enhancers) - gene regulatory network inference using TF binding profiles. Missing TF binding profiles predicted from cis-regulatory enhancer activity (H3K27ac, ATAC-seq, EP300), TF motif scores, average ChIP-seq signal of REMAP peaks in enhancers (logistic regression). Influence score - how well the expression differences between two cell types can be explained by a TF. Python implementation. [Jupyter notebooks](https://github.com/vanheeringen-lab/ANANSE-manuscript). ANANSE-inferred [Tissue-specific networks](https://doi.org/10.5281/zenodo.4814016), [cell type-specific networks](https://doi.org/10.5281/zenodo.4809062), GRNBoost2-inferred [tissue-specific networks](https://doi.org/10.5281/zenodo.4814015). [Tweet](https://twitter.com/svheeringen/status/1324251851223666689?s=20&t=uOf-0h6htdUlkWe1WVdj3A). <details>
    <summary>Paper</summary>
    Xu, Quan, Georgios Georgiou, Siebren Frölich, Maarten van der Sande, Gert Jan C Veenstra, Huiqing Zhou, and Simon J van Heeringen. “ANANSE: An Enhancer Network-Based Computational Approach for Predicting Key Transcription Factors in Cell Fate Determination.” Nucleic Acids Research 49, no. 14 (August 20, 2021): 7966–85. https://doi.org/10.1093/nar/gkab598.
</details>

- [MODifieR](https://gitlab.com/Gustafsson-lab/MODifieR) - R package wrapping 9 gene module inference methods from transcriptomics networks (WGCNA, DIAMOnD, DiffCoEx, MCODE, MODA, Module Discoverer, Clique-Sum, Correlation-Clique). Some methods include differential expression analysis. Consensus module detection. [Docker](https://hub.docker.com/r/ddeweerd/modifier). [Vignette](https://gustafsson-lab.gitlab.io/MODifieR/). <details>
    <summary>Paper</summary>
    Weerd, Hendrik A de, Tejaswi V S Badam, David Martínez-Enguita, Julia Åkesson, Daniel Muthas, Mika Gustafsson, and Zelmina Lubovac-Pilav. “MODifieR: An Ensemble R Package for Inference of Disease Modules from Transcriptomics Networks.” Edited by Lenore Cowen. Bioinformatics 36, no. 12 (June 1, 2020): 3918–19. https://doi.org/10.1093/bioinformatics/btaa235.
</details>

- [corto](https://CRAN.R-project.org/package=corto) - R package for correlation-based gene network and master regulator analysis. Can correct for CNVs. Uses RNA-seq or ATAC-seq data. Benchmarked against ARACNE-AP, minet, RTN. [GitHub](https://github.com/federicogiorgi/corto). <details>
    <summary>Paper</summary>
    Mercatelli, Daniele, Gonzalo Lopez-Garcia, and Federico M Giorgi. “[Corto: A Lightweight R Package for Gene Network Inference and Master Regulator Analysis](https://doi.org/10.1093/bioinformatics/btaa223),” Bioinformatics, Volume 36, Issue 12, 15 June 2020
</details>

- [GENIE3](https://github.com/aertslab/GENIE3) - random forest regression detection of gene modules. Input - expression matrix, output - gene x gene square co-regulation matrix

- [PANDA networks](https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/gtex-networks) - Tissue-Specific Gene Regulatory Networks constructed using [PANDA](https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/panda)

- [SCENIC networks](http://www.grndb.com/) - Tissue-specific networks inferred from single cell data using SCENIC



### Transcription regulators

- [BARTcancer](https://zanglab.github.io/bartcancer/) - predicting the transcriptional regulators from the differentially expressed genes in cancer vs. normal comparison. Database for analyzing TCGA cancers, results separately for up- and downregulated genes. Cross-referenced with [Cistrome Cancer](http://cistrome.org/CistromeCancer/)
    - Thomas, Zachary V., Zhenjia Wang, and Chongzhi Zang. “[BART Cancer: A Web Resource for Transcriptional Regulators in Cancer Genomes](https://doi.org/10.1101/2020.12.17.423327).” Preprint. Cancer Biology, December 18, 2020. 

- [BARTweb](http://bartweb.org/) - identifying TFs whose genomic binding patterns associate with input genomic features. Input - a gene set, a ChIP-seq mapped read dataset, as scored genomic region set. Adaptive lasso to prioritize TFs. hg38 and mm10 supported. Evaluated on [KnockTF data](http://www.licpathway.net/KnockTF/). Similar tools: [TFEA.ChIP](https://bioconductor.org/packages/TFEA.ChIP/), [ChEA3](https://maayanlab.cloud/chea3/)
    - Ma, Wenjing, Zhenjia Wang, Yifan Zhang, Neal E Magee, Yang Chen, and Chongzhi Zang. “[BARTweb: A Web Server for Transcription Factor Association Analysis](https://doi.org/10.1101/2020.02.17.952838).” Preprint. Bioinformatics, February 17, 2020


## Integrative

- Review of tools and methods for the integrative analysis of multiple omics data, cancer-oriented. [Table 1](https://journals.sagepub.com/na101/home/literatum/publisher/sage/journals/content/bbia/2020/bbia_14/1177932219899051/20200130/images/large/10.1177_1177932219899051-table1.jpeg) - multi-omics data repositories (TCGA, CPTAC, ICGC, CCLE, METABRIC, TARGET, Omics Discovery Index). Three broad areas of multi-omics analysis: 1. Disease subtyping and classification based on multi-omics profiles; 2. Prediction of biomarkers for various applications including diagnostics and driver genes for diseases; 3. Deriving insights into disease biology. [Table 2](https://journals.sagepub.com/na101/home/literatum/publisher/sage/journals/content/bbia/2020/bbia_14/1177932219899051/20200130/images/large/10.1177_1177932219899051-table2.jpeg) - software categorized by use case (PARADIGM, iClusterPlus, PSDF, BCC, MDI, SNF, PFA, PINSPlus, NEMO, mixOmics, moCluster, MCIA, JIVE, MFA, sMBPLS, T-SVD, Joint NMF). Brief description of each tool, links, exemplary publications. [Table 3](https://journals.sagepub.com/na101/home/literatum/publisher/sage/journals/content/bbia/2020/bbia_14/1177932219899051/20200130/images/large/10.1177_1177932219899051-table3.jpeg) - visualization portals (cBioPortal, Firebrowse, UCSC Xena, LinkedOmics, 3Omics, NetGestalt, OASIS, Paintomics, MethHC). Description of each, data types, analysis examples.
    - Subramanian, Indhupriya, Srikant Verma, Shiva Kumar, Abhay Jere, and Krishanpal Anamika. “[Multi-Omics Data Integration, Interpretation, and Its Application](https://doi.org/10.1177/1177932219899051).” Bioinformatics and Biology Insights, January 31, 2020 

- [SmCCNet](https://github.com/KechrisLab/SmCCNet) - R package for integrating one or multiple types of omics data with a quantitative or binary phenotype. Based on the concept of sparse multiple canonical analysis (SmCCA) and sparse partial least squared discriminant analysis (SPLSDA) and aims to find relationships between omics data and a specific phenotype. The framework uses LASSO (Least Absolute Shrinkage and Selection Operator) for sparsity constraints, allowing it to identify significant features within the data. Weighted and unweighted modes, to prioritize different data types. Visualization via [SmCCNet Visualization Shiny app](https://smccnet.shinyapps.io/smccnetnetwork/), Cytoscape interface. <details>
    <summary>Paper</summary>
    Liu, Weixuan, Thao Vu, Iain Konigsberg, Katherine Pratte, Yonghua Zhuang, and Katerina Kechris. “SmCCNet 2.0: A Comprehensive Tool for Multi-Omics Network Inference with Shiny Visualization,” https://doi.org/10.1101%2F2023.11.20.567893
</details>

- [MOGONET](https://github.com/txWang/MOGONET) (Multi-Omics Graph cOnvolutional NETworks) - multi-omics integration method using graph convolutional network for classification and biomarker detection. Utilizes View Correlation Discovery Network (VCDN) to explore the cross-omics correlations at the label space. Weighted similarity networks using cosine similarity. Tested various configurations, outperforms baseline methods (KNN, SVM, etc.). Applied to classifying Alzhemier's disease (ROSMAP dataset), LGG, KIRPAN, BRCA normal/cancer and PAM50 subtypes. [Supplementary tables](https://www.nature.com/articles/s41467-021-23774-w#Sec21) S9-11 rank biomarkers. Python implementation. <details>
    <summary>Paper</summary>
    Wang, Tongxin. “MOGONET Integrates Multi-Omics Data Using Graph Convolutional Networks Allowing Patient Classification and Biomarker Identification,”  https://doi.org/10.1038/s41467-021-23774-w 
</details>


- [AJIVE](https://github.com/idc9/py_jive) (Angle-based Joint and Individual Variation Explained) - joint dimension reduction to horizontally (across cohorts) integrate gene expression data across PDXs and human tumor cohorts. Identifies the shared and individual variation of genes across model systems and human cohorts. An extension of [JIVE](https://cran.r-project.org/web/packages/r.jive/vignettes/BRCA_Example.html), employs thresholded SVD, then Principal Angle Analysis to create joint and individual representations of the original inputs. Applied to integrating TCGA and CCLE, TCGA and genetically engineered mouse models (GEMMs), other examples of using joint variation to identify lapatinib response factors. Python, with [R implementation](https://github.com/idc9/r_jive). <details>
    <summary>Paper</summary>
    Price, Brandon A., J. S. Marron, Lisle E. Mose, Charles M. Perou, and Joel S. Parker. “Translating Transcriptomic Findings from Cancer Model Systems to Humans through Joint Dimension Reduction.” Communications Biology 6, no. 1 (February 16, 2023): 179. https://doi.org/10.1038/s42003-023-04529-3.
</details>

- [DIABLO](http://mixomics.org/) - multi-omics analysis method. Overview of previous methods (SNF, Bayesian Consensus Clustering, NMF, JIVE, sGCCA, MOFA, others). Method extends sGCCA multivariate dimensionality reduction that uses SVD and selects co-expressed (correlated) variables from several omics datasets. Methods, model, iterative solution. Design matrix specifies which omics datasets are connected. Variable selection for biomarkers identification.  Visualization options. Part of [mixOmics R package](http://mixomics.org/), [Documentation](https://mixomicsteam.github.io/Bookdown/intro.html)
    - Singh, Amrit, Casey P Shannon, Benoît Gautier, Florian Rohart, Michaël Vacher, Scott J Tebbutt, and Kim-Anh Lê Cao. “DIABLO: An Integrative Approach for Identifying Key Molecular Drivers from Multi-Omics Assays.” Edited by Inanc Birol. Bioinformatics 35, no. 17 (September 1, 2019): 3055–62. https://doi.org/10.1093/bioinformatics/bty1054.

- [MANCIE](https://cran.r-project.org/web/packages/MANCIE/) - matrix analysis and normalization by concordant information enhancement. Bias correction and data integration of distinct genomic profiles on the same samples. Match matrices by rows, run correlation for each row, replace the associated row with modified values using a PCA procedure, Methods. Tested on integration of DHS and gene expression data, TCGA and METABRIC data. R package
    - Zang, Chongzhi, Tao Wang, Ke Deng, Bo Li, Sheng’en Hu, Qian Qin, Tengfei Xiao, et al. “[High-Dimensional Genomic Data Bias Correction and Data Integration Using MANCIE](https://doi.org/10.1038/ncomms11305).” Nature Communications 7, no. 1 (December 2016). 

- [JIVE](https://genome.unc.edu/jive/) - Joint and Individual Variation Explained. Decomposition of (X) multiple (i) omics datasets into three terms: low-rank (constrained) matrices capturing joint variation (J), plus structured variation (A_i) and residual noise. Data are row-centered and scaled by its total variation. Main constrain: the rows of joint and individual matrices should be orthogonal. Estimate matrices by iteratively minimizing ||R||^2 (R=X-J-A). Relationship to PCA, CCA, PLS. Illustrated on TCGA GBM gene expression, methylation, and miRNA data, with interpretation. [Matlab code](https://genome.unc.edu/jive/), [r.jive package](https://cran.r-project.org/web/packages/r.jive/vignettes/BRCA_Example.html)
    - Lock, Eric F., Katherine A. Hoadley, J. S. Marron, and Andrew B. Nobel. “[JOINT AND INDIVIDUAL VARIATION EXPLAINED (JIVE) FOR INTEGRATED ANALYSIS OF MULTIPLE DATA TYPES](https://doi.org/10.1214/12-AOAS597).” The Annals of Applied Statistics 7, no. 1 (March 1, 2013)

- [List of software packages for multi-omics analysis](https://github.com/mikelove/awesome-multi-omics), by Mike Love. Slides for the talk "[Assessing consistency of unsupervised multi-omics methods](https://docs.google.com/presentation/d/1QAaweEc32JzhWHl7YenLdT9w8JUjwaTExe_uve2s22U/edit#slide=id.p)". 


## Classification

- [MLSeq](https://www.bioconductor.org/packages/release/bioc/html/MLSeq.html) - Machine learning interface for RNA-Seq data
    - Zararsız, Gökmen, Dincer Goksuluk, Selcuk Korkmaz, Vahap Eldem, Gozde Erturk Zararsiz, Izzet Parug Duru, and Ahmet Ozturk. “[A Comprehensive Simulation Study on Classification of RNA-Seq Data](https://doi.org/10.1371/journal.pone.0182507).” PLOS ONE 12, no. 8 (August 23, 2017)


## Visualization

- [aPEAR](https://github.com/ievaKer/aPEAR) (Advanced Pathway Enrichment Analysis Representation) - R package for pathway enrichment results visualization as a network. Input: results from clusterProfiler or gprofiler2. Pairwise similarity between pathways is estimated using Jaccard (or, cosine, correlation), the similarity matrix is clustered using Markov (or, hierarchical, spectral) clustering, the most representative pathway name is determined using PageRank (or, HITS, highest NES) algorithm. Output: ggplot2 object, modifiable, can be plotted interactive with plotly. [CRAN](https://CRAN.R-project.org/package=aPEAR). <details>
    <summary>Paper</summary>
    Kerseviciute, Ieva, and Juozas Gordevicius. “APEAR: An R Package for Autonomous Visualisation of Pathway Enrichment Networks.” Preprint. Bioinformatics, March 29, 2023. https://doi.org/10.1101/2023.03.28.534514.
</details>

- [chromoMap](https://CRAN.R-project.org/package=chromoMap) - an R package/function for visualizing BED-like data across chromosomes. Static and interactive (Shiny embeddable) plots, segment-, point-, barplot-, scatterplot visualization. Filters to color visualization by criteria. Colors, spacing, height/width - all customizable. Supports multiple organisms. [Documentation, tutorial](https://lakshay-anand.github.io/chromoMap/index.html)
    - Anand, Lakshay, and Carlos M Rodriguez Lopez. “[ChromoMap: An R Package for Interactive Visualization and Annotation of Chromosomes](https://doi.org/10.1101/605600)” bioRxiv, January 23, 2020

- [Hiplot](https://hiplot-academic.com/) - a web tool for publication-ready biomedical data visualization (command line interpreter Hctl written in Go is available). From heatmaps, correlograms to dimensionaliry reduction, other plots (over 240). Basic statistics, multi-omics, regression, clustering, dimensionality reduction, meta-analysis, survival analysis, risk modelling, etc., highly adjustable. Wraps other packages, such as immundeconv, circos. JSON-based plugin system, [hiplotlib](https://github.com/hiplot/hiplotlib) - development library, [plugin-preview](https://hiplot-academic.com/developer/plugin-preview) - JSON-based previewer to debug plugins. Demo data, data spreadsheet/upload/remote access, SVG download. Created in collaboration with [OpenBiox](https://openbiox.org/#/projects) consortium. Table 1 - comparison with competitors, [ImageGP](http://www.ehbio.com/ImageGP/) and [Galaxy](https://usegalaxy.org/). <details>
    <summary>Paper</summary>
    Li, Jianfeng, Benben Miao, Shixiang Wang, Wei Dong, Houshi Xu, Chenchen Si, Wei Wang, et al. “Hiplot: A Comprehensive and Easy-to-Use Web Service for Boosting Publication-Ready Biomedical Data Visualization.” Briefings in Bioinformatics 23, no. 4 (July 18, 2022): bbac261. https://doi.org/10.1093/bib/bbac261.
</details>

- [genomation](https://bioconductor.org/packages/devel/bioc/vignettes/genomation/inst/doc/GenomationManual.html) - a toolkit for annotation and visualization of genomic data, R package

- [karyoploteR](https://github.com/bernatgel/karyoploteR) - An R/Bioconductor package to plot arbitrary data along the genome

- [pcaExplorer](https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html) - Interactive Visualization of RNA-seq Data Using a Principal Components Approach, R package

- [WIlsON](https://github.com/loosolab/wilson) - Web-based Interactive Omics VisualizatioN, accepts, text files, SummarizedExperiment datasets. R/Shiny, installs as a package. Docker image available. [Web demo](http://loosolab.mpi-bn.mpg.de/wilson/). <details>
    <summary>Paper</summary>
    H. Schultheis, C. Kuenne, J. Preussner, R. Wiegandt, A. Fust, M. Bentsen and M. Looso. WIlsON: Webbased Interactive Omics VisualizatioN. Bioinformatics 35(6) 2018, doi: https://doi.org/10.1093/bioinformatics/bty711
</details>

## Data

- [ARCHS4](http://amp.pharm.mssm.edu/archs4/download.html) - Massive Mining of Publicly Available RNA-seq Data from Human and Mouse

- [cBioPortalData](https://github.com/waldronlab/cBioPortalData) - cBioPortal data as MultiAssayExperiment objects, by Waldron Lab

- [curatedTCGAData](https://github.com/waldronlab/curatedTCGAData) - Curated Data From The Cancer Genome Atlas (TCGA) as MultiAssayExperiment objects, by Waldron Lab

- [gtexRNA](https://github.com/sigven/gtexRNA) - R package for retrieval of tissue-specific expression data from GTEx. By Sigve Nakken, [website](https://sigven.github.io/gtexRNA/)

- [GTEx Visualizations](https://github.com/broadinstitute/gtex-viz) - web-based visualization tools for exploring tissue-specific gene expression and regulation

- [PINS](http://www.cs.wayne.edu/tinnguyen/PINS/PINS.html) - A novel method for data integration and disease subtyping

- [refine.bio](https://www.ccdatalab.org/projects/refinebio) - harmonized microarray and RNA-seq data for various organisms and conditions

- [recount2](https://bioconductor.org/help/workflows/recountWorkflow/) - an R workflow to work with recount2 data

- [GREIN](http://www.ilincs.org/apps/grein/) - re-analysis of RNA-seq datasets from GEO. Download processed data, visualization, power analysis, differential expression, functional enrichment analysis, connectivity analysis with LINCS L1000 data. [GitHub](https://github.com/uc-bd2k/grein), [Docker image](https://hub.docker.com/r/ucbd2k/grein/)
    - Al Mahi, Naim, Mehdi Fazel Najafabadi, Marcin Pilarczyk, Michal Kouril, and Mario Medvedovic. “[GREIN: An Interactive Web Platform for Re-Analyzing GEO RNA-Seq Data](https://doi.org/10.1101/326223),” October 27, 2018

- [DEE2](http://dee2.io) - Digital Expression Explorer - gene- and transcript-level processed data from multiple organisms, amenable for downstream analysis in R etc. [getDEE2](https://github.com/markziemann/dee2/blob/master/AccessDEEfromR.md) R package to get the data
    - Ziemann, Mark, Antony Kaspi, and Assam El-Osta. “[Digital Expression Explorer 2: A Repository of Uniformly Processed RNA Sequencing Data](https://doi.org/10.1093/gigascience/giz022).” GigaScience 8, no. 4 (April 1, 2019)

- [GEMMA](https://gemma.msl.ubc.ca/home.html) - curated transcriptomic database, >10,000 studies, \~34% are brain-related. Query genes, phenotypes, experiments, search for coexpression, differential expression. Processing methods, batch correction. Online access, API, R package. [GitHub](https://github.com/PavlidisLab/Gemma/)
    - Lim, Nathaniel, Stepan Tesar, Manuel Belmadani, Guillaume Poirier-Morency, Burak Ogan Mancarci, Jordan Sicherman, Matthew Jacobson, Justin Leong, Patrick Tan, and Paul Pavlidis. “[Curation of over 10,000 Transcriptomic Studies to Enable Data Reuse](https://doi.org/10.1101/2020.07.13.201442).” Preprint. Bioinformatics, July 14, 2020. 

### Genes

- Impact of hg19, hg38, CHM13 genome build on RNA-seq gene expression. GENCODEv35 gene annotations, lifted over. STAR alignment, RSEM quantification, regressing out batch, RIN, sex, splicing (LeafCutter). Effect on clinically relevant genes. [Supplementary tables](https://www.cell.com/ajhg/fulltext/S0002-9297(24)00168-X#supplementaryMaterial) S2 - Annotation-specific genes, S3 - differentially quantified genes, S4 - build-exclusive genes ("exclude" column flagging genes). Genome-based alignment is recommended. [GitHub](https://github.com/raungar/build_rnaseq_paper_public). <details>
    <summary>Paper</summary>
    Ungar, Rachel A., Pagé C. Goddard, Tanner D. Jensen, Fabien Degalez, Kevin S. Smith, Christopher A. Jin, Devon E. Bonner, Jonathan A. Bernstein, Matthew T. Wheeler, and Stephen B. Montgomery. “Impact of Genome Build on RNA-Seq Interpretation and Diagnostics.” The American Journal of Human Genetics 111, no. 7 (July 2024): 1282–1300. https://doi.org/10.1016/j.ajhg.2024.05.005.
</details>

- [Enrichr](https://maayanlab.cloud/Enrichr/) - enrichment analysis, gene search, term search. Libraries for various signatures are available for [download](https://maayanlab.cloud/Enrichr/#stats)

- [CellMarker](http://biocc.hrbmu.edu.cn/CellMarker/download.jsp) - Cell markers of different cell types from different tissues in human and mouse.

- [A list of updated 1439 DNA-binding transcription factors](https://www.ebi.ac.uk/QuickGO/targetset/dbTF) from re-annotation study of transcription factors in Gene Ontology annotations
    - Lovering, Ruth C., Pascale Gaudet, Marcio L. Acencio, Alex Ignatchenko, Arttu Jolma, Oriol Fornes, Martin Kuiper, et al. “[A GO Catalogue of Human DNA-Binding Transcription Factors](https://doi.org/10.1101/2020.10.28.359232).” BioRxiv, January 1, 2020

- [List of gene lists for genomic analyses](https://github.com/macarthur-lab/gene_lists) - GitHub repo with tab-separated annotated lists

- [CREEDS](https://amp.pharm.mssm.edu/creeds/) - database of manually (and automatically) extracted gene signatures. Single gene perturbations, disease signatures, single drug perturbations. Batch effect correction, when necessary. Overall, good agreement with MSigDb C2. Characteristic Direction (CD) method to detect differential genes. [API access in R](http://rpubs.com/wangz10/177826)
    - Wang, Zichen, Caroline D. Monteiro, Kathleen M. Jagodnik, Nicolas F. Fernandez, Gregory W. Gundersen, Andrew D. Rouillard, Sherry L. Jenkins, et al. “[Extraction and Analysis of Signatures from the Gene Expression Omnibus by the Crowd](https://doi.org/10.1038/ncomms12846).” Nature Communications 7, no. 1 (November 2016)

- [DIOPT](https://www.flyrnai.org/cgi-bin/DRSC_orthologs.pl) - finding ortholog genes among human, mouse, zebrafish, C. elegans, Drozophila, S. cerevisiae. Integration with human GWAS allows to search for orthologs for diseases and traits. Batch conversion, filtering. [DIOPT-DIST](https://www.flyrnai.org/diopt-dist) - DIOPT Diseases and Traits. <details>
    <summary>Paper</summary>
    Hu, Yanhui, Ian Flockhart, Arunachalam Vinayagam, Clemens Bergwitz, Bonnie Berger, Norbert Perrimon, and Stephanie E Mohr. “An Integrative Approach to Ortholog Prediction for Disease-Focused and Other Functional Studies.” BMC Bioinformatics 12, no. 1 (December 2011): 357. https://doi.org/10.1186/1471-2105-12-357.
</details>

- [OGEE v3](https://v3.ogee.info/#/home) - online gene essentiality database. Experimentally validated using CRISPR and RNAi. All data (essential genes, their properties, expression data and more) are downloadable 
    - Gurumayum, Sanathoi, Puzi Jiang, Xiaowen Hao, Tulio L Campos, Neil D Young, Pasi K Korhonen, Robin B Gasser, et al. “[OGEE v3: Online GEne Essentiality Database with Increased Coverage of Organisms and Human Cell Lines](https://doi.org/10.1093/nar/gkaa884)”

- [HGNChelper](https://waldronlab.io/HGNChelper/) - Identify and Correct Invalid HGNC Human Gene Symbols and MGI Mouse Gene Symbols

- [gencode_regions](https://github.com/saketkc/gencode_regions) - Extract 3'UTR, 5'UTR, CDS, Promoter, Genes, Introns etc from GTF files, by Saket Choudhary

- [Extract intron boundaries per transcript](https://gist.github.com/hiraksarkar/ce8a71a6953cb4e9823d868c283bf99d)

- [CHESS](http://ccb.jhu.edu/chess/) - Comprehensive Human Expressed SequenceS, database of novel genes, identified from GTeX data, protein-coding and lncRNA

- [ideogram](https://github.com/eweitz/ideogram) - Chromosome visualization with D3.js. [Examples](https://eweitz.github.io/ideogram/). [ideogRam R wrapper](https://github.com/freestatman/ideogRam)

- [Enhancer-promoter (EP) pairs from Thurman et al., (2012)](ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz)

- [pysradb](https://github.com/saketkc/pysradb) - Python package for interacting with SRAdb and downloading datasets from SRA, by Saket Choudhary. [Documentation](https://www.saket-choudhary.me/pysradb/)


## Misc

- [GEOparse](https://github.com/guma44/GEOparse) - Python library to access Gene Expression Omnibus Database (GEO). [Documentation](https://geoparse.readthedocs.io/en/latest/)

- [ffq](https://github.com/pachterlab/ffq) - A tool to find sequencing data and metadata from public databases (GEO, SRA, EMBL-EBI, others).

- [RRHO](https://systems.crump.ucla.edu/rankrank/) - Rank–rank Hypergeometric Overlap between two gene lists ranked by the degree of differential expression (e.g., signed -log10 p-value). 2D analog of GSEA.  Identifies and visualizes areas of significant overlap by determining the degree of statistical enrichment using the hypergeometric distribution while sliding across all possible thresholds through the two ranked lists. Multiple testing correction (FWER, Benjamini-Yekutieli). [Website](https://systems.crump.ucla.edu/rankrank/) and [R/Bioconductor RRHO package](https://bioconductor.org/packages/RRHO/). <details>
    <summary>Paper</summary>
    Plaisier, Seema B, Richard Taschereau, Justin A Wong, and G Graeber. “Rank–Rank Hypergeometric Overlap: Identification of Statistically Significant Overlap between Gene-Expression Signatures.” Nucleic Acids Research 38, no. 17 (2010): 17. https://doi.org/10.1093/nar/gkq636
</details>

- [gffio](https://github.com/lh3/gffio) - tool to process GFF3 and GTF files. It can convert between GFF3 and GTF, generate 12-column BED, extract CDS/transcript/protein sequences, reorder features and select the longest CDS/transcript. By [Heng Li](https://github.com/lh3)

- [gtftk](https://github.com/dputhier/pygtftk) - A python package and a set of shell commands to handle GTF files. Subcommands for editing GTF files, getting information and summary statistics, selecting by various criteria, converting BED to gtf and other formats, annotating by closest genes and more, getting sequences, coordinates of specific gene elements, coverage profile and other bigWig operations. 
    - Ferre, Q, G Charbonnier, N Sadouni, F Lopez, Y Kermezli, S Spicuglia, C Capponi, B Ghattas, and D Puthier. “[OLOGRAM: Determining Significance of Total Overlap Length between Genomic Regions Sets](https://doi.org/10.1093/bioinformatics/btz810),” Bioinformatics, 15 March 2020 - [OLOGRAM](https://github.com/dputhier/pygtftk) - overlap statistics between sets of genomic regions (BED or GTF). Uses Monte Carlo simulation to model (negative binomial) both the distributions of region length and inter-region length, fast. [Supplementary data](https://doi.org/10.1093/bioinformatics/btz810): Table 1 - comparison of features with Genomic Hyperbrowser, GREAT, CEAS, Bedtools fisher, LOLA. Note - details of the method, test results on experimental data. [Supplementary GitHub repo](https://github.com/dputhier/ologram_supp_mat) - Snakefiles to generate figures.

- [Recommended Coverage and Read Depth for NGS Applications](https://genohub.com/recommended-sequencing-coverage-by-application/), by GenoHub. And, their [NGS Handbook](https://genohub.com/next-generation-sequencing-handbook/)

- [BioJupies](https://amp.pharm.mssm.edu/biojupies/) - analysis of GEO/GTEx data or your own gene expression table/FASTQ in autogenerated Jupyter notebook. Rich set of tools for EDA (PCA, Clustergrammer, Library size analysis), Differential expression analysis (Volcano, MA plots), Enrichment analysis (Enrichr, GO, Pathway, TF, Kinase, miRNA enrichments), L1000 signatures. Best suited for two-group analysis. Includes Methods for the selected tools

- [HGNChelper](https://cran.r-project.org/web/packages/HGNChelper/index.html) - Handy Functions for Working with HGNC Gene Symbols and Affymetrix Probeset Identifiers

- `tximport` - importing transcript abundance datasets from Salmon, Sailfish, kallisto, RSEM, and differential analysis

- [rpkmforgenes](http://sandberg.cmb.ki.se/media/data/rnaseq/instructions-rpkmforgenes.html) - a Python script for calculating gene expression for RNA-Seq data

- [Python interface to access reference genome features from Ensembl](https://github.com/openvax/pyensembl), e.g., genes, transcripts, and exons

- [TPMcalculator](https://github.com/ncbi/TPMCalculator) - converts gene counts to TPM using transcript information from a GTF file. TPM vs. FPKM correlation for validation. C/C++ command line tool, Docker image, CWL workflow
    - Vera Alvarez, Roberto, Lorinc Sandor Pongor, Leonardo Mariño-Ramírez, and David Landsman. “[TPMCalculator: One-Step Software to Quantify MRNA Abundance of Genomic Features](https://doi.org/10.1093/bioinformatics/bty896).” Bioinformatics 35, no. 11 (June 1, 2019)

- Multi-omics madness picture, [Tweet](https://twitter.com/AntoBeck/status/1461478948106170374?s=20), [download](multiomics.jpeg)

![Multi-omics madness](https://pbs.twimg.com/media/FEg3MGmVQAIJnoQ?format=jpg&name=small)

- STAR aligner parameters, from https://doi.org/10.1038/s41467-020-18035-1
    - **ATAC-seq:** --alignEndsType End- ToEnd --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignIntronMax 1 --alignSJDBoverhangMin 999 --alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 2000 --outMultimapperOrder Random --outFilterMultimapNmax 999 --outSAMmultNmax 1
    - **RNA-seq:** --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.9 --out- FilterMatchNminOverLread 0.9 --outFilterMatchNmin 20 --alignIntronMax 200000 --alignMatesGapMax 2000 --alignEndsProtrude 10 ConcordantPair --outMultimapperOrder Random --outFilterMultimapNmax 999
    - **ChIP-seq:** --outFilterMismatchNoverLmax 0.2 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignIntronMax 1 --alignSJDBoverhangMin 999 --outFilterMultimapNmax 1 --alignEndsProtrude 10 ConcordantPair