# A continually expanding collection of RNA-seq tools

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

RNA-seq related tools and genomics data analysis resources. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Pipelines](#pipelines)
  - [Preprocessing](#preprocessing)
  - [Analysis](#analysis)
- [Quality control](#quality-control)
- [Batch effect](#batch-effect)
- [Clustering](#clustering)
- [Timecourse](#timecourse)
- [Differential expression](#differential-expression)
- [Functional enrichment](#functional-enrichment)
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
- [Integrative](#integrative)
- [Classification](#classification)
- [Visualization](#visualization)
- [Data](#data)
  - [Genes](#genes)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Pipelines

- [RNA-seq pipeline developed by NASA GeneLab](https://github.com/nasa/GeneLab_Data_Processing). FastQC/MultiQC, TrimGalore, STAR (two-pass mode), RSEM (to quantify isiforms), DESeq2. Can normalize using ERCCs. produces unnormalized and normalized counts, more. PCA, heatmaps, other visualization. [Step-by-step instructions and commands](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/GL-DPPD-7101-C.md)
    - Galazka, Jonathan. “[NASA GeneLab RNA-Seq Consensus Pipeline: Standardized Processing of Short-Read RNA-Seq Data](https://doi.org/10.1101/2020.11.06.371724),” bioRxiv, November 10, 2020

- [GEO2RNAseq](https://anaconda.org/xentrics/r-geo2rnaseq) - an R-based RNA-seq processing pipeline, from FASTQ files or directly from GEO. Reads meta-data, oriented for two-group differential analysis. Competitors: Galaxy, the Total RNA-seq Analysis Package for R (TRAP), EasyRNASeq, READemption. Installable through Conda https://anaconda.org/xentrics/r-geo2rnaseq
    - Seelbinder, Bastian, Thomas Wolf, Steffen Priebe, Sylvie McNamara, Silvia Gerber, Reinhard Guthke, and Joerg Linde. “[GEO2RNAseq: An Easy-to-Use R Pipeline for Complete Pre-Processing of RNA-Seq Data](https://doi.org/10.1101/771063).” Preprint. Bioinformatics, September 16, 2019

### Preprocessing

- [Check strandedness of RNA-Seq fastq files](https://github.com/betsig/how_are_we_stranded_here)

- [Illumina Instrument Type from fastq](https://www.biostars.org/p/198143/)

- [adapters](https://github.com/stephenturner/adapters) - Adapter sequences for trimming, by Stephen Turner

- [bamcov](https://github.com/fbreitwieser/bamcov) - Quickly calculate and visualize sequence coverage in alignment files in command line

- [bamcount](https://github.com/BenLangmead/bamcount) - BigWig and BAM utilities, coverage, by Ben Langmead

- [bigwig-nim](https://github.com/brentp/bigwig-nim) - single static binary + liberal license of tool to convert bed to bigwig (and back) and get fast coverage stats from bigwig, by Brent Pedersen, [Twitter](https://twitter.com/brent_p/status/1162765386548305923?s=03)

- [cgpBigWig](https://github.com/cancerit/cgpBigWig) - Package of C scripts for generation of BigWig coverage files. 

- [covviz](https://github.com/brwnj/covviz) - calculate and view coverage based variation

- [indexcov](https://github.com/brentp/goleft/tree/master/indexcov) - Quickly estimate coverage from a whole-genome bam or cram index

- [fastq-pair](https://github.com/linsalrob/fastq-pair) - Match up paired end fastq files quickly and efficiently



### Analysis

- [Omics Playground](https://github.com/bigomics/omicsplayground) - self-service bioinformatics platform to analyze, visualize, and integrate omics data. Descriptive statistics, differential gene expression, clustering, GSEA/KEGG, signature, cell type prediction. Input - text-based files. R/Shiny implementation, Docker image. [GitHub](https://github.com/bigomics/omicsplayground), [Documentation](https://omicsplayground.readthedocs.io/)
    - Akhmedov, Murodzhon, Axel Martinelli, Roger Geiger, and Ivo Kwee. “[Omics Playground: A Comprehensive Self-Service Platform for Visualization, Analytics and Exploration of Big Omics Data](https://doi.org/10.1093/nargab/lqz019).” NAR Genomics and Bioinformatics 2, no. 1 (March 1, 2020)

- [Shiny-Seq](https://github.com/schultzelab/Shiny-Seq) - Shiny app, and Docker image. Input - a count table, or kallisto output, and annotations. Normalization (DESeq2), batch effect removal (limma or SVA),  differential expression analysis (DeSeq2), co-expression network analysis (WGCNA), functional enrichment analysis (clusterprofiler), TFBS motif overrepresentation (pcaGopromoter). Visualization (heatmaps, volcano plots). [GitHub](https://github.com/schultzelab/Shiny-Seq). [RNA-seq blog post](https://www.rna-seqblog.com/shiny-seq-advanced-guided-transcriptome-analysis/)
    - Sundararajan, Zenitha, Rainer Knoll, Peter Hombach, Matthias Becker, Joachim L. Schultze, and Thomas Ulas. “[Shiny-Seq: Advanced Guided Transcriptome Analysis](https://doi.org/10.1186/s13104-019-4471-1).” BMC Research Notes 12, no. 1 (December 2019). 

- [3D_RNA-seq]( https://github.com/wyguo/ThreeDRNAseq) - R package and Shiny app, and Docker image, for differential expression, differential alternative splicing, and differential Transcript Usage. Two- and mutliple group analysis, time course. Input - Salmon/Kallisto transcript quantification files, or .csv. Diagnostic plots, PCA, batch removal using RUVseq, limma-voom for differential expression, iso-kTSP and TSIS for isoform switching between groups and in time course, respectively. Heatmaps, barplots, volcano plots, Venn diagrams. 3D stands for three days. [Web](https://ics.hutton.ac.uk/3drnaseq/),  [Manual](https://github.com/wyguo/ThreeDRNAseq/blob/master/vignettes/user_manuals/3D_RNA-seq_App_manual.md),  [GitHub](https://github.com/wyguo/ThreeDRNAseq)
    - Guo, Wenbin, Nikoleta Tzioutziou, Gordon Stephen, Iain Milne, Cristiane Calixto, Robbie Waugh, John WS Brown, and Runxuan Zhang. “[3D RNA-Seq - a Powerful and Flexible Tool for Rapid and Accurate Differential Expression and Alternative Splicing Analysis of RNA-Seq Data for Biologists](https://doi.org/10.1101/656686).” Preprint. Bioinformatics, May 31, 2019.

- [DrEdGE](https://github.com/ptgolden/dredge) - Differential Expression Gene Explorer, Takes in a table of transcript abundance counts, experimental design, other input can be generated in R. Plot, table, heatmap visualization options. Web: http://dredge.bio.unc.edu/ , GitHub: https://github.com/ptgolden/dredge
    - Tintori, Sophia C, Patrick Golden, and Bob Goldstein. “[Differential Expression Gene Explorer (DrEdGE): A Tool for Generating Interactive Online Data Visualizations for Exploration of Quantitative Transcript Abundance Datasets](https://doi.org/10.1101/618439).” Preprint. Genomics, April 25, 2019. 

- [Phantasus](https://github.com/ctlab/phantasus) - interactive exploratory analyses of genomic data, from clustering, PCA, to enrichment, network, and pathway analyses. Works on user data, ARCHS4, TCGA. https://github.com/ctlab/phantasus. [Bioconductor R package](http://www.bioconductor.org/packages/release/bioc/html/phantasus.html) and [Docker image](https://hub.docker.com/r/dzenkova/phantasus)

## Quality control

- [fastp](https://github.com/OpenGene/fastp) - fast C++ parallelized tool for FASTQ quality control, adapter trimming, quality filtering, pruning, polyX (polyG) trimming, works with single- and paired-end data. [GitHub](https://github.com/OpenGene/fastp)
    - Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. “[Fastp: An Ultra-Fast All-in-One FASTQ Preprocessor](https://doi.org/10.1093/bioinformatics/bty560).” Bioinformatics 34, no. 17 (September 1, 2018)

- [MultiQC](http://multiqc.info/) - Summarization and visualization QC results for multiple samples in one report. Recognizes multiple QC tools

- [ngsReports](https://github.com/UofABioinformaticsHub/ngsReports) - An R Package for managing FastQC reports and other NGS related log files. [biorXiv](http://biorxiv.org/content/early/2018/05/02/313148.abstract)

- [sickle](https://github.com/najoshi/sickle) - A windowed adaptive trimming tool for FASTQ files using quality. Post-adapter trimming step

- [fastqc](https://github.com/Malarkey73/fastqc) - an R package for quality control (QC) of short read fastq files, analog of the [original FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- [FastQt](https://github.com/labsquare/FastQt) - FastQC port to Qt5: A quality control tool for high throughput sequence data

- [fastqcheck](https://github.com/VertebrateResequencing/fastqcheck) - Generate statistics on and validate fastq files

## Batch effect

- [ComBat-seq](https://github.com/zhangyuqing/sva-devel) - batch effect correction for RNA-seq data using negative binomial regression. Maintains count nature of RNA-seq data. Tested on simulated data (polyester package), and experimental data. Achieves the highest true positive rate. [Code to reproduce paper](https://github.com/zhangyuqing/ComBat-seq)
    - Zhang, Yuqing, Giovanni Parmigiani, and W. Evan Johnson. “[ComBat-Seq: Batch Effect Adjustment for RNA-Seq Count Data](https://doi.org/10.1101/2020.01.13.904730).” Preprint. Bioinformatics, January 14, 2020. 

- [combat.py](https://github.com/brentp/combat.py) - python / numpy / pandas / patsy version of ComBat for removing batch effects


## Clustering

- [clust](https://github.com/BaselAbujamous/clust) - Python package for identification of consistently coexpressed clusters of genes, within and across datasets. Consensus clustering principle. Aproximately 50% of genes do not cluster well and thus shouldn't be considered. Compared with seven tools (cross-clustering, k-means, SOMs, MCL, HC, Click, WGCNA) using seven different cluster validation metrics. Outperforms all, produces more focused and significant functional enrichment results. 
    - Abu-Jamous, Basel, and Steven Kelly. “[Clust: Automatic Extraction of Optimal Co-Expressed Gene Clusters from Gene Expression Data](https://doi.org/10.1186/s13059-018-1536-8).” Genome Biology 19, no. 1 (December 2018)

## Timecourse

- [LPWC](https://gitter-lab.github.io/LPWC/) - Lag Penalized Weighted Correlation, a similarity measure to group pairs of time series that are not perfectly synchronized. Review of previous approaches (hierarchical clustering, partition-based, Bayesian models). Correlation-based, with the lag penalty for shift, two options to select it. Best for 5 or more time points, for shorter time course - either use high penalty or another tool STEM. Tested on simulated data (ImpulsDE data) using the adjusted Rand index
    - Chandereng, Thevaa, and Anthony Gitter. “[Lag Penalized Weighted Correlation for Time Series Clustering](https://doi.org/10.1186/s12859-019-3324-1).” BMC Bioinformatics 21, no. 1 (December 2020)

- Time course gene expression analysis review. Biological scenarios requiring a time-course, analytical approaches, Table 1 - software for time course analysis (EDGE, BETR, clustering tools, network analysis).
    - Bar-Joseph, Ziv, Anthony Gitter, and Itamar Simon. “[Studying and Modelling Dynamic Biological Processes Using Time-Series Gene Expression Data](https://doi.org/10.1038/nrg3244).” Nature Reviews. Genetics 13, no. 8 (July 18, 2012)

- [DREM 2.0](http://sb.cs.cmu.edu/drem/) - time course analysis of gene expression data. Detects patterns of gene expression changes and the corresponding transcription factors driving them, motif discovery using protein-DNA (ChIP-seq, ChIP-chip, computational) data, differential motif analysis (DECOD method). Hidden Markov Model-based algorithm. Java tool, GUI and command line interface. 
    - Schulz, Marcel H, William E Devanny, Anthony Gitter, Shan Zhong, Jason Ernst, and Ziv Bar-Joseph. “[DREM 2.0: Improved Reconstruction of Dynamic Regulatory Networks from Time-Series Expression Data](https://doi.org/10.1186/1752-0509-6-104).” BMC Systems Biology 6, no. 1 (2012)


## Differential expression

- [Topconfect](https://bioconductor.org/packages/release/bioc/html/topconfects.html) - differential gene expression using confidence intervals of log fold changes. Confect units are interpretable as effect sizes. Uses TREAT (limma) functionality, Methods, also wraps edgeR and DESeq2 functionality. Tested using simulated and tumor-normal data. R package, [GitHub](https://github.com/pfh/topconfects), [Bioconductor](https://bioconductor.org/packages/release/bioc/html/topconfects.html)
    - Harrison, Paul F., Andrew D. Pattison, David R. Powell, and Traude H. Beilharz. “[Topconfects: A Package for Confident Effect Sizes in Differential Expression Analysis Provides a More Biologically Useful Ranked Gene List](https://doi.org/10.1186/s13059-019-1674-7).” Genome Biology 20, no. 1 (December 2019)

- [Degust](http://degust.erc.monash.edu/) -  interactive RNA-seq analysis, RNA-seq exploration, analysis and visualisation

- [ideal](http://bioconductor.org/packages/release/bioc/html/ideal.html) - Interactive Differential Expression AnaLysis. 

- [EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano) - Publication-ready volcano plots with enhanced colouring and labeling, by Kevin Blighe

- [ALDEx](https://bioconductor.org/packages/release/bioc/html/ALDEx2.html) - differential abundance analysis of compositional data. ANOVA-like approach, partitioning within-condition to between-condition variation. Uses a Dirichlet-multinomial model to infer abundance from counts. Methods description for transforming proportional data into independent components
    - Fernandes, Andrew D., Jean M. Macklaim, Thomas G. Linn, Gregor Reid, and Gregory B. Gloor. “[ANOVA-like Differential Expression (ALDEx) Analysis for Mixed Population RNA-Seq](https://doi.org/10.1371/journal.pone.0067019).” PloS One 8, no. 7 (2013)


## Functional enrichment

- [GSEABenchmarking](https://github.com/waldronlab/GSEABenchmarking) - an R package for systematic testing of gene set enrichment analyses. 10 major enrichment analyses tested on runtime, % significant sets, type I error rate, relevance to phenotype. Microarray and TCGA RNA-seq data. Best performing - overrepresentation analysis, aka Fisher's, hypergeometric test. [Tweet by Levi Waldron](https://twitter.com/LeviWaldron1/status/1142092301403115521?s=03)
    - Geistlinger, Ludwig, Gergely Csaba, Mara Santarelli, Marcel Ramos, Lucas Schiffer, Charity W Law, Nitesh Turaga, et al. “[Towards a Gold Standard for Benchmarking Gene Set Enrichment Analysis](https://doi.org/10.1101/674267).” Preprint. Bioinformatics, June 19, 2019. 

- [PaintOmics 3](http://www.paintomics.org/) - web tool for KEGG pathway enrichment analysis and visualization of gene expression (also, metabolite, protein, region-based data) over pathway diagrams. Competitors: MapMan, KaPPA-View, Pathview Web. Auto-detection of IDs. Analyzes fold change, time course
    - Hernández-de-Diego, Rafael, Sonia Tarazona, Carlos Martínez-Mira, Leandro Balzano-Nogueira, Pedro Furió-Tarí, Georgios J. Pappas, and Ana Conesa. “[PaintOmics 3: A Web Resource for the Pathway Analysis and Visualization of Multi-Omics Data](https://doi.org/10.1093/nar/gky466).” Nucleic Acids Research 46, no. W1 (July 2, 2018)

- [CEMiTool](https://bioconductor.org/packages/release/bioc/html/CEMiTool.html) - gene co-expression analysis, reimplements WGCNA, includes selection of a soft-thresholding power using Cauchi distribution, gene enrichment analysis and, optionally, PPI network. Good overview of WGCNA algorithm. 
    - Russo, Pedro S. T., Gustavo R. Ferreira, Lucas E. Cardozo, Matheus C. Bürger, Raul Arias-Carrasco, Sandra R. Maruyama, Thiago D. C. Hirata, et al. “[CEMiTool: A Bioconductor Package for Performing Comprehensive Modular Co-Expression Analyses](https://doi.org/10.1186/s12859-018-2053-1).” BMC Bioinformatics 19, no. 1 (20 2018)

- [EnrichmentBrowser](https://bioconductor.org/packages/release/bioc/html/EnrichmentBrowser.html) - R package for microarray/RNA-seq normalization, ID mapping, differential analysis, functional enrichment (many methods) and network analyses and visualization
     - Geistlinger, Ludwig, Gergely Csaba, and Ralf Zimmer. “[Bioconductor’s EnrichmentBrowser: Seamless Navigation through Combined Results of Set- & Network-Based Enrichment Analysis](https://doi.org/10.1186/s12859-016-0884-1).” BMC Bioinformatics 17 (January 20, 2016)

- `gProfileR` - enrichment of gene lists, from GO to KEGG and others, organism-specific, ortholog conversion

- [Metascape](http://metascape.org/gp/index.html#/main/step1) - multi-gene lists functional analysis, auto gene ID recognition, >40 databases. [data_analysis_portals.xlsx]((https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09234-6/MediaObjects/41467_2019_9234_MOESM3_ESM.xlsx)) - 25 data analysis portals
    - Zhou, Yingyao, Bin Zhou, Lars Pache, Max Chang, Alireza Hadj Khodabakhshi, Olga Tanaseichuk, Christopher Benner, and Sumit K. Chanda. “[Metascape Provides a Biologist-Oriented Resource for the Analysis of Systems-Level Datasets](https://doi.org/10.1038/s41467-019-09234-6).” Nature Communications 10, no. 1 (December 2019)

- [Positional Gene Enrichment analysis](http://silico.biotoul.fr/pge/) of gene sets for high resolution identification of overrepresented chromosomal regions.  
    - [De Preter, K., Barriot, R., Speleman, F., Vandesompele, J., Moreau, Y., 2008, Nucleic Acids Res.](https://academic.oup.com/nar/article/36/7/e43/2410597)

### Transcription regulators

- [ChEA3](https://amp.pharm.mssm.edu/chea3/) - predicting regulatory TFs for sets of user-provided genes. Improved backend reference gene set data (six datasets), ranking of the most significantly enriched TFs. Benchmarking against several other TF prioritization tools (overviewed in intro). [Docker](https://hub.docker.com/r/maayanlab/chea3), [API web-interface and downloadable data](https://amp.pharm.mssm.edu/chea3/)
    - Keenan, Alexandra B., Denis Torre, Alexander Lachmann, Ariel K. Leong, Megan L. Wojciechowicz, Vivian Utti, Kathleen M. Jagodnik, Eryk Kropiwnicki, Zichen Wang, and Avi Ma’ayan. “[ChEA3: Transcription Factor Enrichment Analysis by Orthogonal Omics Integration](https://doi.org/10.1093/nar/gkz446).” Nucleic Acids Research, May 22, 2019. 

- [RABIT](http://rabit.dfci.harvard.edu/) - find TFs regulating a list of genes. Integrated ChIP-seq and gene expression data, regression framework. Tested in experimental KO data, tumor-profiling cohorts. 
    - Jiang, Peng, Matthew L. Freedman, Jun S. Liu, and Xiaole Shirley Liu. “[Inference of Transcriptional Regulation in Cancers](https://doi.org/10.1073/pnas.1424272112).” Proceedings of the National Academy of Sciences 112, no. 25 (June 23, 2015)

- [RcisTarget](https://github.com/aertslab/RcisTarget) - finding enriched motifs in cis-regulatory regions in a gene list

## Non-canonical RNAs

### Alternative splicing

- Benchmarking review, two types of alternative splicing analysis: differential splicing and differential isoform detection. DESeq2, DEXSeq, Limma and NOISeq perform well overall. [GitHub](https://github.com/gamerino/benchmarkingDiffExprAndSpl)
    - Merino, Gabriela A, Ana Conesa, and Elmer A Fernández. “[A Benchmarking of Workflows for Detecting Differential Splicing and Differential Expression at Isoform Level in Human RNA-Seq Studies](https://doi.org/10.1093/bib/bbx122).” Briefings in Bioinformatics 20, no. 2 (March 25, 2019)

- [MAJIQ](https://majiq.biociphers.org/) - local splicing variation analysis. Detects canonical and alternative splicing events. Quantifies as Percent Selected In (PSI). Differential splicing as delta PSI. Visualization using VOILA package. Python 3.
    - Vaquero-Garcia, Jorge, Alejandro Barrera, Matthew R. Gazzara, Juan González-Vallinas, Nicholas F. Lahens, John B. Hogenesch, Kristen W. Lynch, and Yoseph Barash. “[A New View of Transcriptome Complexity and Regulation through the Lens of Local Splicing Variations](https://doi.org/10.7554/eLife.11752).” ELife 5 (February 1, 2016)

- [MISO](https://miso.readthedocs.io/en/fastmiso/) (Mixture-of-Isoforms) is a probabilistic framework that quantitates the expression level of alternatively spliced genes from RNA-Seq data, and identifies differentially regulated isoforms or exons across samples. - By modeling the generative process by which reads are produced from isoforms in RNA-Seq, the MISO model uses Bayesian inference to compute the probability that a read originated from a particular isoform.- MISO treats the expression level of a set of isoforms as a random variable and estimates a distribution over the values of this variable. - The estimation algorithm is based on sampling, and falls in the family of techniques known as Markov Chain Monte Carlo (“MCMC”). 
    - Katz, Yarden, Eric T. Wang, Edoardo M. Airoldi, and Christopher B. Burge. “[Analysis and Design of RNA Sequencing Experiments for Identifying Isoform Regulation](https://doi.org/10.1038/nmeth.1528).” Nature Methods 7, no. 12 (December 2010)

- [RegTools](https://regtools.readthedocs.io/en/latest/) - integration of somatic variants from DNA-seq and splice junctions from RNA-seq data to identify variants causing aberrant splicing in cancer. 
    - Feng, Yang-Yang, Avinash Ramu, Kelsy C Cotto, Zachary L Skidmore, Jason Kunisaki, Donald F Conrad, Yiing Lin, et al. “[RegTools: Integrated Analysis of Genomic and Transcriptomic Data for Discovery of Splicing Variants in Cancer](https://doi.org/10.1101/436634),” November 25, 2018. 

- [rMATS](http://rnaseq-mats.sourceforge.net/) alternative splicing detection tool. Using paired samples.RNA-seq depth and alternative splicing power - 200M reads minimum. [rMATS-turbo](https://github.com/Xinglab/rmats-turbo/) - 100X faster implementation, [Tweet](https://twitter.com/YiXing77/status/1267834667715235843?s=20)
    - Shen, Shihao, Juw Won Park, Zhi-xiang Lu, Lan Lin, Michael D. Henry, Ying Nian Wu, Qing Zhou, and Yi Xing. “[RMATS: Robust and Flexible Detection of Differential Alternative Splicing from Replicate RNA-Seq Data](https://doi.org/10.1073/pnas.1419161111).” Proceedings of the National Academy of Sciences of the United States of America 111, no. 51 (December 23, 2014)

- [tappAS](http://tappas.org/) - functional impact of alternative splicing. Input - transcript-level count matrix

- [vast-tools](https://github.com/vastgroup/vast-tools) - A toolset for profiling alternative splicing events in RNA-Seq data
    - Irimia, Manuel, Robert J. Weatheritt, Jonathan D. Ellis, Neelroop N. Parikshak, Thomas Gonatopoulos-Pournatzis, Mariana Babor, Mathieu Quesnel-Vallières, et al. “[A Highly Conserved Program of Neuronal Microexons Is Misregulated in Autistic Brains](https://doi.org/10.1016/j.cell.2014.11.035).” Cell 159, no. 7 (December 18, 2014)

### miRNAs

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

- [lncRNAKB](http://psychiatry.som.jhmi.edu/lncrnakb/) - database of long noncoding RNAs. lncRNAs are typically less conserved, expressed low on average and highly tissue-specific. Combines six resources (CHESS, LNCipedia, NONCODE, FANTOM, MiTranscriptome, BIGTranscriptome). Information about tissue-specific expression, eQTL, WGCNA co-expression to predict functions in a tissue-specific manner, random forest prediction of protein-coding score. Data: GTF gene annotation, tissue-specific expression (TPM, counts, eQTL). [RNA-seq blog post](https://www.rna-seqblog.com/lncrnakb-a-comprehensive-knowledgebase-of-long-non-coding-rnas/)
    - Seifuddin, Fayaz, Komudi Singh, Abhilash Suresh, Yun-Ching Chen, Vijender Chaitankar, Ilker Tunc, Xiangbo Ruan, et al. “[LncRNAKB: A Comprehensive Knowledgebase of Long Non-Coding RNAs](https://doi.org/10.1101/669994).” Preprint. Bioinformatics, June 13, 2019. 

- [UClncR](http://bioinformaticstools.mayo.edu/research/uclncr-pipeline/) - detecting and quantifying expression of unknown and known lncRNAs. Works for unstranded and stranded RNA-seq. Incorporates StringTie, Sebnif for novel lncRNA detection, iSeeRNA for assessing noncoding potential. Annotates lncRNAs by the nearby protein-coding genes. Tested on real data using Gencode annotations with parts of lncRNA annotations removed. 
    - Sun, Zhifu, Asha Nair, Xianfeng Chen, Naresh Prodduturi, Junwen Wang, and Jean-Pierre Kocher. “[UClncR: Ultrafast and Comprehensive Long Non-Coding RNA Detection from RNA-Seq](https://doi.org/10.1038/s41598-017-14595-3).” Scientific Reports 7, no. 1 (December 2017)


### circRNAs

- [CIRCpedia](http://www.picb.ac.cn/rnomics/circpedia/) database of cornRNAs from human, mouse, and some model organisms. Ribo-, poly(A)-, RNAse R methods for enriching for circRNAs. [CIRCexplorer2](https://circexplorer2.readthedocs.io/en/latest/) for the analysis of such experiments
    - Zhang et al., “Diverse Alternative Back-Splicing and Alternative Splicing Landscape of Circular RNAs.”

### Gene fusion

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
    - Flensburg, Christoffer, Alicia Oshlack, and Ian J Majewski. “[Detecting Copy Number Alterations in RNA-Seq Using SuperFreq](https://www.biorxiv.org/content/10.1101/2020.05.31.126888v1)” n.d., 22.

- [CaSpER](https://github.com/akdess/CaSpER) - identification of CNVs from  RNA-seq data, bulk and single-cell (full-transcript only, like SMART-seq). Utilized multi-scale smoothed global gene expression profile and B-allele frequency (BAF) signal profile, detects concordant shifts in signal using a 5-state HMM (homozygous deletion, heterozygous deletion, neutral, one-copy-amplification, high-copy-amplification). Reconstructs subclonal CNV architecture for scRNA-seq data. Tested on GBM scRNA-seq, TCGA, other. Compared with HoneyBADGER. R code and tutorials
    - Serin Harmanci, Akdes, Arif O. Harmanci, and Xiaobo Zhou. “[CaSpER Identifies and Visualizes CNV Events by Integrative Analysis of Single-Cell or Bulk RNA-Sequencing Data.](https://doi.org/10.1038/s41467-019-13779-x)” Nature Communications 11, no. 1 (December 2020)

- [CNVkit-RNA](https://github.com/etal/cnvkit) - CNV estimation from RNA-seq data. Improved moving average approach, corrects for GC content, gene expression level, gene length, correlation of gene expression and CNV (estimated from TCGA). [Docs](https://cnvkit.readthedocs.io/en/stable/), [Video tutorial](https://youtu.be/JGOIXYoLG6w)
    - Talevich, Eric, and A. Hunter Shain. “[CNVkit-RNA: Copy Number Inference from RNA-Sequencing Data.](https://doi.org/10.1101/408534)” Preprint. Bioinformatics, September 4, 2018

- [InferCNV](https://www.bioconductor.org/packages/release/bioc/html/infercnv.html) - Inferring copy number alterations from tumor single cell RNA-Seq data. R package. [GitHub wiki](https://github.com/broadinstitute/inferCNV/wiki). Part of [Trinity Cancer Transcriptome Analysis Toolkit](https://github.com/NCIP/Trinity_CTAT/wiki)

- [SQUID](https://github.com/Kingsford-Group/squid) - transcriptomic structural variation caller. Genome segment graph, then rearrange segments so that as many read alignments as possible are concordant with the rearranged sequence. Compared with MUMmer3, DELLY2, LUMPY in simulated settings, and with SOAPfuse, deFuse, FusionCatcher, JAFFA, INTEGRATE tools using real data
    - Ma, Cong, Mingfu Shao, and Carl Kingsford. “[SQUID: Transcriptomic Structural Variation Detection from RNA-Seq](https://doi.org/10.1186/s13059-018-1421-5).” Genome Biology 19, no. 1 (12 2018)

- [transindel](https://github.com/cauyrd/transIndel) - Indel caller for DNA-seq or RNA-seq


## Networks

- [GENIE3](https://github.com/aertslab/GENIE3) - random forest regression detection of gene modules. Input - expression matrix, output - gene x gene square co-regulation matrix

## Integrative

- [DIABLO](http://mixomics.org/) - multi-omics analysis method. Overview of previous methods (SNF, Bayesian Consensus Clustering, NMF, JIVE, sGCCA, MOFA, others). Method extends sGCCA multivariate dimensionality reduction that uses SVD and selects co-expressed (correlated) variables from several omics datasets. Methods, model, iterative solution. Design matrix specifies which omics datasets are connected. Variable selection for biomarkers identification.  Visualization options. Part of [mixOmics R package](http://mixomics.org/), [Documentation](https://mixomicsteam.github.io/Bookdown/intro.html)
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

- [chromoMap](https://cran.r-project.org/web/packages/chromoMap/index.html) - Interactive Visualization and Mapping of Human Chromosomes, R package

- [genomation](https://bioconductor.org/packages/devel/bioc/vignettes/genomation/inst/doc/GenomationManual.html) - a toolkit for annotation and visualization of genomic data, R package

- [karyoploteR](https://github.com/bernatgel/karyoploteR) - An R/Bioconductor package to plot arbitrary data along the genome

- [pcaExplorer](https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html) - Interactive Visualization of RNA-seq Data Using a Principal Components Approach, R package

- `WIlsON` - Web-based Interactive Omics VisualizatioN, accepts, text files, SummarizedExperiment datasets. 


## Data

- [refine.bio](https://www.ccdatalab.org/projects/refinebio) - harmonized microarray and RNA-seq data for various organisms and conditions

- [cBioPortalData](https://github.com/waldronlab/cBioPortalData) - cBioPortal data as MultiAssayExperiment objects, by Waldron Lab

- [curatedTCGAData](https://github.com/waldronlab/curatedTCGAData) - Curated Data From The Cancer Genome Atlas (TCGA) as MultiAssayExperiment objects, by Waldron Lab

- [ARCHS4](http://amp.pharm.mssm.edu/archs4/download.html) - Massive Mining of Publicly Available RNA-seq Data from Human and Mouse

- [PINS](http://www.cs.wayne.edu/tinnguyen/PINS/PINS.html) - A novel method for data integration and disease subtyping

- [recount2](https://bioconductor.org/help/workflows/recountWorkflow/) - an R workflow to work with recount2 data

- [GREIN](http://www.ilincs.org/apps/grein/) - re-analysis of RNA-seq datasets from GEO. Download processed data, visualization, power analysis, differential expression, functional enrichment analysis, connectivity analysis with LINCS L1000 data. [GitHub](https://github.com/uc-bd2k/grein), [Docker image](https://hub.docker.com/r/ucbd2k/grein/)
    - Al Mahi, Naim, Mehdi Fazel Najafabadi, Marcin Pilarczyk, Michal Kouril, and Mario Medvedovic. “[GREIN: An Interactive Web Platform for Re-Analyzing GEO RNA-Seq Data](https://doi.org/10.1101/326223),” October 27, 2018

- [DEE2](http://dee2.io) - Digital Expression Explorer - gene- and transcript-level processed data from multiple organisms, amenable for downstream analysis in R etc. [getDEE2](https://github.com/markziemann/dee2/blob/master/AccessDEEfromR.md) R package to get the data
    - Ziemann, Mark, Antony Kaspi, and Assam El-Osta. “[Digital Expression Explorer 2: A Repository of Uniformly Processed RNA Sequencing Data](https://doi.org/10.1093/gigascience/giz022).” GigaScience 8, no. 4 (April 1, 2019)

- [GEMMA](https://gemma.msl.ubc.ca/home.html) - curated transcriptomic database, >10,000 studies, \~34% are brain-related. Query genes, phenotypes, experiments, search for coexpression, differential expression. Processing methods, batch correction. Online access, API, R package. [GitHub](https://github.com/PavlidisLab/Gemma/)
    - Lim, Nathaniel, Stepan Tesar, Manuel Belmadani, Guillaume Poirier-Morency, Burak Ogan Mancarci, Jordan Sicherman, Matthew Jacobson, Justin Leong, Patrick Tan, and Paul Pavlidis. “[Curation of over 10,000 Transcriptomic Studies to Enable Data Reuse](https://doi.org/10.1101/2020.07.13.201442).” Preprint. Bioinformatics, July 14, 2020. 

### Genes

- [List of gene lists for genomic analyses](https://github.com/macarthur-lab/gene_lists) - GitHub repo with tab-separated annotated lists

- [CREEDS](https://amp.pharm.mssm.edu/creeds/) - database of manually (and automatically) extracted gene signatures. Single gene perturbations, disease signatures, single drug perturbations. Batch effect correction, when necessary. Overall, good agreement with MSigDb C2. Characteristic Direction (CD) method to detect differential genes. [API access in R](http://rpubs.com/wangz10/177826)
    - Wang, Zichen, Caroline D. Monteiro, Kathleen M. Jagodnik, Nicolas F. Fernandez, Gregory W. Gundersen, Andrew D. Rouillard, Sherry L. Jenkins, et al. “[Extraction and Analysis of Signatures from the Gene Expression Omnibus by the Crowd](https://doi.org/10.1038/ncomms12846).” Nature Communications 7, no. 1 (November 2016)

- [OGEE v3](https://v3.ogee.info/#/home) - online gene essentiality database. Experimentally validated using CRISPR and RNAi. All data (essential genes, their properties, expression data and more) are downloadable 
    - Gurumayum, Sanathoi, Puzi Jiang, Xiaowen Hao, Tulio L Campos, Neil D Young, Pasi K Korhonen, Robin B Gasser, et al. “[OGEE v3: Online GEne Essentiality Database with Increased Coverage of Organisms and Human Cell Lines](https://doi.org/10.1093/nar/gkaa884)”

- [HGNChelper](https://waldronlab.io/HGNChelper/) - Identify and Correct Invalid HGNC Human Gene Symbols and MGI Mouse Gene Symbols

- [gencode_regions](https://github.com/saketkc/gencode_regions) - Extract 3'UTR, 5'UTR, CDS, Promoter, Genes, Introns etc from GTF files, by Saket Choudhary

- [Extract intron boundaries per transcript](https://gist.github.com/hiraksarkar/ce8a71a6953cb4e9823d868c283bf99d)

- [CHESS](http://ccb.jhu.edu/chess/) - Comprehensive Human Expressed SequenceS, database of novel genes, identified from GTeX data, protein-coding and lncRNA

- [GTEx Visualizations](https://github.com/broadinstitute/gtex-viz) - web-based visualization tools for exploring tissue-specific gene expression and regulation

- [ideogram](https://github.com/eweitz/ideogram) - Chromosome visualization with D3.js. [Examples](https://eweitz.github.io/ideogram/). [ideogRam R wrapper](https://github.com/freestatman/ideogRam)

- [Enhancer-promoter (EP) pairs from Thurman et al., (2012)](ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz)

- [pysradb](https://github.com/saketkc/pysradb) - Python package for interacting with SRAdb and downloading datasets from SRA, by Saket Choudhary. [Documentation](https://www.saket-choudhary.me/pysradb/)


## Misc

- [Recommended Coverage and Read Depth for NGS Applications](https://genohub.com/recommended-sequencing-coverage-by-application/), by GenoHub. And, their [NGS Handbook](https://genohub.com/next-generation-sequencing-handbook/)

- [BioJupies](https://amp.pharm.mssm.edu/biojupies/) - analysis of GEO/GTEx data or your own gene expression table/FASTQ in autogenerated Jupyter notebook. Rich set of tools for EDA (PCA, Clustergrammer, Library size analysis), Differential expression analysis (Volcano, MA plots), Enrichment analysis (Enrichr, GO, Pathway, TF, Kinase, miRNA enrichments), L1000 signatures. Best suited for two-group analysis. Includes Methods for the selected tools

- [HGNChelper](https://cran.r-project.org/web/packages/HGNChelper/index.html) - Handy Functions for Working with HGNC Gene Symbols and Affymetrix Probeset Identifiers

- `tximport` - importing transcript abundance datasets from Salmon, Sailfish, kallisto, RSEM, and differential analysis

- [rpkmforgenes](http://sandberg.cmb.ki.se/media/data/rnaseq/instructions-rpkmforgenes.html) - a Python script for calculating gene expression for RNA-Seq data

- [Python interface to access reference genome features from Ensembl](https://github.com/openvax/pyensembl), e.g., genes, transcripts, and exons

- [TPMcalculator](https://github.com/ncbi/TPMCalculator) - converts gene counts to TPM using transcript information from a GTF file. TPM vs. FPKM correlation for validation. C/C++ command line tool, Docker image, CWL workflow
    - Vera Alvarez, Roberto, Lorinc Sandor Pongor, Leonardo Mariño-Ramírez, and David Landsman. “[TPMCalculator: One-Step Software to Quantify MRNA Abundance of Genomic Features](https://doi.org/10.1093/bioinformatics/bty896).” Bioinformatics 35, no. 11 (June 1, 2019)

