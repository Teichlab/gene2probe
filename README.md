# gene2probe

Probe-based spatial transcriptomics technologies offer great flexibility in terms of the transcripts that can be profiled. 
While technologies such as VisiumHD and VisiumFFPE enable the capture of most protein-coding genes, there are many applications that might require the design of custom probes.
A few examples include:

1. Adding probes for genes not currently included in the probe-set (e.g., XIST).
2. Adding probes to distinguish between specific isoforms of the same gene.
3. Adding more probes to increase the detection of a gene of specific interest in a particular experiment.

Such probes often need to fulfill specific requirements (GC content, specific nucleotides in specific positions, avoid polymorphism and repeats) and need to be specific for their targets. 
**gene2probe** aims to help the user to design such probes. The default parameters are tailored towards the current [recommendations by 10x Genomics for VisiumHD and VisiumFFPE](https://cdn.10xgenomics.com/image/upload/v1697739385/support-documents/CG000621_CustomProbeDesign_TechNote_RevC.pdf), but in principle **gene2probe** can be used to design probes of any length and nucleotide requirement.

**We always strongly recommend additionally manually BLASTing the selected probes before proceeding with ordering them.**

We also note that **gene2probe** is currently tailored towards designing probes for the same species as those covered by the core probe set (and provided input files specifically cover the human genome). 
While the design of custom probes against bacterial or viral genomes is an exciting application, they pose many additional considerations that gene2probe doesn't currently take into account.

## Installation

**gene2probe** builds on many awesome bioinformatic tools that need to be preinstalled. Luckily, they should all be easy to install via conda/pip.

We recommend setting up a dedicated conda environment.

```
conda create -n gene2probe_env python=3.11
```

Install bedtools and blast

```
conda activate gene2probe_env
conda install bioconda::blast
conda install bioconda::bedtools
```

Then install gene2probe via pip

```
pip install git+https://github.com/Teichlab/gene2probe.git
```

## Probe design recommendations and main principle

Depending on the assay you are planning to use, probe requirements will differ in terms of length, required GC content, whether the probes are split or not, and required nucleotides in specific positions.

Regardless of the specific requirements, you typically don't want to design probes agains repeats/low complexity regions or parts of the genome with common polymorphism in the species of interest. Additionally, you want to make sure that your probe specifically targets your gene of interest and will not lead to the detection of off-targets.

Our default parameters provided in the tutorials are tailored around the current recommendations by 10x Genomis for VisiumHD. These are the following:
1. Probe length of 50 nts, split at the 25th position.
2. GC content of each probe side should range between 44% and 72%.
3. The 25th nucleotide (last nucleotide of the left hand side probe should be a T).
4. Not overlapping common polymorpshism, repeats or low complexity regions (ideally nowhere in the probe, as a minimal requirement not within 5 nts from the ligation junction).
5. A minimum of 5 mismatches between at least the left or right hand side of the probe and any off-targets.
6. Multiple probes for the same gene should not overlap each other.
7. Prioritise probes with shorter homopolymer lengths.

Due to its modular nature, gene2probe allows the user to finetune these criteria. This gives flexibility in terms of using gene2probe to design probes for other assays, but also allows the user to become more strict or lenient depending on the number of probes that are available for a given gene.

For example, you can start by removing all probes overlapping a common polymorphism in any position - if that leads to too few probes, you can relax the requirement to +/- 5 nts from the ligation junction.

## Resources

We rely on several publicly available resources, most of which can be directly obtained via the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables).

The following resources are required:
1. Genome file in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format). For the provided tutorials, we use the most popular human genome version (hg38), which can be downloaded from [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz). Other genomes can be obtained via the [UCSC genome browser](https://hgdownload.soe.ucsc.edu/downloads.html).
2. Genome annotation in [GTF format](https://genome.ucsc.edu/FAQ/FAQformat.html#format4). There are many different annotation versions per genome, here we recommend the [RefSeq annotation for hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz), which is curated and only includes high-confidence transcripts. 
You typically don't want to design a probe against some obscure exon only found in a rare isoform. However, depending on your experimental design you might want to do just that, so choosing the right GTF file is an essential part of this process. 
More GTF files can be accessed from [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables) and [Ensembl](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/).
3. Gaps in the human genome assembly in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). As long as you stick to hg38, you can use the ones provided in `hg38_resources/hg38_gaps.bed`. Otherwise, you can download the ones for your genome assembly from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables) (select Mapping and Sequencing/Gaps and download in BED format for the whole genome).
4. Repeats and low complexity regions in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). As long as you stick to hg38, you can use the ones provided in `hg38_resources/hg38_rmsk.bed`. Otherwise, you can download the ones for your genome assembly from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables) (select Repeats/RepeatMasker and download in BED format for the whole genome).
5. SNPs and small indels in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). As long as you stick to hg38, you can use the ones provided in `hg38_resources/hg38_snp151Common.bed`. Additional options can be downloaded from from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables) (we used Variation/Common SNPs(151)/snp151Common but you might also want to consider other options).

 **In all cases, please make sure that your annotations match your genome assembly!**

Additionally, you will need one or more blast databases to check for off-target effects. We provide instructions to construct your database in `notebooks/001_make_blast_database.ipynb`.

## Usage and documentation

Please refer to our demo notebooks for examples on how to design your own probes.
We also provide examples on how to make your own blast databases.

**Remember that we always strongly recommend additionally manually BLASTing the selected probes before proceeding with ordering them.**

## Acknowledgements

We are thankful to 10x Genomics, Cecilia Kyanya, Sam Dougan and members of the Wellcome Sanger PAM informatics team for useful discussion. We also acknowledge being aided by ChatGPT4 when writing and documenting this pipeline.
