### gene2probe

Probe-based spatial transcriptomics technologies offer great flexibility in terms of the transcripts that can be profiled. 
While technologies such as VisiumHD and VisiumFFPE enable the capture of most protein-coding genes, there are many applications that might require the design of custom probes.
A few examples include:

1. Adding probes for genes not currently included in the probe-set (e.g., XIST).
2. Adding probes to distinguish between specific isoforms of the same gene.
3. Adding more probes to increase the detection of a gene of specific interest in a particular experiment.

Such probes often need to fulfill specific requirements (GC content, specific nucleotides in specific positions, avoid polymorphism and repeats) and need to be specific for their targets. 
gene2probe aims to help the user to design such probes. The default parameters are tailored towards the current recommendations by 10x Genomics for VisiumHD and VisiumFFPE, but in principle gene2probe can be used to design probes of any length and nucleotide requirement.

We always strongly recommend additionally manually BLASTing the selected probes before proceeding with ordering them.

We also note that gene2probe is currently tailored towards designing probes for the same species as those covered by the core probe set (and provided input files specifically cover the human genome). 
While the design of custom probes against bacterial or viral genomes is an exciting application, they pose many additional considerations that gene2probe doesn't currently take into account.
