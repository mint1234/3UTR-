priUTR pipeline

This pipeline implements methods to analyze alternative 3’ untranslated regions (3’UTRs) using standard RNA-seq data in non-model plant species. 

Authors

Min Tu, Yin Li, Yaping Feng, Dibyendu Kumar
Waksman Institute of Microbiology, Rutgers, The State University of New Jersey, USA.

Citation

Please cite the following artile when using “priUTR”:
Min Tu, Yin Li, Yaping Feng, Dibyendu Kumar. (2020) priUTR pipeline identifies alternative 3’UTRs and provides insights into their regulation in maize endosperm. (in submission)

License

GNU General Public License v3.0

1.	Introduction of the priUTR pipeline -- Detect alternative 3' UTR in plants using RNA-seq data 

priUTR pipeline is developed in Scientific Linux release 6.10 (Carbon). This pipeline is based on Linux shell scripts, it requires a 64-bit CPU computer running on Linux system. 5-50GB running RAM is recommended depending on data size. To make full use of standard RNA-seq data and to provide user-friendly 3’UTR analysis function for many non-model plant species and plant biologists who may have limited bioinformatic skills, the priUTR pipeline can detect alternative 3’UTRs in comparison to the reference genome annotation and provide information regarding the identified alternative 3’UTRs. 
2.	Workflow of priUTR and upstream RNA-seq analysis pipeline
 
Figure 1. Workflow of the priUTR and upstream RNA-seq analysis pipeline. 


2.  Environment and dependencies

-- samtools Version 1.3.1, http://www.htslib.org 
-- bedtools v2.24.0, https://bedtools.readthedocs.io/en/latest/ 
-- TransDecoder Version 2.0.1, https://github.com/TransDecoder/TransDecoder/releases 
-- gffread Version 0.9.8c, https://github.com/gpertea/gffread 
-- ViennaRNA Version 2.4.14, https://www.tbi.univie.ac.at/RNA/RNAfold.1.html 
-- ToGenePred: https://bioconda-recipes-demo.readthedocs.io/en/docs/recipes/ucsc- gtftogenepred/README.html 
-- exteractFromFasta.pl: https://github.com/jonbra/NGSAbel/blob/master/scripts/extractFromFasta.pl 
-- qsub Version 3.0.4: http://docs.adaptivecomputing.com/torque/2-5- 12/help.htm#topics/commands/qsub.htm 

3.	Input and Output

Example input and output files can be downloaded from Waksman Institute of Microbiology, Rugters, The State University of New jersey (https://data.waksman.rutgers.edu/200106-3utr/).

4.	Usage

4.1	Download priUTR pipeline from the following Github website (https://github.com/mint1234/3UTR-).

4.2	Create two folders, namely “inputs” and “softwares”. Install all of the above-mentioned software (see “Environment and dependencies”) in the “softwares” folder. 

4.3	Then move the following files into the “inputs” folder: (1) genome.gff3 (genome annotation); (2) genome.fa (assembled genome sequence); (3) the output files from RNA-seq read aligner and transcriptome assembler (Tophat and Cufflinks, etc.), including three “.bam” files corresponding to the triplicated RNA-seq samples respectively generated by Tophat, isoforms_fpkm.tracking (Cuffdiff output file to show isoform expression); merged.gtf (Cuffmerge output file of a transcript assembly merged from multiple RNA-seq samples). Note: example input files can be downloaded from: https://data.waksman.rutgers.edu/200106-3utr/. 

4.4	Please set local path to “$PRE” in the shell file, and then execute the “genome.sh” file using command “sh genome.sh”; The script “genome.sh” is to extract 3’UTR information of all annotated transcripts in the reference genome. 

4.5	Please set local path to “$PRE” in the ”qsub-3r.sh“ and ”3r.sh“ files. Then set condition names to “${SN}”, and set how many RNA-seq samples (or experiment treatments or conditions) to “i” in the ”qsub-3r.sh“ file. The priUTR pipeline can automatically loop execute for “i” times to proceed multiple RNA- seq samples (“i” samples, in this case) in one run. In the “3r.sh” file, the default setting of RPKM filter is “A=1”, and the default setting of 3’UTR-reads filter is “B=30”, they can be set or adjusted in the pipeline manually. Then execute the “qsub-3r.sh” using command “sh qsub-3r.sh” (which loop executes “3r.sh”). 

4.6	The output results will be generated in two folders: “results” and “2nd-SN” folder. The “results” folder will contain three result files for RNA-seq sample: (1) the “SN-alter-3UTR.fasta” file is a fasta format file containing all the alternative 3’UTR sequences identified; (2) the “SN.change.txt” file contains the following information in each column: 1) TCONS ID: assembled transcript ID from Cufflinks;2) transcripts ID: corresponding transcript ID in the annotation; 3) gene ID: corresponding gene ID in the annotation; 4) class code: the type of match between the Cufflinks transcripts and the reference transcript. A full description of class code meaning can be obtained from the link: http://cole-trapnell-lab.github.io/cufflinks/ cuffcompare/ . 5) alternative 3UTR length; 6) annotated 3UTR length; 7) transdecoder; 8)+/- transdecoder; 9) 3utr length transdecoder: columns 7, 8, and 9 represent the source, strand and length of predicted 3’UTR using transdecoder; 10) FPKM-replicate-1; 11) FPKM-replicate-2; 12) FPKM-replicate-3: columns 10, 11, 12 represent the expression levels of the alternative 3’UTR for RNA-seq replicate1, 2, and 3, respectively. (3) the “SN-R1/2/3-3utr-counts-30.txt” file contains the following information for each column: 1) chromosome: chromosomal location of the alternative 3’UTR; 2) alternative 3utr_start; 3) alternative 3utr_end: columns 2 and 3 represent the start and end positions of the alternative 3’UTR. 4) transcript ID: the corresponding transcript ID in the annotation; 5) addition information; 6) +/-; 7) number of reads mapped to the 3'UTR region, those transcripts with less than 30 reads have been filtered out; 

4.7	The “2nd-SN” folder contains all the secondary structure plots that are predicted for the priUTR- identified alternative 3’UTR sequences by RNAfold from the ViennaRNA package with minimum free energy. For the predicted secondary structures of alternative 3’UTRs, please also see the example output files available at Waksman Institute (folder namely “2nd- Endosperm_16DAP” and “2nd-Endosperm_20DAP” at link: https://data.waksman.rutgers.edu/200106-3utr/) 


5.	Issues and bug reports

Please use https://github.com/mint1234/3UTR-/issues to submit issues, bug reports, and comments.

For further information, please contact with mintu@waksman.rutgers.edu.

6.	References

Trapnell C, Roberts A, Goff L, Pertea G, Kim D, Kelley DR, Pimentel H, Salzberg SL, Rinn JL, Pachter L. (2012) Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. Nature Protocol. 7(3): 562-578. 

Trapnell C, Hendrickson DG, Sauvageau M, Goff L, Rinn JL, Pachter L. (2013) Differential analysis of gene regulation at transcript resolution with RNA-seq. Nature Biotechnology. 31, 46–53. 

Tophat website: http://ccb.jhu.edu/software/tophat/index.shtml 

Cufflinks website: http://cole-trapnell-lab.github.io/cufflinks/manual/ 

Cufflinks github website: https://github.com/cole-trapnell-lab/cufflinks

Wu X, Liu M, Downie B, Liang C, Ji G, Li QQ, Hunt AG. (2011) Genome-wide landscape of polyadenylation in Arabidopsis provides evidence for extensive alternative polyadenylation. Proc. Natl. Acad. Sci. USA，108, 12533- 12538. 

Jafar Z, Tarig S, Sadiq I, Nawaz T, Akhtar MN. (2019) Genome-Wide Profiling of Polyadenylation Events in Maize Using High-Throughput Transcriptomic Sequences. G3: Genes, Genomes, Genetics. 9, 2749-2760. 

Ye C, Long Y, Ji G, Li QQ, Wu X. (2018) APAtrap: identification and quantification of alternative polyadenylation sites from RNA-seq data. Bioinformatics, 34, 1841-1849. 

Arefeen A, Liu J, Xiao X, Jiang T. (2018) TAPAS: tool for alternative polyadenylation site analysis. Bioinformatics, 34, 2521- 2529. 
