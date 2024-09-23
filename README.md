# Analysis-T2T
Human-associated metagenomic data often contain human nucleic acid information, which can affect the accuracy of microbial classification or raise ethical concerns. These reads are typically removed by aligning them to the human genome and filtering them out before metagenomic analysis, using various metagenomic mapping tools or human reference genomes. In this study, we conducted a comprehensive analysis to identify the optimal combination of alignment software and human reference genomes using benchmarking data. Our findings show that the combination of bwa-mem and the telomere-to-telomere (T2T) human genome is the most effective at removing human reads in simulated data. We also analyzed T2T-derived sequences in RefSeq to understand how T2T reduces false positive results. Finally, we assessed clinical samples and found that T2T effectively reduces host-derived contamination, particularly in low microbial biomass samples. This study provides a thorough overview of the application of T2T in metagenomic analysis and highlights its significance in improving microbial classification accuracy.

## Description of the scripts

removehost.py:  script for removing host sequences.

get_taxon_name_and_category.py: get taxon name and category.

get_refseq_blast_stats.py: get the reference genome information for BLAST alignment.

get_microbe_report_for_TN_reads.py: get microbe report for TN reads

get_performance_stats.py: get performance data and statistical results.

get_TaxID_for_FN_reads.py: get TaxID for FN reads.

get_GC_stats_v5.py: get_GC_stats: calculate the GC content of different reference genomes

analysis.R: 
The main script in the project includes statistical analysis results and drawing functions.

normality_test.R: Make MCC statistics and draw qqplot for the results of different strategies.

performance.R: 
Performance results by employing different genome as the reference genome using different algorithms.

benchmarking.R: 
draw benchmarking plot

consensus.R: 
draw the consistent results
