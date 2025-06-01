This repository contains two scripts designed for detecting and visualizing Loss of Heterozygosity (LOH) across the whole genome using individual-level SNP data from gVCF files.

1. detect_LOH_from_individual_gvcf.pl
Function:
Detects genome-wide LOH regions based on SNP data from a single individual's gVCF file. This script does not rely on joint-called VCFs and is suitable for analyzing each sample independently.

Input:

A gVCF file for a single individual
Optional parameters for window size, heterozygosity thresholds, etc.
Output:

A tab-delimited file listing detected LOH regions with chromosome, position, and metrics
Can be used directly as input for visualization

2. plot_genomewide_LOH.R
Function:
Generates a genome-wide visualization of LOH regions, based on the output from the Perl script. Useful for identifying patterns or clusters of LOH across chromosomes.

Input:

Output file from detect_LOH_from_individual_gvcf.pl

Output:

Manhattan-style plot of LOH across the genome
