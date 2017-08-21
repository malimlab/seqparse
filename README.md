## seqparse
Software to analyse RT termination sites (as reported in Pollpeter et al. 2017)

## Prequisites

bowtie-1.1.2 (in same directory as seqparse)
fastx-0.0.13 (in same directory as seqparse)
samtools (in /usr/local)
bam-readcount (in /usr/local)

## Usage

Place the MiSeq .fastq.gz files into a directory

./parse.sh <directory containing .fastq.gz files>

## Output

A directory named "parse_results" containing:

readCount.csv - counts of RT product lengths in each sample
readCount_normalised.csv - RT product lengths normalised to total reads per sample
readCount_length_adjusted.csv - as above, but using oligo control sample to adjust for bias in read lengths
baseVariation_<sample_name>.csv - Per-base mutation summary

## Files

parse.sh - main script to parse data
parse_sam.pl - takes .sam files and extracts RT product lengths
rc_extract.pl - extracts mutation information from readcount files