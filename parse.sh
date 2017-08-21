#!/bin/bash

echo "Decompressing files..."
gunzip -v $1/*_001.fastq.gz

echo "Collapsing duplicates and trimming adaptors..."

mkdir $1/001_collapsed
mkdir $1/002_trimmed

for f in $1/*_R2_001.fastq
do
	tempname=${f##*/}
	name=`basename $tempname .fastq`
	echo "* $name"
	../fastx-0.0.13/fastx_collapser -Q33 -i $1/$name.fastq -o $1/001_collapsed/$name.fasta
	../fastx-0.0.13/fastx_trimmer -f 41 -i $1/001_collapsed/$name.fasta -o $1/002_trimmed/$name.fasta
done

echo "Aligning with bowtie..."

mkdir $1/003_sam

for f in $1/002_trimmed/*.fasta
do
	tempname=${f##*/}
	name=`basename $tempname .fasta`
	echo "* $name"
	../bowtie-1.1.2/bowtie -v 3 --best --quiet -m 1 -f -S /Users/andrew/Lab\ data/Darja\ sequencing/bowtie-1.1.2/indexes/HIVss $1/002_trimmed/$name.fasta $1/003_sam/$name.sam
done

echo "Generating and sorting BAM files..."

for f in $1/003_sam/*.sam
do
	tempname=${f##*/}
	name=`basename $tempname .sam`
	echo "* $name"
	/usr/local/bin/samtools view -bS -o $1/003_sam/$name.bam $1/003_sam/$name.sam
	/usr/local/bin/samtools sort $1/003_sam/$name.bam $1/003_sam/$name.sorted
	/usr/local/bin/samtools index -b $1/003_sam/$name.sorted.bam
done

echo "Extracting base variation..."

mkdir $1/004_rc

for f in $1/003_sam/*.sorted.bam
do
	tempname=${f##*/}
	name=`basename $tempname .sorted.bam`
	echo "* $name"
	/usr/local/bin/bam-readcount -w 0 -f /Users/andrew/Lab\ data/Darja\ sequencing/bowtie-1.1.2/genomes/HIVss.fna $1/003_sam/$name.sorted.bam 2>/dev/null > $1/004_rc/$name.rc
done

mkdir $1/parse_results

echo "Calculating pause sites..."

./parse_sam.pl $1/003_sam/*.sam > $1/parse_results/regression.txt

echo "Writing base variation files..."
for f in $1/004_rc/*.rc
do
	tempname=${f##*/}
	name=`basename $tempname .rc`
	newname=${name%%_*}
	echo "* $name"
	./rc_extract.pl $1/004_rc/$name.rc > $1/parse_results/baseVariation_$newname.csv
done

echo "Process complete."
