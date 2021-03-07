#!/bin/bash

## get_consensus_sequence_Illumina.sh ##
# It's a simple script to qc Illumina data and generate virus consensus sequences 
# Written by Filipe Moreira

# Usage:$ bash get_consensus_sequence_Illumina.sh ~/path/to/lib/dir ~/path/to/reference.fasta ~/path/to/trimmomatic/adapters.fasta	min_coverage number_of_threads pattern_dir

LIBDIR=$1
REF=$2
ADAPTERS=$3
MINCOV=$4
THREADS=$5
PATTERN=$6

cd $LIBDIR

mkdir -p RAW_QC
mkdir -p FILTERED_QC
mkdir -p CONSENSUS

echo "Sample,Raw,Paired_filtered,Unpaired_filtered,Mapped,Efficiency,Coverage,Cov_10x,Cov_100x,Genome_span,Genome_span_raw" > stats_report.csv

# For each sample
for i in *$PATTERN*
do
	
	
	# Print sample directory name
	echo "$i"
	
	
	# Enter sample diretory
	cd $i
	
	
	# Decompressing raw data
	echo "Decompressing raw data..."
	gunzip *gz 
	
	
	# Get name of raw data 
	R1="$(ls *fastq | head -1)"
	R2="$(ls *fastq | tail -1)"
	
	
	# QC report for raw data 
	fastqc -t $THREADS *fastq
	mv *fastqc.zip ../RAW_QC
	mv *fastqc.html ../RAW_QC
	
	
	# Filter data with trimmomatic
	# Spec HEADCROP:30 for amplicon sequencing
	echo "Performing strict QC..."
	trimmomatic PE -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 HEADCROP:30 MINLEN:50
	
	
	# QC report for filtered data
	fastqc -t $THREADS trim*fastq
	mv *fastqc.zip ../FILTERED_QC
	mv *fastqc.html ../FILTERED_QC
	
	# Concatenated unpaired reads
	cat trim.u.$R1 trim.u.$R2 > trim.uni.$R1
	
	# Get name of qc data 
	U=trim.uni.$R1
	R1=trim.p.$R1
	R2=trim.p.$R2
	NAME=$(echo $R1 | sed -E 's/_.+//g' | sed 's/trim.p.//')
	
	
	
	# Map reads agains ref genome
	echo "Mapping..."
	bowtie2-build $REF reference
	bowtie2 --very-sensitive --no-unal -p $THREADS -x  reference -1 $R1 -2 $R2 -U $U -S $NAME.sam
	
	
	# Manipulate mapping file
	echo "Indexing and sorting bam file..."
	samtools faidx $REF 
	samtools view -S -b $NAME.sam > $NAME.bam
	samtools sort $NAME.bam -o $NAME.sorted.bam
	samtools index $NAME.sorted.bam 
	
	
	# Call variants 
	echo "Calling variants..."
	bcftools mpileup -E -Ou -f $REF $NAME.sorted.bam  | bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz
	bcftools index calls.vcf.gz
	bcftools norm -f $REF calls.vcf.gz -Ob -o calls.norm.bcf
	bcftools filter --IndelGap 5 calls.norm.bcf -Ob -o calls.norm.flt-indels.bcf
	
	
	# Get consensus sequences
	echo "Inferring consensus sequences."
	cat $REF | bcftools consensus calls.vcf.gz > $NAME.temp.consensus.fa
	
	bedtools genomecov -bga -ibam  $NAME.sorted.bam > table_cov.txt
	bedtools genomecov -d -ibam  $NAME.sorted.bam > table_cov_basewise.txt
	
	grep -w 0$ table_cov.txt > table_coverage_zero.txt
	awk  '$4 < '$MINCOV'' table_cov.txt > table_coverage_min-cov-$MINCOV.txt
	
	bedtools maskfasta -fi  $NAME.temp.consensus.fa -fo masked.$NAME.raw-consensus.fa -bed table_coverage_zero.txt
	bedtools maskfasta -fi  $NAME.temp.consensus.fa -fo masked.$NAME.$MINCOV.consensus.fa -bed table_coverage_min-cov-$MINCOV.txt
	
	
	# Copy mincov sequence to consensus dir
	cp masked.$NAME.$MINCOV.consensus.fa ../CONSENSUS/
	
	
	# Write stats report
	calc() { awk "BEGIN{print $*}"; }
	
	# genome extension covered
	NBASES=$(grep -v '>' masked.$NAME.$MINCOV.consensus.fa | tr -cd 'N' | wc -c)
	NBASESRAW=$(grep -v '>' masked.$NAME.raw-consensus.fa | tr -cd 'N' | wc -c)
	LEN=$(wc -l table_cov_basewise.txt | awk '{print $1}')
	GENCOV=$(calc 1-$NBASES/$LEN)
	GENCOVRAW=$(calc 1-$NBASESRAW/$LEN)
	
	# number of raw reads
	RAW=$(grep -cE "^\+$" *fastq | head -n 1 | sed -E 's/.+\://g')
	RAW=$(calc $RAW*2)
	
	# number of paired end filtered reads
	PAIRED=$(grep -cE "^\+$" trim.p*fastq | sed -E 's/.+\://g' | head -n 1)
	PAIRED=$(calc $PAIRED*2)
	
	# number of unpaired filtered reads
	UNPAIRED=$(grep -cE "^\+$" $U | sed -E 's/.+\://g')
	
	# number of nicely mapped reads
	MAPPED=$(samtools view -c -F 260 $NAME.sorted.bam)
	
	# frequency of reads actually used for consensus sequence inference
	USAGE=$(calc $MAPPED/$RAW)
	
	# mean and SD depth/coverage
	SUM=$(cat  table_cov_basewise.txt | awk '{sum+=$3; print sum}' | tail -n 1)
	COV=$(calc $SUM/$LEN)
	
	# bases 10x>
	COV10=$(awk  '$3 > 10' table_cov_basewise.txt | wc -l)
	
	# bases 100x>
	COV100=$(awk  '$3 > 100' table_cov_basewise.txt| wc -l)
	
	#Printing to report...
	echo "$NAME,$RAW,$PAIRED,$UNPAIRED,$MAPPED,$USAGE,$COV,$COV10,$COV100,$GENCOV,$GENCOVRAW" >> ../stats_report.csv
	
	echo ""
	
	# Go to the previous directory
	cd ../

done

cd CONSENSUS/
for i in *fa
do
	NEWNAME=$(echo $i | sed  's/masked.//' | sed  's/.consensus.fa//g')
	REFNAME=$(cat $REF | head -n 1 | sed 's/>//')
	cat $i | sed "s/$REFNAME/$NEWNAME/" > renamed.$i
done
cd ..



# Multiqc plot for raw data
multiqc RAW_QC/
mv multiqc_report.html multiqc_report_RAW.html
mv multiqc_data/ RAW_QC/
mv multiqc_report_RAW.html RAW_QC/

# Multiqc plot for filtered data
multiqc FILTERED_QC/
mv multiqc_report.html multiqc_report_FILTERED.html
mv multiqc_data/ FILTERED_QC/
mv multiqc_report_FILTERED.html FILTERED_QC/


echo "#######################"
echo "########The end########"
echo "#######################"


cd

exit

