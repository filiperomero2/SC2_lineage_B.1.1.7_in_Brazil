#!/bin/bash

## GetConSeqVir.sh ##

# It is a simple pipeline for inference of viral consensus genome sequences from Illumina data.
# The pipeline performs quality control, mapping, variant calling and consensus sequence inference using a standard bioinformatics stack.

# It requires 3 positional arguments:
#	1 - Absolute path for libraries root directory;
#	2 - Path for a reference genome in fasta format;
#	3 - Path for a fasta file harboring Illumina adapter sequences.

# In addition, minimum sequencing coverage and number of threads may also be specified:
#	4 - Minimum sequencing coverage to call a base on the consensus sequence (default = 100x)
#	5 - Number of threads available for processing (default = 1)

# Example usage:$ ./GetConSeqVir.sh ~/LIBRARIES/RUN_1/ ~/REFERENCE_GENOMES/reference.fasta ~/ASSETS/adapters.fasta 200 6

##############################################################################################################################

# Set variables...
LIBDIR=$1
REF=$2
ADAPTERS=$3
MINCOV=$4
THREADS=$5


##############################################################################################################################
# Define functions and validate arguments...

checkfile() { if [[ -f $1 ]] ; then echo "$1 ->  file identified.";  else echo "- File $1 not found!"; exit; fi; }
checkdir() { if [[ -d $1 ]] ; then echo "$1 -> Directory identified";  else echo "- Directory $1 not found!"; exit; fi; }
calc() { awk "BEGIN{print $*}"; }

checkdir $LIBDIR
checkfile $REF
checkfile $ADAPTERS


if [ -z "$MINCOV" ]
then
      echo "Minimum coverage is not defined, using the default (100x)"
	  MINCOV=100
else
      echo "Minimum coverage threshold identified: $MINCOV"
fi

if [ -z "$THREADS" ]
then
      echo "Number of threads not specified, starting with the default (1)"
	  THREADS=1
else
      echo "Number of threads available for processing: $THREADS"
fi


##############################################################################################################################
# Start the pipeline!

FILE=$(echo $REF | sed 's/.*\///g')
REFNAME=$(cat $REF | head -n 1 | sed 's/>//')
REFPATH=$(echo $REF | sed "s/$FILE//g")
cd $REFPATH
bowtie2-build $REF reference
REF2=$(pwd)/reference

cd $LIBDIR
mkdir -p RESULTS/RAW_QC/ RESULTS/FILTERED_QC/ RESULTS/CONSENSUS/
echo "sample,number_of_raw_reads,number_of_paired_filtered_reads,number_of_unpaired_filtered_reads,number_of_mapped_reads,efficiency,average_depth,coverage_10x,coverage_100x,coverage_1000x,genome_coverage" > stats_report.csv

# For each sample
for i in */
do
	if [ $i == "RESULTS/" ]
	then
		continue
	fi

	# Print sample directory name
	echo "$i"
		
	# Enter sample diretory
	cd $i
	
	# Decompress raw data
	echo "Decompressing raw data..."
	gunzip *gz 
		
	# Get the names of raw data files
	R1="$(ls *fastq | head -1)"
	R2="$(ls *fastq | tail -1)"
		
	# QC report for raw data 
	fastqc -t $THREADS *fastq
	mv *fastqc.zip ../RESULTS/RAW_QC
	mv *fastqc.html ../RESULTS/RAW_QC
		
	# Filter data with trimmomatic
	echo "Performing strict QC..."
	trimmomatic PE -phred33 $R1 $R2 trim.p.$R1 trim.u.$R1 trim.p.$R2 trim.u.$R2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 HEADCROP:30 MINLEN:50
	
	# QC report for filtered data
	fastqc -t $THREADS trim*fastq
	mv *fastqc.zip ../RESULTS/FILTERED_QC
	mv *fastqc.html ../RESULTS/FILTERED_QC
	
	# Concatenate unpaired reads
	cat trim.u.$R1 trim.u.$R2 > trim.uni.$R1
	
	# Get the names of qc reads
	U=trim.uni.$R1
	R1=trim.p.$R1
	R2=trim.p.$R2
	NAME=$(echo $R1 | sed -E 's/_.+//g' | sed 's/trim.p.//')
	
	# Map reads against reference genome and handle the bam file
	echo "Mapping..."
	bowtie2 --very-sensitive --no-unal -p $THREADS -x $REF2 -1 $R1 -2 $R2 -U $U | samtools view -bS - > $NAME.bam   
    samtools faidx $REF 
	samtools sort $NAME.bam -o $NAME.sorted.bam
	samtools index $NAME.sorted.bam 
	
	# Call variants 
	echo "Calling variants..."
	bcftools mpileup --max-depth 10000 -E -Ou -f $REF $NAME.sorted.bam  | bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz
	bcftools norm -f $REF calls.vcf.gz -Oz -o calls.norm.vcf.gz
	bcftools filter --IndelGap 5 calls.norm.vcf.gz -Oz -o calls.norm.flt-indels.vcf.gz
    bcftools index calls.norm.flt-indels.vcf.gz
	
	
	# Get consensus sequences
	echo "Inferring consensus sequences..."
	cat $REF | bcftools consensus calls.norm.flt-indels.vcf.gz > $NAME.temp.consensus.fa
	bedtools genomecov -bga -ibam  $NAME.sorted.bam > table_cov.txt
	bedtools genomecov -d -ibam  $NAME.sorted.bam > table_cov_basewise.txt
	awk  '$4 < '$MINCOV'' table_cov.txt > table_coverage_min-cov-$MINCOV.txt
	bedtools maskfasta -fi  $NAME.temp.consensus.fa -fo masked.$NAME.consensus.fa -bed table_coverage_min-cov-$MINCOV.txt
	sed -i "s/$REFNAME/$NAME/"  masked.$NAME.consensus.fa
    
    # Copy mincov sequence to consensus dir
	cp masked.$NAME.consensus.fa ../RESULTS/CONSENSUS/

	# Write to the stats report
	# genome extension covered
	NBASES=$(grep -v '>' masked.$NAME.consensus.fa | tr -cd 'N' | wc -c)
	LEN=$(wc -l table_cov_basewise.txt | awk '{print $1}')
	GENCOV=$(calc 1-$NBASES/$LEN)
	
	# number of raw reads
	RAW=$(grep -cE "^\+$" *fastq | head -n 1 | sed -E 's/.+\://g')
	RAW=$(calc $RAW*2)
	
	# number of paired end filtered reads
	PAIRED=$(grep -cE "^\+$" trim.p*fastq | sed -E 's/.+\://g' | head -n 1)
	PAIRED=$(calc $PAIRED*2)
	
	# number of unpaired filtered reads
	UNPAIRED=$(grep -cE "^\+$" $U | sed -E 's/.+\://g')
	
	# number of properly mapped reads
	MAPPED=$(samtools view -c -F 260 $NAME.sorted.bam)
	
	# frequency of reads actually used for consensus sequence inference
	USAGE=$(calc $MAPPED/$RAW)
	
	# mean depth
	SUM=$(cat  table_cov_basewise.txt | awk '{sum+=$3; print sum}' | tail -n 1)
	COV=$(calc $SUM/$LEN)
	
	# bases 10x>
	COV10=$(awk  '$3 > 10' table_cov_basewise.txt | wc -l)
	
	# bases 100x>
	COV100=$(awk  '$3 > 100' table_cov_basewise.txt | wc -l)

    # bases 1000x>
	COV1000=$(awk  '$3 > 1000' table_cov_basewise.txt | wc -l)
	
	#Printing to report...
	echo "$NAME,$RAW,$PAIRED,$UNPAIRED,$MAPPED,$USAGE,$COV,$COV10,$COV100,$COV1000,$GENCOV" >> ../stats_report.csv
	
	echo ""
	
	# Go to the previous directory
	cd ../

done

cd RESULTS/

cat CONSENSUS/*fa > CONSENSUS/CONSENSUS.fasta

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

cd ..


echo "#######################"
echo "####### The end #######"
echo "#######################"


exit
