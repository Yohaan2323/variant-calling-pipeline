FULL CODE

1. Download FASTQ files

Wget
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR101/070/SRR101470/SRR101470_1.fastq.gz

wget
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR101/070/SRR101470/SRR101470_2.fastq.gz

2. Index the reference genome (chr1.fa)

bwa index chr1.fa

3. Align reads to the reference using BWA-MEM

bwa mem chr1.fa SRR101470_1.fastq.gz SRR101470_2.fastq.gz >
SRR101470_aligned.sam

4. Convert SAM to BAM (keep SAM)

samtools view -Sb SRR101470_aligned.sam > SRR101470_aligned.bam

5. Sort and index the BAM for IGV

samtools sort SRR101470_aligned.bam -o SRR101470_sorted.bam

samtools index SRR101470_sorted.bam

6. Mark duplicates with Picard

java -jar picard.jar MarkDuplicates \

I=SRR101470_sorted.bam \

O=SRR101470_marked.bam \

M=marked_dup_metrics.txt

7. Index marked BAM file

samtools index SRR101470_marked.bam

8. Add read groups (required for GATK tools)

/home/yohaan/practice/ref_genome/chr1/gatk-4.5.0.0/gatk
AddOrReplaceReadGroups \

-I SRR101470_marked.bam \

-O SRR101470_rg.bam \

--RGID SRR101470 \

--RGLB Lib1 \

--RGPL ILLUMINA \

--RGPU SRR101470 \

--RGSM SRR101470

9. Index BAM with read groups

samtools index SRR101470_rg.bam

10. Prepare known-sites files for BQSR

Download and extract only chromosome 1:

# Download

wget
https://rcs.bu.edu/examples/bioinformatics/gatk/ref/Homo_sapiens_assembly38.dbsnp138.vcf

wget
https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz

# Compress and index if needed

bgzip Homo_sapiens_assembly38.dbsnp138.vcf

bcftools index Homo_sapiens_assembly38.dbsnp138.vcf.gz

# Extract chr1 and index

bcftools view -r chr1 Homo_sapiens_assembly38.dbsnp138.vcf.gz -Oz -o
chr1_dbsnp.vcf.gz

bcftools index chr1_dbsnp.vcf.gz

bcftools view -r chr1 Homo_sapiens_assembly38.known_indels.vcf.gz -Oz -o
chr1_known_indels.vcf.gz

bcftools index chr1_known_indels.vcf.gz

11. Create FASTA dictionary (required by GATK)

gatk CreateSequenceDictionary -R chr1.fa -O chr1.dict

12. Perform Base Quality Score Recalibration (BQSR)

gatk BaseRecalibrator \

-R chr1.fa \

-I SRR101470_rg.bam \

--known-sites chr1_dbsnp.vcf.gz \

--known-sites chr1_known_indels.vcf.gz \

-O recal_data.table

13. Apply BQSR

gatk ApplyBQSR \

-R chr1.fa \

-I SRR101470_rg.bam \

--bqsr-recal-file recal_data.table \

-O SRR101470_bqsr.bam

14. Index recalibrated BAM

samtools index SRR101470_bqsr.bam

15. Call variants using HaplotypeCaller

gatk HaplotypeCaller \

-R chr1.fa \

-I SRR101470_bqsr.bam \

-O SRR101470_raw_variants.vcf.gz
