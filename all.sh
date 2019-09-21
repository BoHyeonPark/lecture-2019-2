#!/bin/bash

if [ $# -ne 1 ];then
  echo "#usage: sh $0 [samplename]"
  exit
fi

SAMPLE=$1
BWA="/data/etc/bwa/bwa"
SAMTOOLS="/data/etc/samtools/bin/samtools"
REFERENCE="/data/reference/ucsc.hg19.fasta"
JAVA="/usr/bin/java"
PICARD="/data/etc/picard/picard.jar"
GATK="/data/etc/gatk/GenomeAnalysisTK.jar"
SNPEFF="/data/etc/snpEff/snpEff.jar"
MILLS="/data/etc/bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
A1KG="/data/etc/bundle/1000G_phase1.indels.hg19.sites.vcf"
DBSNP138="/data/etc/bundle/dbsnp_138.hg19.vcf"

##1. Map to Reference
#1.1 BWA mem: FASTQ to SAM
#${BWA} mem -R "@RG\tID:test\tSM:${SAMPLE}\tPL:ILLUMINA" ${REFERENCE} ${SAMPLE}_1.filt.fastq.gz ${SAMPLE}_2.filt.fastq.gz > ${SAMPLE}.mapped.sam

#1.2 Samtools: SAM to BAM
#samtools view –Sb ${SAMPLE}.mapped.sam > ${SAMPLE}.mapped.bam

#1.3 Samtools sort: Make Sorted BAM
#samtools sort –o ${SAMPLE}.sorted.bam ${SAMPLE}.mapped.bam
#samtools view ${SAMPLE}.mapped.bam | head
#samtools view ${SAMPLE}.sorted.bam | head

#1.4 (Opional) Run from 1.1 to 1.3 in one command
#${BWA} mem -R "@RG\tID:test\tSM:${SAMPLE}\tPL:ILLUMINA" ${REFERENCE} ${SAMPLE}_1.filt.fastq.gz ${SAMPLE}_2.filt.fastq.gz | ${SAMTOOLS} view -Sb - | ${SAMTOOLS} sort - > ${SAMPLE}.sorted.bam


##2. Mark Duplicate
#2.1 Picard MarkDuplicate: Sorted BAM to Markdup BAM
#${JAVA} -jar ${PICARD} MarkDuplicates I=${SAMPLE}.sorted.bam O=${SAMPLE}.markdup.bam M=${SAMPLE}.markdup.metrics.txt

#2.2 Samtools index: Make BAM index
#${SAMTOOLS} index ${SAMPLE}.markdup.bam


##3. GATK
#3.1 GATK Target Realign
#${JAVA} -jar ${GATK} -T RealignerTargetCreator -R ${REFERENCE} -I ${SAMPLE}.markdup.bam -known ${MILLS} -known ${A1KG} -o ${SAMPLE}.intervals #-L bed
#${JAVA} -jar ${GATK} -T IndelRealigner -R ${REFERENCE} -I ${SAMPLE}.markdup.bam -known ${MILLS} -known ${A1KG} -targetIntervals ${SAMPLE}.intervals -o ${SAMPLE}.realign.bam #-L bed

#3.2 GATK BaseRecalibrator
#${JAVA} -jar ${GATK} -T BaseRecalibrator -R ${REFERENCE} -I ${SAMPLE}.realign.bam -knownSites ${MILLS} -knownSites ${A1KG} -knownSites ${DBSNP138} -o ${SAMPLE}.table #-L bed
#${JAVA} -jar ${GATK} -T PrintReads -R ${REFERENCE} -I ${SAMPLE}.realign.bam -o ${SAMPLE}.recal.bam -BQSR ${SAMPLE}.table #-L bed

#3.3 GATK HaplotypeCaller, GenotypeGVCFs
${JAVA} -jar ${GATK} -T HaplotypeCaller -R ${REFERENCE} -I ${SAMPLE}.recal.bam --emitRefConfidence GVCF --dbsnp ${DBSNP138} -o ${SAMPLE}.g.vcf #-L bed
#{JAVA} -jar ${GATK} -T GenotypeGVCFs -R ${REFERENCE} -V ${SAMPLE}.g.vcf -o ${SAMPLE}.raw.vcf #-L bed

#3.4 GATK Variant Filter
# Select SNP
#{JAVA} -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${SAMPLE}.raw.vcf -o ${SAMPLE}.raw.snp.vcf --selectTypeToInclude SNP

# Filter SNP
#{JAVA} -jar ${GATK} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLE}.raw.snp.vcf -o ${SAMPLE}.filtered.snp.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 ||  MQRankSum < -12.5  ||  ReadPosRankSum < -8.0" --filterName SNP_FILTER

# Select INDEL
#{JAVA} -jar ${GATK} -T SelectVariants -R ${REFERENCE} -V ${SAMPLE}.raw.vcf -o ${SAMPLE}.raw.indel.vcf --selectTypeToInclude INDEL

# Filter INDEL
#{JAVA} -jar ${GATK} -T VariantFiltration -R ${REFERENCE} -V ${SAMPLE}.raw.indel.vcf -o ${SAMPLE}.filtered.indel.vcf --filterExpression "QD < 2.0  || FS > 200.0 || ReadPosRankSum < -20.0" --filterName INDEL_FILTER

# Combine SNPs and INDELs
#{JAVA} -jar ${GATK} -T CombineVariants -R ${REFERENCE} -V ${SAMPLE}.filtered.snp.vcf -V ${SAMPLE}.filtered.indel.vcf -o ${SAMPLE}.filtered.vcf -genotypeMergeOptions UNSORTED


##4. SnpEff: Annotation
#4.1 Annotate Variants
#{JAVA} -jar -Xmx4g ${SNPEFF} -v hg19 ${SAMPLE}.filtered.vcf > ${SAMPLE}.filtered.annotated.vcf
