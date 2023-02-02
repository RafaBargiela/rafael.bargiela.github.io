---
title: "Sequence Variants searching using reference genome"
author: "Rafael Bargiela"
date: "18-February-2022"
output: html_document
---
## Tools needed
- bwa or Bowtie2
- samtools
- iVar or Freyja

## Indexing referece genome
```{bash}
bwa index $REF2
samtools faidx $REF2
```

## Aligning Pair-end fastq sample file to reference
Using Bowtie2
```{bash}
bowtie2 -x "${REF}" -1 "${DIR}/${line[0]}.${line[1]}.R1.001.fastq.gz" -2 "${DIR}/${line[0]}.${line[1]}.R2.001.fastq.gz" -S "SAM/${line[0]}.sam"
```

Using BWA
```{bash}
bwa mem -o "ALIGNMENTS/${line[0]}.sam" "${REF2}" "${DIR}/${line[0]}.${line[1]}.R1.001.fastq.gz" "${DIR}/${line[0]}.${line[1]}.R2.001.fastq.gz"
```
Converting sam file to sorted bam fil
```{bash}
samtools view -b -o "ALIGNMENTS/${line[0]}.bam" "ALIGNMENTS/${line[0]}.sam"
samtools sort -o "ALIGNMENTS/${line[0]}.sorted.bam" "ALIGNMENTS/${line[0]}.bam"
```
## Merging alingments when we have more than one replicate
```{bash}
   samtools view -b -o "ALIGNMENTS/${line[0]}.bam" "ALIGNMENTS/${line[0]}.sam"
	 samtools sort -o "ALIGNMENTS/${line[0]}.sorted.bam" "ALIGNMENTS/${line[0]}.bam"
```

## Variants calling
Using samtools and iVar
```{bash}
samtools mpileup -aa -A -d 0 -B -Q 0 -o "VARIANTCALL/${SAMPLENAME}.consensus.merged.samtools.mpileup.tsv" --reference "${REF2}" "ALIGNMENTS/${SAMPLENAME}.merged.bam"
ivar variants -p "${SAMPLENAME}.IvarVariants" -q 20 -t 0.03 -r "$REF2" -g "${REF}.gff" < "VARIANTCALL/${SAMPLENAME}.consensus.merged.samtools.mpileup.tsv"
```
Using Freyja (recommended for SARS-CoV-2 variants seach)
```{bash}
freyja variants "ALIGNMENTS/${SAMPLENAME}.merged.bam" --variants "${SAMPLENAME}.FreyjaVariants" --depths "${SAMPLENAME}.FreyjaDepths" --ref "$REF2" 
freyja demix "${SAMPLENAME}.FreyjaVariants.tsv" "${SAMPLENAME}.FreyjaDepths" --output "${SAMPLENAME}.FreyjaDemix"
```

## Summarizing output from Freyja
```{bash}
Freyja.Demix2TxT.pl <DEMIX_FILE> > <OUTPUT.FreyjaDemix.txt>
```
```{bash}
If we have more than one replicate, to make combined tables with results:
```
```{bash}
Freyja.CombineTxT.pl <DIR_with_Demix.txt_files> > <Output>
```

