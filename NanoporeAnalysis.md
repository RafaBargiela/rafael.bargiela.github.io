# Post-sequencing analysis of Nanopore reads
Description of the steps to follow in order to process Oxford Nanopore Technologies (ONT) reads after sequencing.

# Index
- [Dependencies](#Dependencies)
- [Concatenation](1.Concatenation-of-fastq-files)
- [Demultiplexing](#2.Demultiplexing)
- [Remomving adapters and filtering](#3.Adapters-removal-and-reads-filtering)
- [Classifying reads with Qiime2](#4.Analysis-of-filtered-reads-using-Qiime2)

## Dependencies
- PoreChop (https://github.com/rrwick/Porechop) or Deepbinner (https://github.com/rrwick/Deepbinner)
- PoreChop_ABI (https://github.com/bonsai-team/Porechop_ABI)
- Nanopack (https://github.com/wdecoster/nanopack), for using NanoFilt/Chopper
- Canu (https://github.com/marbl/canu)
- filtlong (https://github.com/rrwick/Filtlong)
- Qiime2 (https://qiime2.org/)
- Minimap2 (https://github.com/lh3/minimap2)
- Samtools (https://github.com/samtools/samtools)


## 1.Concatenation of fastq files
After sequencing, several or many fastq files will be produced from the same run. If sequencing process has already separated the different barcode files, you will have different folders corresponding to each barcode, with a group of fastq files inside. Otherwise, all files will be in the same folder.

Multiplexed or demultiplexed, all fastq files in the folder/s must be concatenated in one fastq file. Assuming we have files from 0 to 19, we will concatenate them following natural ordering:
```bash
cat File_{0..19}.fastq.gz > Concatenated.fastq.gz
```

## 2.Demultiplexing
Currently, demultiplexing process is already included on Nanopore sequencing process, as the tool _Gumpy_ is included in the pipeline. However, you can download this tool by yourself and use it if your run haven't been demultiplexed. Also, tools like Porechop or Deepbinner did the demultiplexing previously. Porechop is not longer maintained, but is still used since its efficiency and also for other tasks like adapter removal.

PoreChop can be used to demultiplex the reads according to barcodes, in case they haven't been demultiplexed before.
```bash
porechop -i input_reads.fastq.gz -b output_dir --threads "${CORES}"
```

## 3.Adapters removal, correction and reads filtering
On this part we first need to removed any possible adapter in the reads using _PoreChop_ABI	_. Then, reads will be fitered by length and quality with _filtlong_. Finally, selected reads will be corrected using _canu_.

### 3.1.Trimming adapters with PoreChop_ABI
PoreChop_ABI is an extension of old PoreChop which is still maintained. Its advantage over the old version is that it doesn't need external databases for the adapters. It combines the use of an internal adapters database with an algorithm which infers the adapters sequence on the reads. Additionally, you can add custom adapter to this database (you must search for the file _adapters.py_) or use a custom adapters text file which will be combined with the adapters database by using _-cap_ option.

In our case, the custom adapters file will include primers for 16s rRNA V4 region amplification, Nextera adapters, 18S rRNA, and amplification of capsid proteins for the analysis of viral sequences. The file will have the following format:
```
V416S
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGBCAGCMGCCGCGGTAA
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACHVGGGTWTCTAAT
18S
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAACCTGGTTGATCCTGCCAGT
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTGATCCTCCTGCAGGTTCACCTAC
Nextera
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGCCAGCMGCCGCGGTAA
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACHVGGGTWTCTAAT
```
First line is the name of the adapter, then the forward sequence and then the reverse. If there is only one sequence, leave one line emptly.

We will run PoreChop_ABI as follows:
```bash
porechop_abi -abi -i "${IN}".fastq.gz -o "${OUT}".trim.fastq.gz --format fastq.gz -cap CustomAdapters.txt --threads "${CORES}"
```

### 3.2.Filtering the reads
Once adapters have been removed, reads have to be filtered. 

#### 3.2.1 Reads Correction
Optionally, we can try to correct the reads using programs like the assembler _canu_:
```bash
canu -corrected -p "${OUT}" -d "${OUTDIR}" genomeSize="${SIZE}" -nanopore "${OUT}".trim.fasta.gz
```
Note that $SIZE must be changed according the origin of the reads to correct.
However, take in mind that its output is in fasta format, so you won't be able to use quality values on downstream analysis. Also, using only the corrected reads will shorten the number of reads to analysis.

#### 3.2.2 Reads filtering
Now adapters are gone, we will discard those sequences with length < 500 bp and/or mean quality < 10. We can do it by using __NanoFilt__ (currently being changed to _Chopper_) from the _nanopack_ python package.
```bash
gunzip -c "${IN}".fastq.gz | NanoFilt -q 10 -l 500 | gzip > "${IN}".Filtered.fastq.gz
```
Alternatively, we can also use the program __filtlong__:
```bash
filtlong --min_length 1000 --keep_percent 80 --min_mean_q 10 "${IN}".fastq.gz | gzip > "${OUT}".Filtered.fastq.gz
```
With this command line we will make sure that the program doesn't make a massive discarding of reads, keeping at least the 80% of the reads.

This program also allows to use reference sequence to do the trimming, either Illumina PE sequence or an assembly in fasta format:
```bash
filtlong -1 Illumina.R1.fastq -2 Illumina.R2.fastq --min_length 1000 "${IN}".fastq.gz | gzip > "${OUT}".Filtered.fastq.gz
filtlong -a Assembly.fasta --min_length 1000 "${IN}".fastq.gz | gzip > "${OUT}".Filtered.fastq.gz
```
In this case, PhreD quality would be ignored and it would be evaluated by _K_-mer matching to the reference sequences.

## 4.Analysis of filtered reads
Downstream analysis depends on the origin of our samples. If they are amplicons of some marker genes like 16S or 18S rRNA genes, we will try to classify them into OTUs or ASVs and then classify their taxonomy againts a reference database like Greengenes or SILVA.
Another possibility is that we have marker genes but from viral capsid, where we will need to get a database of this proteins and map the reads to the different capsid genes to finally get a summary of reads counts and get a viral classification.
Finally, the third alternative is having genomic or metagenomic reads, where we will need to assembly the sequences. On this point, we can also be interested on getting Metagenomic Assembled Genomes (MAGs), so maybe reads will have to be mapped over reference genomes to get the final independent MAGs.

### 4.1. Clusterization of reads into OTUs/ASVs and taxonomic classification
On this step we will use _DADA2_ and _QIIME2_ in case the ASVs option. There is also the alternative of clusterized the reads by OTUs using usearch or similar programs. Finally, clusters will be compare against last SILVA database for taxonomic classification.

### 4.2. Analysis of amplicons from Capsid VPs (viruses)

#### 4.2.1. Creating a custom database of Capsid VPs proteins.
First, we have to make our own custom Capsid VPs database. To do this we have downloaded all the metadata of proteins related to "Capsid VP" from the NIAID Virus Pathogen Database and Analysis Resource (ViPR), which is currently part of the Bacterial and Viral Bioinformatics Resource Center (BV-BRC, https://www.bv-brc.org/).
From the metadata, we extracted all the accession numbers and downloaded the aminoacid sequences corresponding to each protein sequence from the NCBI by using a custom Perl script and _E-Utilities_(https://www.ncbi.nlm.nih.gov/books/NBK25500/).
Since this database combine data from different gene databases (PATRIC,RefSeq, ...), some genes can be repeated, so we can filter the produced fasta file to ensure there are more than one sequence with the same accession number.
Finally, we obtained a custom database on fasta format with 97,665 sequences of Capsid VPs which we from now on will call _BVBRC.Capsid_VPs.fna_ file. Using _minimap2_ we are going to create and index of this sequences:
```bash
minimap2 -d BVBRC.Capsid_VPs.nr.mmi BVBRC.Capsid_VPs.nr.fna
```
With the _-d_ option we send the indexed sequences to _BVBRC.Capsid_VPs.nr.mmi_.

#### 4.2.2. Mapping reads over the Capsid VPs database and creating summary file
Using again _minimap2_ we will map the corrected reads produce by _canu_ against our custom database, resulting in a SAM alignment file.
```bash
minimap2 -ax map-ont BVBRC.Capsid_VPs.nr.mmi "${OUT}".trim.fasta.gz > "${OUT}".BVBRC.sam
```
Regard here in the _-x map-ont_ option, which is a specific setting for nanopore reads.
Now, we will create a BAM alignment from the SAM file and then a sorted BAM file.
```bash
samtools view -b "${OUT}".BVBRC.sam > "${OUT}".BVBRC.bam
samtools sort "${OUT}".BVBRC.bam -o "${OUT}".BVBRC.sorted.bam
```
Now, we will build a summary table of the proteins that have succesfully mapped with some reads and the count of them. To do this, first we will discard those unmapped reads. These reads are flagged with _4_ in the second column of the line, and the reference sequence name (3rd column) will be with _*_. We can clean that directly with the command:
```bash
samtools view -F 4 "${OUT}".BVBRC.sorted.bam > "${OUT}".BVBRC.sorted.filtered.txt
```
Now, we can patch the text file to to create the final summary table of counts.
