# Virus prospection on Parys Mountain metagenome
## 1.Assembling Metagenome from RAW fastq files
Before to search viruses on the metagenome we need to create an assembly for the raw sequencing reads. Two approaches will be follow. First, an assembly will be performed by MegaHit with previous reads curation using bbtools. Second approach will be using metaSPAdes.
### 1.1 Assembling using Megahit
First, we order the pair-end reads, in case files could be desordered:
```bash
repair.sh in="${R1}" in2="${R2}" out="${OUT}/${NAME}.R1.repaired.fastq.gz" out2="${OUT}/${NAME}.R2.repaired.fastq.gz" repair=t outs="${OUT}/${NAME}.singletons.fastq.gz" 
```
Then, we use bbdul to remove possible artifacts and adapters sequences:
```bash
bbduk.sh in="${OUT}/${NAME}.R1.repaired.fastq.gz" in2="${OUT}/${NAME}.R2.repaired.fastq.gz" out="${OUT}/${NAME}.R1.artifacts.fastq.gz" outm="${OUT}/${NAME}.R1.artifacts.discard.fastq.gz" out2="${OUT}/${NAME}.R2.artifacts.fastq.gz" outm2="${OUT}/${NAME}.R2.artifacts.discard.fastq.gz" hdist=1 k=31 ftm=5 ref=artifacts threads="${CORES}"
bbduk.sh in="${OUT}/${NAME}.R1.artifacts.fastq.gz" in2="${OUT}/${NAME}.R2.artifacts.fastq.gz" out="${OUT}/${NAME}.R1.clean.fastq.gz" out2="${OUT}/${NAME}.R2.clean.fastq.gz"  ktrim=r k=23 mink=11 hdist=1 tbo=t qtrim=r trimq=20 ref=adapters threads="${CORES}"
```
Once raw reads are curated, paired files are merged:
```bash
bbmerge.sh in="${OUT}/${NAME}.R1.clean.fastq.gz" in2="${OUT}/${NAME}.R2.clean.fastq.gz" out="${OUT}/${NAME}.merged.fastq.gz" outu1="${OUT}/${NAME}.R1.unmerged.fastq.gz" outu2="${OUT}/${NAME}.R2.unmerged.fq.gz" ihist=insert_size.txt usejni=t
```
Finally, assembly is performed using MegaHit with metagenomic preset and contigs of 1,000 bp of minimum length:
```bash
megahit --read "${OUT}/${NAME}.merged.fastq.gz" -o "${OUT}/MEGAHIT" --out-prefix "${NAME}" -t "${CORES}" --presets meta-sensitive --min-contig-len 1000
```
This methods produced an assembly of 69,271,856 bp in 29,787 contigs.

### 1.2 Assembling using metaSPAdes
For this method we start from the last step previously to run MEGAHIT on the last section, using outputs from bbmerge. To run SPADES for metagenomic assembling we used the following command set up:
```bash
spades.py --meta -s "${OUT}/${NAME}.merged.fastq.gz" -1 "${OUT}/${NAME}.R1.unmerged.fq.gz" -2 "${OUT}/${NAME}.R2.unmerged.fq.gz"  -t "${CORES}" -o "${OUT}/SPADES"
```
On this case, the assembly produced is much longer, with 244,233,985 bp on 160,964 scaffolds.

## 2.Search of viruses using virSorter
To search virus genomes on our metagenomic assembly we will use virSorter, comparing results from two available version (1 and 2).

### 2.1 Using virSorter (version 1)
Database was previously downloaded and located on virsorter-data folder.
```bash
wrapper_phage_contigs_sorter_iPlant.pl -f $FAS --db 1 --wdir $OUT --ncpu 40 --data-dir $VIRDATA --no_c
```
### 2.2 Using virSorter2
First, as same as virSorter 1, we need to download the database, what is made after the installation with the command below:
```bash
virsorter setup -d db -j 4
```
where -d specifies the path and name of the folder for the database.
It seems that virSorter2 provides a more complete database and its output is more ellaborated than v1. So, for now on we will reffer to virSorter2 output.

VirSorter provides a set of viral sequences detected on your metagenome, but not a list of taxonomic assignations. In order to do this, we have to use the sequences detected by virSorter and assign them a taxonomy by futher pipeline.

### 2.3 Taxonomic assignation of viral sequences detected with virSorter
On this step we will use the output from VirSorter2. Firstly, we need to extract the putative genes/proteins from the viral scaffolds provided by virSorter. For this prediction we will use *Prodigal* with the following commands:
```bash
prodigal -i $FAS -a "${OUT}.faa" -f gff -m -p meta -o "${OUT}.gff"
```
Note that we specified *-p meta* option to perform the prediction over a metagenome, and the option *-m* to mas the N's on the sequences and avoid the prediction of genes across them.
## 3. Direct genome virus search using MetaviralSPAdes (ONGOING)

```bash
spades.py --metaviral -s "${MERGED}" -1 "${UNMERGED1}" -2 "${UNMERGED2}" -t 40 -o "${OUT}_SPADESmetaviral"
```
 
## 4. Classification of viral sequences using VCONTACT2
Using previously produced predictions for viral sequences from metagenomes, in order to get a taxonomic classfication of the viral diversity, we need first to predict genes/proteins using Prodigal:
```bash
prodigal -i $FAS -a "${OUT}.faa" -f gff -m -p meta -o "${OUT}.gff"
```
Once we have all the predicted genes/proteins, to get the classification we will use the program VconTact2. First, we will use the command gene2genome to format the Prodigal output to vcontact format:
```bash
vcontact2_gene2genome -p "${OUT}.faa" -o "${OUT}.g2g.csv" --source-type Prodigal-FAA
```
Then, the *.g2g.csv* files produced is used by vcontact, combined with protein/genes fasta file, to finally classify the viral content.
```bash
vcontact2 --raw-proteins "${OUT}.faa" --rel-mode Diamond --proteins-fp "${OUT}.g2g.csv" --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /scratch/b.chsa18/CONDA/conda_envs/VirSorter2/bin/cluster_one-1.0.jar --output-dir "${OUT}_VCONTACT"
```
