# Assembling Whole Genome Shotgun reads

## Dependencies

- BBtools 
- SPADES 

## 1. Quality control and filtering

Before running the assembling we will check the raw reads and make a filtering in order to keep only the best of them. For this we will make use of the different tools on _BBtools_ suite.
First, in order to re-pair some reads that could be disordered or had some mate eliminated we will run _repair.sh_:
```bash
repair.sh in="${FILENAME}.R1.fastq.gz" in2="${FILENAME}.R2.fastq.gz" out="${FILENAME}.R1.repaired.fastq.gz" out2="${FILENAME}.R2.repaired.fastq.gz" repair=t outs="${FILENAME}.singletons.fastq.gz"
```
After this, we ensure we have all pairs of reads sorted.

Then, we will get rid of contaminant sequences or artifacts and remove adapters. To do this, we will use the script _bbduk.sh_ in two differernt steps, using the output files from _repair.sh_:

- First, removing the artifacts or contaminant sequenes:
```bash
bbduk.sh in="${FILENAME}.R1.repaired.fastq.gz" in2="${FILENAME}.R2.repaired.fastq.gz" out="${FILENAME}.R1.noArtifacts.fq.gz" outm="${FILENAME}.R1.discarded.fq.gz" out2="${FILENAME}.R2.noArtifacts.fq.gz" outm2="${FILENAME}.R2.discarded.fq.gz" hdist=1 k=31 ftm=5 ref=artifacts threads=40
```
- Here, removing the adapters, using the output from the step above as input:
```bash
bbduk.sh in="${FILENAME}.R1.noArtifacts.fq.gz" in2="${FILENAME}.R2.noArtifacts.fq.gz" out="${FILENAME}.R1.noArtifacts.Trim.fq.gz" out2="${FILENAME}.21.noArtifacts.Trim.fq.gz" ktrim=r k=23 mink=11 hdist=1 tbo=t qtrim=r trimq=20 ref=adapters threads=40
```

Finally, pair-end reads were merged using bbmerge.sh:
```bash
bbmerge.sh in="${FILENAME}.R1.noArtifacts.Trim.fq.gz" in2="${FILENAME}.21.noArtifacts.Trim.fq.gz" out="${FILENAME}.merged.fq.gz" outu1="${FILENAME}.R1.unmerged.fq.gz" outu2="${FILENAME}.R2.unmerged.fq.gz" ihist=insert_size.txt usejni=t
```
Both, merged and unmerged files will be used as input on Spades.

## 2. Assembling metagenomic reads

Once filtering and trimming of reads are done, we will proceed with the assembling, for what we will use Spades with all output files produces by bbmerge, on -meta mode for metagenomic data:

```bash
spades.py --meta -s "${FILENAME}.merged.fastq.gz" -1 "${FILENAME}.R1.unmerged.fq.gz" -2 "${FILENAME}.R2.unmerged.fq.gz"  -t "${CORES}" -m "${MEM}" -o "${OUTPUT_DIR}" --checkpoints all
```
Here we can use multi-threading with -t setting the number of cores ($CORES) and the dedicated memory with $MEM on -m option. For large metagenomes is convenient to use checkpoints, so in case of crashing we can resume the assembling using the -continue or -restart-from options.

