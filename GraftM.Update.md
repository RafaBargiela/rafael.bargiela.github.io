# Using GraftM for taxonomic classification of metagenomes
Steps to get a taxonomic classification of a Whole Genome Sequencing (WGS) from metagenomic samples. Classification is based on marker genes (e.g., 50S Ribosomal protein L2), so databases provided by GraftM can be also updated with additional sequences in case to need to enrich specific groups.

## Get GraftM and GraftM database
Instruction for installation of the program are available on:
https://github.com/geronimp/graftM
Basically, you can install it as a python3 package:
```bash
pip install graftm
```
Databases for the analysis are available on:
https://data.ace.uq.edu.au/public/graftm/
We will use the database based on 50S ribosomal protein L2: https://data.ace.uq.edu.au/public/graftm/7/7.07.ribosomal_protein_L2_rplB.gpkg.tar.gz

## Dependencies
orfm v. >= 0.2.0 (https://github.com/wwood/OrfM)
hmmer v. >= 3.1b1 (http://hmmer.janelia.org/)
mfqe v. >= 0.5.0 (https://github.com/wwood/mfqe)
pplacer v. >= 2.6.32 (http://matsen.fhcrc.org/pplacer/)
krona v. >= 2.4 (http://sourceforge.net/p/krona/home/krona/)
mafft v. >= 7.22 (http://mafft.cbrc.jp/)
diamond v. >= 0.9 (https://github.com/bbuchfink/diamond) 
FastTreeMP (http://www.microbesonline.org/fasttree/)

Best way to install all dependencies is to create a conda environment and call it before running GraftM, otherwise maybe you could have some troubles with the installation of some of them:
```bash
git clone https://github.com/geronimp/graftM
cd graftM
conda env create -n graftM -f graftm.yml
conda activate graftM
```
This will also download the full program. If you haven't installed via pip as explained above, you can installed from source as follows from the code above:
```bash
cd bin
export PATH=$PWD:$PATH
graftM -h
```
However, I recommend to use the installation from pip.

## 1. Updating the packages provided by GraftM
We don't need necessarily to update its database, but in case we have some specific sequences from microorganisms which we would like to trace, it would be a good idea to produce a new database including these additional sequences.
 
### 1.1. Get a list of new sequences to update
Make a list of genomes to look for markers, the easiest way is to filter metadata from GTDB (https://data.gtdb.ecogenomic.org/releases
/release207/207.0/ar53_metadata_r207.tar.gz for archaea), GTDB identifiers are prefixed with GB_ or RS_ (GenBank and RefSeq respectively), the prefixes should be removed and the assemblies may be downloaded from NCBI via Entrez or batch download.

Assemblies can also be dowloaded using NCBI datasets program (https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/), using the following command line:
```bash
      	datasets download genome accession Assembly_accession --exclude-rna --filename output.zip
```
You can create an script to do this in a loop.

### 1.2. Searching for the marker protein
Predict proteins (i used prodigal, many assemblies in NCBI lack predicted CDS). 
Search for the markers in predicted proteomes with HMM provided in the GraftM package (i set relatively strict e-value around 1e-10 and took the top hit, the hits may be checked by Pfam search).

### 1.3. Creating taxonomy file
#### 1.3.1 Taxonomy table for GraftM analysis
Make a table of new 9-digit ID, genome assembly and taxonomy string. Rename the proteins with new IDs according to the table and put them to a single file like "seqs.fasta". This IDs are arbitrary.

Make a tab-separated table with IDs and corresponding taxonomy string (example files are here https://github.com/geronimp/graftM/tree/main/example_data/create) an name it e.g. "taxonomy.tsv".

#### 1.3.2 Updating the taxonmy table to GTDB
Since the most accepted taxonomic names are those assigned by the GTDB database, you can update the taxonomies on your table to their GTDB names, when possible. To do this, we can use the function _TabLinToGTDB_, from the in-house Perl module GTDB.pm, in a Perl script processing the table line by line. Remember, lineages must be introduces as tab-separated.

#### 1.4. Checking conflicts with old database
Taxonomy.tsv shoul be checked for conflicts with existing taxonomy in the updating package and one of the files should be corrected if such conflicts exists. All spaces in species names should be replaced to underlines.

### 1.5. Running GraftM to update the package
Finally, GraftM package is updated running the following command:
```bash
Run graftM update --graftm_package old.gpkg --taxonomy taxonomy.tsv --sequences seqs.fasta --output new.gpkg
```
You will see information on any conflicts during package update and return to point 4 to remove them. Regard that new package is not a file but a folder, so to reference the package for analysis you will need to provide the path to this folder.

## 2. Running taxonomic classification with GraftM
Once the package for classification is ready, we can run GraftM _graft_ command. In our case, we will perform the classification for a pair-end sequenced metagenome, previously assembled, using the assembly as reference:
```bash
graftM graft --forward R1.fastq --reverse R2.fastq --expand_search_contigs Ref.assembly.fasta --graftm_package new.gpkg --evalue 1e-10 --placements_cutoff 0.85 --threads "${CORES}" --output_directory "${OUTDIR}"
```




