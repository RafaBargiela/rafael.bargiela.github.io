# PANGENOME ANALYSIS

## Dependencies
- NCBI-datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/)
- CheckM (https://github.com/Ecogenomics/CheckM)
- CheckM-databases (https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz)
- GET_Homologues (https://github.com/eead-csic-compbio/get_homologues)
- TrimAL (https://github.com/inab/trimal)
- ClipKIT (https://github.com/JLSteenwyk/ClipKIT)
- seqlim (https://github.com/kyungtaekLIM/seqlim) or Geneious

## 1. Download a list of genomes of the taxonomic group to study
In order to create the pangenome of a target group, first we need to gather a dataset of genomes belonging to microorganisms classified on this group. This can be done using the accession number of the assemblies from NCBI. We could get all the genomes included in the specific group on GTDB, if it is already included on this database. Otherwise, we could carefully select a list of accession numbers from NCBI, regarding their taxonomy.

To download all the genomes at once, make a list on a text file with an accession per line. Then, we can use a loop on a script to download each of the genomes using NCBI dataset tool:
```bash
LIST=$1
while read line;
do
	echo "Downloading ${line} genome"
	datasets download genome accession "${line}" --exclude-rna --filename "${line}.zip"
done < $LIST
```

## 2. Check quality of genomes using CheckM
We need that the set of genomes selected has minimum completeness percentages and a maximum percentage of contamination retained. To check this we will use CheckM:
```
checkm lineage_wf -x fna -t $CORES "${DIR}" "${DIR}/checkM"
```
On "${DIR}/checkM" directory, on _storage_ folder, you can check the completeness and contamination percentages on the stats.tsv files. If all genomes files are store in the same input directory, checkM will process the analysis for all of them and common output files will be appended. 

## 3. Pangenome Analysis
To analyse the core genome or Pangenome shared among the set of genomes selected we will use the program _get_homologues_. First, we will analyse the protein coding genes shared among all genomes and cluster them in protein families, which will be classified into **core families** (shared among all genomes), **dispensable** (only present in some of strains) and **unique** (present in one or two genomes). Then, we will use those proteins present in all the genomes to develop a phylogenetic tree.

### 3.1 Gathering all the protein fasta files

In order to get the conserved orthologous proteins we classify the protein-coding genes in protein families, using _Get_homologues_. To use the program, we will need to have the predicted proteins from the coding genes of each genome in the same folder. Usually in fasta format with _.faa_ as extension. When you download the genomes with _datasets_ they usually bring this file among all the downloaded dataset. If not, you can use _Prodigal_ to predict the coding genes and get the amino acids fasta file.

The use of _Prodigal_ is of special **good practice before to start the analysis**, since you will have **unique genome IDs** as prefix for each gene/protein on each genome. This will avoid any issues regarding possible duplicated sequence names in the MSAs in further steps, when MSAs trimming (TrimAl rename automatically the duplicated names, but ClipKIT doesn't) or in the final concatenation step. You can use _Prodigal_ as follows:

```bash
prodigal -i $GENOME -a "${OUT}.faa" -f gff -m -p single -o "${OUT}.gff"
```
This will provide the _.faa_ fasta files with the amino acid sequence of the coding genes which we need from each genome to run _get_homologues.pl_.

Another consideration before running the protein family prediction it would be to rename all the sequences on the _.faa_ files to **include the corresponding genome ID in the headers**. This would make easier further processes when analyzing the Core family proteins (**See section 4**).

### 3.2 Protein family prediction with Get_homologues

There are different options of how to run the program, however the easiest is to locate all the amino acid fasta files _.faa_ on the same folder and run the _get_homologues_ as follows:

```bash
get_homologues.pl -d "${DIR}" -m local -M -A -P -t 0 -C 50 -S 50 -n 40
```
Where $DIR is the directory where .faa files are located. Among the options set, regard on the -C and -S options where we configure a 50/50 rule to select proteins with at least 50% of identity and 50% of alignment coverage to build the protein families. A and P options are to calculate average indentity and percentage of conserved proteins among genomes. With -M we set orthoMCL as the algorightm to build the clusters, default is by biderectional best-hit (BDBH) or you can also use -G to set COGS as algorithmn (COGtriangle). It could be interesting to add the **option -e to exclude inparalogues** (check next section).

Output will be all send to a drectory named $DIR_homologues. Here we will find a lot of files, most of them gunzipped files corresponding to the individual blast search of each genome and the files corresponding to each of their blast databases made by _makeblastdb_. Among them, there is a folder which name usually starts with the name of the genome/proteome used as reference and the options used when running get_homologues.pl (e.g. GENOME_f0_0taxa_algOMCL_e0_C50_S50_). By default the program uses as reference the genome with least sequences, but you can set any you like the -r option. Inside there are all the fasta files for each of the protein families, which will be important for further analysis. Also, there is another output file with the same name than the folder above but ended as _clusters_list_, where there is a small description of each of the protein families and the list of the genomes analysed including this family.

Other output files produced by the command line we set above are _Avg_identity.tab_ and _POCP.tab_ files. These are both matrices with the average percentage of identity and the percentage of conserved proteins between each genome to the others.

### 3.3 Classification and annotation of protein families
#### 3.3.1 Classification of the Protein Families
Parsing the _clusters_list_ file, we can get the number of conserved families (present in all of the genomes), dispensable (present in some of them) or unique (only present in 1 or 2 genomes). On this file, first line of each cluster description starts by _cluster_. On this line we can retrieve the cluster number, the name and also the size. This one means the number of proteins in the protein family, but one protein could be repeated on the same genome, so it doesn't necessary mean the number of different genomes of the dataset where the family is present. The genomes where the proteins of the cluster have been found are listed on lines below, starting by \\(:\\). You must regard if any of these genomes is repeated, then there are more than one copy of this protein family type on that genome:

```
cluster 1_ANIHLMIM_00001 size=15 taxa=14 file: 1_ANIHLMIM_00001.faa dnafile: void
: E-plasma_bin.faa
: GCA_002204705.1.faa
: GCA_002498845.1.faa
: GCA_002502705.1.faa
: GCA_002505185.1.faa
: GCA_011334705.1.faa
: GCF_000152265.2.faa
: GCF_000195915.1.faa
: GCF_001402945.1.faa
: GCF_002078355.1.faa
: GCF_003205235.1.faa
: GCF_900083515.1_C.divulgatum_S5.faa
: GCF_900083515.1_C.divulgatum_S5.faa
: GCF_900090055.1_C.divulgatum_PM4.faa
: GCF_900176435.1.faa
```

On the example above, there are two proteins of this cluster repeated on the genome GCF_900083515.1_C.divulgatum_S5.faa. By default, _get_homologues.pl_ include **inparalogues** in the analysis, so they will be included in the protein cluster. Knowing if there are more than one copy of an ortholog in a genome is interesting, but for our further analysis we will need to remove the copies in order to have alignments for each protein cluster with the same number of sequences for concatenation. If you want to directly produce clusters without inparalogues, you must run _get_homologues.pl_ with the option -e.

#### 3.3.2 Annotation for Protein Families
To an annotation of each of the protein clusters assigned, mainly those belonging to the core group family proteins, we will use the command _annotate_cluster.pl_ over each of protein family fasta file:

```bash
annotate_cluster.pl -f "${CLUSTER_FILE}" -P -D > "${CLUSTER_NAME}.ann.txt" 2>&1
```
With this command line we will align the family cluster to the closest Pfam domain. Option -o in the program only saves the alignment produced, but we are also interested on the specific name and annotation of the aligned Pfam domain. To do this, we also need to save the output and also the error messages of the program, there is why the use of \(2>&1\). 
Finally, we can parse each of these files to get there specific domain annotation into a table.

Additionally, it would be interesting to add **annotation for COGs, arCOGs or KEGG ontologies**. Any of the annotations (Pfam domains, COGs, arCOGs or KEGGs) will be used on steps below to create the family size table for the ancestral reconstruction and gene gain/loss estimation.

## 4. Multiple alignment (MSA) of Core Family Proteins and concatenation
In order to produce a phylogenetic tree, we need first to select which are the core families, make an individual MSA for each one, and finally concatenate all of them to produce a concatenated MSA which will be use to build the tree. 

Core families will be those present in all the genomes analysed, but here we nee to pay attention on something. If we don't switch off the option of using **inparalogs** with _Get_homologues_, some of the core families not only will be present in all genomes but also will present more than one sequence from the same genome. Therefore, some families could have higher number of sequences of expected. For those cases, we would need to manually **remove some sequences** on these families, leaving just one sequence per genome in the corresponding fasta file, before performing any MSA. In order to make this process faster, it would be advisable to rename the sequences on each genome _.faa_ file to include in the **header the accession or ID** of the corresponding genome. Hence, it will be easy to find which sequence in a core family fasta file belong to the same genome. 

### 4.1 MSA of each core protein family
We need to create a MSA for each core family protein. As is said above, make sure all fasta files have the same number of sequences. This number must be equal to the number of genomes analysed (CORE families means they are present in all of them) and it there should be only ONE sequence from each genome.

To make the alignments we will use the program mafft with the algorithm L-INS-i:
```bash
	mafft --maxiterate 1000 --thread 40 --localpair "${file}" > "${DIR}/${CLUSTER}.aln.fasta"
```
where $file is the _.faa_ file of each core protein family, $DIR the directory where storing the alignments and $CLUSTER the name or ID of the protein family to give name to the alignment file.

### 4.2 Filtering and trimming each core protein family MSA
Filtering and trimming a MSA is a very common process when making multiple alignments for phylogenetic studies. However, currently is been test that the performance of the filtering programs don't improve or even worsen the alignments on the final phylogenetic results (Tan _et al._ 2015, Steenwick _et al._ 2020). Taking that in mind, below is displayed the command line to use TrimAl, one of the best performing alignment trimming programs, and clipKit, a new trimming tool which takes into account this issue and try to solve it by keeping parsimony informative sites (i.e. sites with at least two characters that each occur at least twice).
```bash
	trimal -in "${DIR}/${CLUSTER}.aln.fasta" -out "${DIR}/${CLUSTER}.trimal.aln.fasta" -automated1
```
or
```bash
	clipkit "${DIR}/${CLUSTER}.aln.fasta" -m kpic-smart-gap -o "${DIR}/${CLUSTER}.clipKit.aln.fasta"
```

### 4.3 Concatenation of all MSA

Geneious program allows concatenation of MLAs. So, if you are a current user of this licenced program, you can just import the folder with the MSAs and then go to Tools and concatenate. Then, you just have to export to Document the alignment. But if you are not so lucky to have the licence for Geneious, you can algo use _seqlim_, which is available on GitHub:

```bash
seqlim -o $Concatenatd.alignment.fasta -infmt fasta -outfmt fasta cath $MSA_FOLDER
```

Regard that option _cath_ is for horizontal concatenation, which is our aim. Other options results in a piled alignment. Once we have the concatenated MSA, we can use it to build a phylogenetic tree.

Finally, we will need to change the names of the sequences on the concatenated MSA. Depending on the program each concatenated alignment will have different type of names. _Geneious_ just concatenate sequences names and seqlim uses the name of the first sequence introduce. For a practical matter, we will rename each sequence on the concatenated MSA for the corresponding genome name, ID or a simple identifier easy to recognize (e.g. GENOME1). To do this, we have to take in mind the order used for the analysis, returned by _Get_homologues_ on the file **_input_order.txt_**. Order shown on this file is the same order followed to introduces sequences on the family protein fasta files and the same shown in the MSA. So, the first sequence on the concatenated MSA corresponds to the first genome on the list shown by input_order.txt. 

### 4.4 Building the phylogenetic tree

To get the phylogenetic tree we will use the R packages _ape_ and _phangorn_. The specific code will be **explained in another tutorial**, but there is one important consideration to have in mind: If you have not yet, rename the genomes on either the concatenated MSA or the final phylogenetic tree for the simpliest genome names as possible. **Without any non-alphanumerical character**. An use these short names on the next step when building the Family size table. Hence, you will avoid any strange results when using _Count_.

Also, it is important to take in mind the **tree labels for nodes**. If we calculate the bootstrap for the tree using _phangorn_, node labels will be changed by the bootstrap values. This will create issues when using _Count_, since it uses the tree node labels to name columns corresponding to ancestors. Therefore, if we don't change the node labels we won't be able to distinguish among these columns. **But, we also need the bootstrap values**, so the best option will be to add the node label a unique text but keeping the bootstrap value (e.g. from "100" to "Node1BS100"). Then, when plotting the tree with all the data we can parse the node labels again and keep only the bootstrap values.

## 5. Ancestral Reconstruction and gene gain/loss estimation

For this step we will use the program _Count_. To use it we need to starting files: a phylogenetic tree and a Family Size table. As tree, we will use the one we have built using the concatenated MSA with Core Famliy proteins. As is said in the previous step, genomes names on the Family Size table and on the tree must match completely and can't have any non-alphanumerical character.

Family Size table is a tab-separated documuent which looks like follows:
```
family Genome1 Genome2 Genome3 Genome4 Genome5 Genome6
arCOG00001 1 2 1 0 0 0
arCOG00002 1 0 1 1 1 1
arCOG00004 0 0 0 0 1 1
```
If there are some families with annotation of the type you want to use (arCOGs in this case), you can try to add the **Pfam annotations** or any other annotation that you have additionally done.

_Count_ on the command line has three different packages: _AsymmetricWagner_, _ML_ and _Posteriors_ and it runs under _Java_. The program will be run by the following structure:
```bash
 java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.$COUNT_PACKAGE -OPTIONS $tree $FamilySizeTable [$RatesTable]
```
Where _-Xmx4096M_ means the CPU memory to use in the process (4Gb in this case, in Megas), \\(ca.umontreal.iro.evolution.genecontent.\\) is the full class name of the packages in java. Tree and Family size table are always needed, but Rates table is only required by _Posteriors_.
 
### 5.1. Gene gain/loss estimation by Wagner Parsimony

First we will apply the package _AsymeticWagner_ to inffer the gain/loss rate over our protein families. Command line will be like this:
```bash
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner \ 
tree.newick FamilySizeTable.txt > WagnerParsimony.Output.txt
```
The output consists of three parts: family sizes at ancestral taxa (lines sarting with “# FAMILY”), genome sizes (lines starting with “# PRESENT”), and lineagespecific gene family size changes (lines starting with “# CHANGE”). A “# FAMILY” line lists the family sizes at each taxon that minimize the parsimony penalty of the reconstruction. Further columns give the number of lineageswhere the family was lost (“Losses”) or newly appeared (“Gains”), or where it expanded (“Expansions”) and reduced (“Reductions”) in size (quoting Count manual, Csurös M. 2010).
```
# FAMILY name Genome1 Genome2 Genome3 ... root Gains Losses Expansions Reductions
# FAMILY arCOG00001 0 0 0 ... 1 0 1 1 0
# FAMILY arCOG00002 1 1 1 ... 2 0 0 0 3
```
Lines starting with “# PRESENT” aggregate the same information across different families: for every taxon, a line gives the number of families with positive size,as well as the total of the family sizes (i.e., number of all genes). Lines starting with “# CHANGE” give aggregate information on lineage-specific changes (quoting Count manual, Csurös M. 2010). Regard that, after the genomes columns are those for the nodes, representing the estimations for the different branches of the tree, including the root or Lowest Common Ancestor (LCA).

### 5.2. Creating Rates table by Maximum Likehood

_Count_ uses the birth and death phylogenetic model to compute the ancestral reconstruction. To do that, we can estimate first the optimized parameters for its calculation using the _ML_ package. An example of calculation would be:

```bash
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.ML \ -opt_rounds 100 -duplication_k 4  -length_k 4 tree.newick FamilySizeTable.txt > Rates.txt
```
It is a good idea to perform the optimization in model hierarchy. Starting with a simple model and use to perform more complex rate models with it (as is suggested in _Count_ manual):

```bash
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.ML \ -uniform_duplication true tree.newick FamilySizeTable.txt > Rates1.txt
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.ML \ tree.newick FamilySizeTable.txt Rates1.txt > Rates2.txt
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.ML \ -max_paralogs 100 -length_k 3 tree.newick FamilySizeTable.txt Rates2.txt > Rates3.txt
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.ML \ -max_paralogs 100 -length_k 3 -duplication_k 3 ex.tre ex.txt Rates3.txt > Rates4.txt
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.ML \ -max_paralogs 10000 -length_k 4 -duplication_k 4 tree.newick FamilySizeTable.txt Rates4.txt > Final.Rates.txt
```

### 5.2. Infering ancestral gene content

The _Posteriors_ application infers ancestral gene content by posterior **probabilities** in a phylogenetic birth-and-death model. You can run it like this:

```bash
java -Xmx4096M -cp "${COMPLETE_PATH}"/Count.jar ca.umontreal.iro.evolution.genecontent.Posteriors \ tree.newick FamilySizeTable.txt Final.Rates.txt > Posteriors.output.txt
```
A **line** from _Posteriors_ output looks like this:
```
Family ... Censy:1 Censy:m Censy:gain Censy:loss Censy:expansion Censy:reduction ...
...
arCOG00001 ... 0.0 0.0 0.0 0.0 0.0 0.9985936233867079
```
The program also reports the posterior probabilities also for the empty phylogenetic profile in a line where family name is ABSENT.
Basically, instead of gain/loss absolute figures as WagnerParsimony application, Posteriors returns the gain/loss, expanssion and reduction probabilities of each family on each taxon and node. Below is the description of this probabilities from the manual, for a better understanding:

- **Family size**. The program computes two probability values: whether a family has/had 1 (p1) or multiple (p>1) members at a given (ancestral or terminal) taxon. The absence probability (family size of 0) can be computed as (1−p1−p>1). For taxon x, the posterior probabilities p1 and p>1 are listed under column headers x:1 and x:m, respectively. (Where x is the ID for the genomes used on the input tree/table). In other words, x:1 and x:m mean **Present families** and **Family with multicopy**, respectfully, on each genome. In case of ancestors, they mean **Probabilities** of being present and in multicopy.
- **Family gain and loss**. The program computes the probabilities that in a given lineage, the family size changed from a positive number to 0 (loss), or from 0 to a positive number (gain). The loss and gain probabilities for the edge leading to taxon x are given under the column headers x:loss and x:gain, respectively.
- **Expansions and reductions**. The program also computes the probabilities for size changes in retained families along each lineage. Expansions (size 45 change from 1 to something larger) and reductions (size change from 2 or more to 1) on the edge leading to taxon x are listed under the headers x:expansion and x:reduction.
- **Rate categories**. If the model has non-constant family-specific rate variations, then the output of Posteriors includes the posterior probabilities for families belonging into discrete categories. The columns for these probabilities have headers in the syntax of Cc/p, where c is the category’s machinegenerated identifier (a positive integer), and p is a point in the lattice defined by the Cartesian product of the discrete category indices. For example,C43/e2,d3,l0,t0, is rate class 43, in wich the the edge length (e), duplication rate (d), loss rate (l) and gain rate (t) category indices are 2,3,0 and 0,respectively.

### 5.3 Understanding the Ancestral reconstruction table

Once all the steps with Count are finished, the main output table for our attention is the table results from _Posteriors_. Here there are three groups of columns: First, columns for rate categories; second, columns for the genomes analysed, corresponding to the terminal taxons on the tree; and finally, the columns with probabilities for ancestors or nodes in the tree. We can avoid get into the rate categories, since we have information enough on the other columns.

Then, for each genome/terminal taxon and ancestor/node we have 6 different columns:
- Size columns: Named x/1 and x/m, which mean, as explained above, posibilities for each family to be present and in multicopy, respectfully. Regard that, in case of the terminal taxons (analysed genomes), possibilities here will be always 1 or 0, since we already know if the families are present and/or in multicopy. These are only probabilities in case of ancestors.
- Gain, loss, expansion and reduction: Here are the posibilities of each genome/ancestor on each family for gene gain, loss, expansion or reduction.

To analyse these probabilities, we can set a threshold of trust, e.g. 0.5, which use to filter the probabilities. Therefore, those posibilities over the threshold we will take then as TRUE and those lower as FALSE. For each genome/ancestor, their total number of TRUE in the gain, loss, expanssion, reduction columns would be taken as theorethical or putative gains, losses, expanssions or reductions. In case of ancestor, this will be also applied to the size columns. So the total number of TRUE families on each ancestor on the column x:1 or Present will be the putative number of present families on that ancestor/node; the same for the x:m or Multicopy column.

## 6. Final remarks

This has been a long way until complete the full analysis. Now, we need still to analysis the obtained results and produce nice figures. To do that, the most important outputs generated here are the **phylogenetic tree produced from the concatenated MSA** and the **Ancestral Reconstruction table** produced at the end of the analysis with _Count_, after _Posteriors_.

Using Perl/R scripts you could turn the Ancestral reconstruction table to a TRUE/FALSE table according with the threshold of probabilities set. Then, with other scripts you can count the TRUEs on each column for ancestors and genomes and make a summary file.

Probably, it would be a good idea to add the family annotation to the Ancestral Reconstruction table, to ease the interpretation of the data. Also, making a draft phylogenetic tree to see the location of each of the nodes will help us to check which are the important ancestors (nodes) for our analysis.

When 