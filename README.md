# Building a Peptide Atlas on Liverpool HPC

[**Introduction**](#introduction)

[**Making a fasta search database**](#making-a-fasta-search-database)

[**Choosing and retrieving peptide data**](choosing-and-retrieving-peptide-data)

[**Processing peptide data**](#processing-peptide-data)

<p style="margin-left: 40px">
<ol type="a">
  <li>
    <a href="https://github.com/VibrantStarling/make-protein-atlas/edit/main/README.md#3a-search-tools-">Search tools</a>
  </li>
 <li>
   <a href="https://github.com/VibrantStarling/make-protein-atlas/edit/main/README.md#3b-instrument-parameters-">Instrument parameters</a>
   </li>
 <li>
   <a href="https://github.com/VibrantStarling/make-protein-atlas/edit/main/README.md#3c-qc-of-data-">QC of data</a></li>
</ol>
</p>

[**Combining and collapsing data**](combining-and-collapsing-data)

[**Visualising peptides on a genome**](visualising-peptides-on-a-genome)

[**Running the pipeline**](running-the-pipeline)

# Introduction <a name="introduction"></a>

For this task, we want to get peptide-level support for genes within a genome, sometimes called proteogenomics. We source mass spectrometry (MS) data from ProteomeXchange, and its underlying repositories \- PRIDE, MassIVE, IProx etc. We then define a search database (see below), and run the Trans Proteomic Pipeline on our linux HPC, to search MS raw data against the database, generate statistics etc. The data is then passed to Eric Deutsch at PeptideAtlas, who loads it into their database, and does a few extra steps of post-processing. We then get data out in TSV format. 

**Several steps of this guide have been converted into a set of scripts and a pipeline (process-pxd.sh) that runs them all together.** The github for all these files can be found here: [https://github.com/VibrantStarling/make-protein-atlas](https://github.com/VibrantStarling/make-protein-atlas). You can also find further instruction on how to use this repository at the bottom of this document, [here](#running-the-pipeline).

But first, a step by step explanation of what we’re doing…

You will need to ensure the following are installed:
  - pridepy
  - pandas
  - numpy
  - matplotlib
  - gffread
  - miniconda

```
# DO THIS IN YOUR HOME DIRECTORY

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
source ~/miniconda3/bin/activate
conda init --all
```

```
# activate python and install the other programs

module load python   
pip3 install pandas matplotlib numpy pridepy
conda install bioconda::gffread
```

# 1\. Making a fasta search database <a name="making-a-fasta-search-database"></a>

The fasta database is a compilation of all the relevant proteins that you want to be able to find in your peptide data, including common contaminants and reverse reads. The makeup of a database will depend on what you are looking for, but generally you want a balance of size and specificity. A large database with an excess of irrelevant reads increases computation time without improving statistical power. Conversely, a small database with highly specific reads could be very efficient, but give misleading results. 

There is no hard and fast rule for what to include in a database. Try out a few variations, or talk to your supervisor and colleagues as to what seem reasonable

## Starting out

Basic principles \- when performing proteogenomics, we generally start from an assumption that the “official” gene models are not completely correct. We want to find peptides that map onto the proteins from those gene models (the “official” proteome), but also find peptides that could  map to other gene models. These gene models may have been predicted by pipelines in the past \- like “Ensembl *ab initio*” models, by running *de novo* predictions e.g. with BRAKER3 / FunAnnotate (instructions on how to do this: [https://github.com/VibrantStarling/genomics-from-the-ground-up](https://github.com/VibrantStarling/genomics-from-the-ground-up) ), or by assembling transcriptomics data into plausible transcripts. For example, if long read (ISO-SEQ) transcripts exist, these can be translated in three frames, taking all open reading frames that run from an ATG to a stop codon. 

It is generally not advisable to perform a six-frame genome translation. This produces a very large search database (losing statistical power), and does not predict intron-spanning peptides (ISPs). ISPs are some of the most useful datatypes we can get from proteomics.

## Making your own proteome fasta

There are several ways to compile fasta sequences for a search database. The end aim is to have a range of protein sequences in fasta format that represent genes you want to find within the mass spec protein data.

### 1\. Establish a base from:

**Preexisting proteome data**

1. Find your reference proteome on Uniprot, VEuPathDB, NCBI or wherever has the best annotation for your species, e.g. [https://www.uniprot.org/uniprotkb?query=proteome:UP000002311](https://www.uniprot.org/uniprotkb?query=proteome:UP000002311)   
2. Download the proteome in fasta format, and make sure it is amino acid sequence

**Amnio acid fasta files**

1. Select your fasta files with your desired genes  
2. If you have multiple files, merge all protein fasta files into one, e.g. using cat \*.fasta \> ../final\_database/full\_species\_merged.fasta. *It does not matter if you have repeat sequences.*

**gff/gtf files**

1. Retrieve the gff or gtf of the sets of proteins you are interested in searching (or the amino acid fasta file with an associated gff). You might have one or several.  
2. Download the nucleotide genome fasta for the gff/gtf  
3. Use something like gffread or AGAT to extract the *protein* amino acid sequence from the gff.   
   1. `gffread -y CDS_aa.faa -g genome.fa input.gff `
   2. `agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t cds -p`
4. Then merge all protein fasta files into one, e.g. using `cat *.fasta > ../final_database/full_species_merged.fasta`

### 2\. Add contaminants and decoys 

With [FragPipe](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html#specify-a-protein-sequence-database) - this quickly and simply adds both contaminants and decoys in one and produces everything as one fasta. Make a note of your decoy prefix. Fragpipe can be downloaded and run as a GUI program as depicted below, or in commandline.

**You can pass a fasta file (combined from the official proteins and your extra annotations) to fragpipe like this:**  
![](https://github.com/VibrantStarling/make-protein-atlas/blob/main/Building%20a%20Peptide%20Atlas%20on%20Liverpool%20HPC/images/image3.png)

# 2\. Choosing and retrieving peptide data <a name="2.-choosing-and-retrieving-peptide-data"></a>

When retrieving your datasets, make sure you target your species of interest. For pathogens in particular, you want to ensure data is NOT from host infected proteomes, although it can be co-cultured datasets if this is common for a given species. For these cases, you should include the host proteome (e.g. human) in the search database.

For general global proteomics, you want to avoid datasets enriched for post translational modifications (PTMs) like phosphorylation, glycosylation or acetylation etc. You also want to avoid interactome and TurboID sets which usually only identify a few proteins. Ideally, source articles for datasets would mention that they have identified \>2K \- 3K proteins, or a minimum say of \>1K. 

## Choosing datasets

1. Go to [https://proteomecentral.proteomexchange.org/](https://proteomecentral.proteomexchange.org/) and search for your species.   
2. Download the results as a TSV, then move in to google sheets to annotate them.  
3. Add a column for:  
   1. **Priority**  \- for scoring how likely you are to use them (0,1,2 \-\> No, Maybe, Yes)  
   2. **Notes** \- keep track of reasons for your scoring, particularly Maybes  
   3. **Parameters** \- track what label type was used (TMT, iTRAQ, label-free, *etc*.). This will determine what parameters we will use for MS-Fragger searches  
   4. **Instrument** \- track what instrument was used for the dataset   
      1. Favour Thermo instruments. Orbitrap Fusion, Orbitrap Lumos, Orbitrap Astral, Q-Exactive, Q-Exactive HF are best. Maybe include LTQ Orbitrap if paper reports it is a big data sets e.g. identifying \>2000 proteins). LTQ Orbitrap have low-res fragment ions, and need this setting change in search parameters  
      2. Exclude Waters (SYNAPT) instruments, and usually exclude Bruker (timsTOF) and SCIEX (TripleTOF), unless you are struggling for data.  
   5. **Data dependent acquisition (DDA)** \- If paper or submission mentions Data independent acquisition (DIA), we exclude.   
      1. Markers of DIA include  (aka. SWATH) include using software DIA-NN, MaxDIA, Spectronaut or Skyline.   
      2. Markers of DDA including using MaxQuant, ProteomeDiscoverer, FragPipe…  
   6. **Source database.** It is easier to get fast download from PRIDE direct to the Linux cluster. If large datasets are in other databases, we can get them via a more manual process e.g. using FileZilla for FTP transfer to a local PC, then copying them to the cluster, but this is slow and annoying\!  
4. Go through the sheet and find the datasets you want to use.   
5. Often the link to a source publication is missing. If you google search for the PXD identifier, you can often find this embedded in a paper. Keep track of this (add manually to the publication column). We can inform PRIDE support to add the publication link to the dataset.

## Downloading datasets

**NOTE \-** this step is included in the pipeline at the end of this document and can handle multiple PXDs at once. This section breaks down what you would need to do to do it yourself for one dataset.

1. On the login/head node on HPC make a directory for your species in  /mnt/hc-storage/groups/peptideatlas\_builds/. Avoid spaces in your directory name, call it something informative like Toxoplasma\_gondii\_ME49, ToxoplasmaGondiiME49, or TgondiiME49  
2. activate [screen]() so you can leave your window to idle without risk of it timing out
```
screen
```
3. make sure pridepy is installed and use it to download your pride accession
```
pridepy download-all-public-raw-files -a PXD###### -o PXD###### -p aspera
```
4. Check that your downloads have worked

# 

# 3\. Processing peptide data <a name="processing-peptide-data"><a/>

The code for these steps can be found in this slurm script: TPP-fragger-wholedir-postprocess-GENERAL.sh

## 3a\. Search tools <a name="search-tools"></a>

A typical pipeline has these steps, all run with SLURM:

1. Convert raw to mzML with msconvert  
2. Run search with Comet or MS-Fragger \- for peptide atlas tasks Andy favours MS-fragger for higher speed  
   1. **Absolute key step is setting the variable modifications and instrument tolerances correctly.**  
   2. There are pros and cons to using array jobs here, Andy prefers a single slurm script that doesn’t use array jobs but locks 80 threads (one machine) per PXD from start to finish (if cluster is quiet, this is perfect, if busy, it is less ideal \- potentially hogging resources in some steps). A single job is default in the pipeline script too.  
3. Run TPP tools \- PeptideProphet (xinteract) and IProphet (InterProphetParser)  
4. Convert pepXML to tsv and do hacky post-processing to check everything has worked (Andy’s python scripts).   
   1. **I would like to re-engineer these final steps to produce a clean “build” for VEuPathDB, including converting peptides to BED or BAM format for genome display, and sorting peptides that are only present in some annotation sets but not others. This is on the TO DO list for someone in the team.**

## 3b\. Instrument parameters <a name="instrument-parameters"><a/> 

MS-fragger requires a file that contains specific parameters for each type of data (label-free, TMT, TMTpro, etc.). In some cases this might be changing the mass of cysteine or lysine to account for different labeling. A summary of MS-fragger parameters can be found [here](https://github-wiki-see.page/m/Nesvilab/MSFragger/wiki/Setting-the-Parameters).   
The most important settings are getting the variable and fixed modifications correct. If these are incorrect, we will not identify many peptides. Assuming we are not searching for enriched PTMs, then we have to focus on learning the quantification reagents (if any) used.

We have several .param files with set parameters for the most common label types we come across.

The main thing you need to remember to change across all .param files  are the num\_threads and database\_name options. This is done automatically in the pipelines, but is something you should bear in mind.

**Other very important parts of the param file:**
```
variable\_mod\_01 .. 16 Sets variable modifications (variable\_mod\_01 to variable\_mod\_16). Space separated values with 1st value being the modification mass and the second being the residues (specified consecutively as a string) it modifies. You can find information on the syntax [here](https://github-wiki-see.page/m/Nesvilab/MSFragger/wiki/Setting-the-Parameters). Here are some examples from one of the .param files:

variable\_mod\_01 \= 15.994900 M 3 (methionine oxidation)  
variable\_mod\_02 \= 42.010600 \[^ 1 (protein N-terminal acetylation)  
\# variable\_mod\_03 \= 79.966330 STY 3 (for phosphorylation, commented out with ‘\#’ because its not needed)  
\# variable\_mod\_04 \= \-17.026500 nQnC 1 (for pyro-Glu or loss of ammonia at peptide N-terminal, commented out because its not needed)  
\# variable\_mod\_05 \= \-18.010600 nE 1 (for glutamic acid at the N-terminal,  commented out because its not needed)  
variable\_mod\_06 \= 229.162930 n^ 1 (TMT on peptide N-terminus)  
variable\_mod\_07 \= 229.162930 S 1 (TMT on serine))
```

Fixed modifications set static numbers for mass for various peptides

```
add\_Q\_glutamine \= 0.000000  
add\_K\_lysine \= 229.162930  
add\_E\_glutamic\_acid \= 0.000000
```

For TMTpro, fixed modifications all become `304.207146`
For label-free, comment out variable mods and put fixed lysine mass to `0.0000`

Another critical parameter is instrument resolution. Assuming a recent Orbitrap / Q-Exactive, we run in high-high resolution mode (precursors and fragments), which looks like this:
```
precursor\_mass\_lower \= \-10        	\# Lower bound of the precursor mass window.  
precursor\_mass\_upper \= 10        	\# Upper bound of the precursor mass window.  
precursor\_mass\_units \= 1        	\# Precursor mass tolerance units (0 for Da, 1 for ppm).  
data\_type \= 0                    	\# Data type (0 for DDA, 1 for DIA, 2 for gas-phase fractionation DIA).  
precursor\_true\_tolerance \= 10    	\# True precursor mass tolerance (window is \+/- this value).  
precursor\_true\_units \= 1        	\# True precursor mass tolerance units (0 for Da, 1 for ppm).  
fragment\_mass\_tolerance \= 10    	\# Fragment mass tolerance (window is \+/- this value).  
fragment\_mass\_units \= 1           	\# Fragment mass tolerance units (0 for Da, 1 for ppm).  
calibrate\_mass \= 2                	\# Perform mass calibration (0 for OFF, 1 for ON, 2 for ON and find optimal parameters).  
use\_all\_mods\_in\_first\_search \= 0	\# Use all variable modifications in the first search (0 for No, 1 for Yes).  
decoy\_prefix \= rev\_             	\# Prefix of the decoy protein entries. Used for parameter optimization only.
```

If we run any LTQ-Orbitrap (these are high-low), so we change fragment mass tolerance to 0.5Da. The rest of the parameters remain the same.
```
fragment\_mass\_tolerance \= 0.5    	\# Fragment mass tolerance (window is \+/- this value).  
fragment\_mass\_units \= 0           	\# Fragment mass tolerance units (0 for Da, 1 for ppm).
```

## 3c\. QC of data <a name="qc-of-data"><a/>

Andy has made some [python QC scripts](https://github.com/VibrantStarling/make-protein-atlas/tree/main/python-qc-scripts) that filter data and convert them into tsv format. In order to run these, you need to logged into the login node have run:

module load python   
pip3 install pandas matplotlib numpy

You should only have to do this once, but sometimes updates to the login node’s Python might reset this. Consider setting up a conda environment.

After running the python scripts, this extracts the peptide-spectrum match (PSM) data from the final pepXML file in the pipeline, into various CSV files, and applying FDR thresholding \- the final file is: PXD0XXXXX\_interact-ipro.pep\_thresholded.tsv. A stats file is also produced PXD0XXXXX\_psm\_stats.tsv, which produces data like this:

| raw\_file | spectra\_count | psm\_count | recovery |
| :---- | :---- | :---- | :---- |
| 20201215_CDPK1TimeCourse_TMTPro1\_WP\_DMSO\_fxn1.raw | 39076 | 13058 | 0.334169311 |
| 20201215_CDPK1TimeCourse_TMTPro1\_WP\_DMSO\_fxn3.raw | 42884 | 13547 | 0.315898703 |
| 20201215_CDPK1TimeCourse_TMTPro1\_WP\_DMSO\_fxn4.raw | 45142 | 13986 | 0.309822338 |
| 20201211_CDPK1TimeCourse_TMTPro1\_WP\_ZAP\_fxn1.raw | 30469 | 9240 | 0.30325905 |
| 20201211_CDPK1TimeCourse_TMTPro1\_WP\_ZAP\_fxn4.raw | 42892 | 12905 | 0.300871957 |

psm\_count is count of PSMs at 1% FDR,   
spectra\_count is count of total spectra,   
recovery is the 1% FDR count / total count

If the analysis has gone well, we are looking at recovery of say 25-60% PSMs per raw file (and ideally a large volume of spectra searched). Sometimes the proportion drops really low (say \<10%) \- this is usually because we have put the search settings wrong, or perhaps it is low quality data. 

# 

# 4\. Combining and collapsing data <a name="combining-and-collapsing-data"></a>

Run [run-combine-and-threshold-BATCH.sh](https://github.com/VibrantStarling/make-protein-atlas/blob/main/slurm-and-shell-scripts/run-combine-and-threshold-BATCH.sh) to combine, collapse, and threshold all the individual PXD files by peptide.

sbatch run-combine-and-threshold-BATCH.sh pxdlist.txt outsuffix

# 5\. Visualising peptides on a genome <a name="visualising-peptides-on-a-genome"></a>

*ANDY: Ideally, I would like to convert peptide-level data to BED or BAM format for visualising directly in genome browsers, like JBrowse. This is technically fairly straightforward, since we could get the GFF of gene models from which proteins are derived. From proteomics results, we get the position of the peptide within a protein, and then need an extra step to figure out the corresponding genomic position \- including mapping some peptides across intron junctions. The algorithm is a little fiddly but not really that complex, but I don’t currently have a code base I am happy with. This is a task for someone to take on…*

This has yet to be fully hashed out. In theory the steps might include:

1. Find a way to translate pepXML to a gff.  
2. Map peptides to the genome and generate a SAM file  
3. Use samtools to convert the SAM file to BAM, then index the BAM file  
4. Use bedtools intersect to identify shared regions between the gff and bam file ([http://quinlanlab.org/tutorials/bedtools.html\#bedtools-intersect](http://quinlanlab.org/tutorials/bedtools.html#bedtools-intersect))  
5. Create a bigwig file and use that to visualise continuous protein data [https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/10\_data\_visualization.html](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/10_data_visualization.html) 

You can load a BED or BAM file in to IGV to visualise the data alongside genome data

proBAMsuite, a possibly useful tool? [https://www.mcponline.org/article/S1535-9476(20)33669-0/pdf](https://www.mcponline.org/article/S1535-9476\(20\)33669-0/pdf) 

# 

# Running the pipeline <a name="running-the-pipeline"></a>

**process-pxd.sh combines all the necessary steps to download, process, and qc a set of PXD files. You still have to assemble the PXD files and fasta database yourself, but the rest should be automated.**

It allows you to run each stage of the pipeline separately, so if you need to add a file, or just rerun the postprocessing step, you can do it without reprocessing every single file.

You can download a repository with all the necessary directories and directory structures from 

git clone https://github.com/VibrantStarling/make-protein-atlas.git

**General usage:**

If you do not specify any of the flags (-D, \-M, or \-F), then it will run the whole pipeline by default.
```
bash process-pxd.sh -l pxd-list [-b database.fasta] [-q python-qc-directory] [-s slurm-script-directory] [-p fragger-params-directory] [-D] [-M] [-F]
```

***This script is functional, but not optimised, so you may need a specific file layout for it (as arranged in this reposiotory). However, in theory, you could have one folder for all the scripts and just supply that one file to each of the directory options (-q, \-s, \-p)***


**You will need:**

1. **A PXD list.** A list of PXD accessions and the type of data in them. Accession and type are delimited by a tab, and each row has a separate entry. Type defines which MS-fragger params file will be used. It will always look for the pattern of `fragger\_${type}.params`, where `${type}` could be filled in with `LM` to make `fragger\_LM.params`. If you make a new .params file, make sure it follows that naming convention.  
```
PXD015269  TMT
PXD033642  TMT
PXD040598  TMTpro
PXD047027  LF
```
2. **A database fasta file.** A fasta file of protein information to search against. This is a collection of your organism of interest, reverse reads, and contaminants.  
3. A directory that contains all the [**python QC scripts**](https://github.com/VibrantStarling/make-protein-atlas/tree/main/python-qc-scripts)
4. A directory containing [**additional slurm and bash scripts.**](https://github.com/VibrantStarling/make-protein-atlas/tree/main/slurm-and-shell-scripts)
5. A directory containing [**MS fragger parameter files,**](https://github.com/VibrantStarling/make-protein-atlas/tree/main/ms-fragger_example_params) with names that correspond to the type specified in the PXD-list 

Feed the two files and three directories to the process-pxd.sh script and it will run the whole pipeline by default. It will download and produce directories in your current directory.

bash process-pxd.sh \-l pxd\_list.txt \-b path/to/database.fasta \-q path/to/python-qc-directory \-s path/to/slurm-and-bash-script-directory \-p path/to/MS-fragger-params-directory

If you need to run the download (-D), mzml conversion (-M), or fragger search and postprocessing steps (-F) separately, add the relevant flag. 

**After you have run the pipeline go to step 4, combining and collapsing data**
