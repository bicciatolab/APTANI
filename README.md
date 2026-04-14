# APTANI and APTANI<sup>2</sup>

**APTANI** and **APTANI<sup>2</sup>** are two methods for selecting potentially relevant aptamers through a sequence-structure analysis. Both tools contain modules to predict specific secondary structures in each selection round and to rank aptamers by motifs embedded in their predicted structures.

#### Contact:

silvio.bicciato@unipd.it; mattia.forcato.1@unipd.it

#### Citations:

- **APTANI**: Caroli J, Taccioli C, De La Fuente A, Serafini P, Bicciato S. APTANI: a computational tool to select aptamers through sequence-structure motif analysis of HT-SELEX data. _Bioinformatics_. 2016 Jan 15;32(2):161-4  doi: [10.1093/bioinformatics/btv545](https://doi.org/10.1093/bioinformatics/btv545)
- **APTANI<sup>2</sup>**: Caroli J, Forcato M, Bicciato S. APTANI2: update of aptamer selection through sequence-structure analysis. _Bioinformatics_ 2020 Apr 1;36(7):2266-2268.  doi: [10.1093/bioinformatics/btz897](https://doi.org/10.1093/bioinformatics/btz897)

# Table of Contents

- [APTANI](https://github.com/bicciatolab/APTANI#APTANI)
- [APTANI<sup>2</sup>](https://github.com/bicciatolab/APTANI#APTANI2)
- [APTANI and APTANI<sup>2</sup> sample data, scripts and output](https://github.com/bicciatolab/APTANI#Example)

## APTANI

**APTANI** consists of four major steps. The **first step** calculates the frequencies of the aptamer sequences produced by the SELEX process. Given the input file, APTANI calculates the frequency of each individual aptamer sequence by dividing the number of occurrences of each sequence by the total number of analyzed sequences. Frequency counts can be quantified either using the whole aptamer sequence or only its variable region, e.g., to count only those sequences that change during the SELEX process. Setting a threshold on the frequency counts, the pool of original aptamers can be filtered to select only high frequent aptamers for the further steps. In the **second step**, APTANI predicts the secondary structure of each aptamer that passes the frequency filter, and extracts all motifs represented in these structures. As secondary structures, we consider four loop structures: hairpins, bulge (either right or left) and intra-strands. By default, this step is performed on the whole aptamer sequence, but the user can choose to investigate only the variable region. In this case, the variable region will be merged with the flanking constant regions to include motifs that are across that specific portion. In the **third step**, we assume that if a motif is shared by a large fraction of the sampled aptamer pool, it is highly probable that it will emerge as highly populated, even if we do not use the entire pool of aptamer sequences. Thus, a number of aptamers and their secondary structure motifs are randomly picked from the output of the previous step and their sequences are aligned to each other to obtain a consensus sequence. The step is iteratively repeated for a number of times defined by the user. In **the last step**, each consensus sequence generated in the previous step is matched with the relative loop structure pool. Specifically, each nucleotide of the consensus sequence is matched with the respective nucleotide of the pool sequence to calculate an identity matching score based on an identity/non-identity scoring system (with a fixed gap penalty of -0.5). Finally, APTANI returns an output table, containing all results of the various steps and a clustered tree of sequences that passed the frequency filter.

### Installation
* **APTANI** is fully coded in Python version 3.3. It is supported on any Unix system and can be run also in Windows environment
* Python version 3.3 or later are mandatory to run **APTANI**
* In order to run **APTANI** it is recommend to download also Clustal Omega, available [here](https://www.clustal.org)
* Finally, for the secondary investigation of aptamers, **APTANI** exploits the software RNASubopt, from the Vienna Package, available [here](https://www.tbi.univie.ac.at/RNA/index.html#download)
* **APTANI** can be downloaded from the APTANI\_v1 directory or as a zipped file (_APTANI\_v1.zip_)

### Help and documentation
**APTANI** has been developed to cope the need for bioinformatics tools able to investigate and analyze data derived both from SELEX and HT-SELEX, focusing on aptamers. It is a single instruction software, which utilizes the Unix shell to function and APTANI provides a lot of different parameters in order to be a flexible software and to adapt to a large specter of analysis. 

#### Input parameters
The complete list and discussion of the available parameters of APTANI can be found here:

* [APTANI_v1 help](https://github.com/bicciatolab/APTANI/blob/main/APTANI_v1_help.pdf) file
* Supplementary Information material of Caroli J, Taccioli C, De La Fuente A, Serafini P, Bicciato S. APTANI: a computational tool to select aptamers through sequence-structure motif analysis of HT-SELEX data. _Bioinformatics_. 2016 Jan 15;32(2):161-4  doi: [10.1093/bioinformatics/btv545](https://doi.org/10.1093/bioinformatics/btv545)

#### Output
APTANI returns in output a comma-separated value (.csv) file for each secondary structure motif:

* Hairpins_data.csv
* Intra_Strand_data.csv
* Left_Bulges_data.csv
* Right_Bulges_data.csv

Each file contains:

1. the aptamer sequence
1. the aptamer frequency
1. the id of the aptamer cluster
1. the aptamer sub-sequence of the secondary structure motif
1. the consensus sequence
1. the alignment score between the consensus and the secondary structure motif sequences
1. the frequency of the secondary structure motif sub-sequence

## APTANI<sup>2</sup>

**APTANI<sup>2</sup>** builds on the evidence that sequence motifs with a significant binding potential are enriched at a given SELEX cycle and have an intrinsic high stability across the secondary structures of different aptamers. Therefore, in APTANI2 potentially relevant aptamers are selected using a scoring function that integrates the frequency and the stability of each motif retrieved in any aptamer sequence. The core of APTANI2 consists of three major steps:

1. quantification of aptamer frequency 
2. identification of motifs associated to secondary structures
3. scoring of sequence motifs and aptamers.

The **first step** calculates the relative frequency of each aptamer sequence produced by the SELEX process, as described in the original APTANI. In the **second step**, aptamers that pass the frequency filter are processed with RNASubopt to predict all possible secondary structures within a range of 3 Kcal/mol above the minimum free energy (MFE). Each predicted secondary structure is then investigated to retrieve all the associated hairpin, intra-strand, bulge (left and right), and G-quadruplex sequence motifs. In the **third step**, first each motif is characterized in terms of frequency and structural stability of the related secondary structure through a motif score (_MtfScore_) calculated using each motif frequency and relative stability (defined as the ratio between the median energy value of all secondary structures containing a given motif and the lowest Minimum Free Energy of all predicted secondary structures). Then, the distribution of the _MtfScore_ values is used to define the top scoring motifs as those motifs with an _MtfScore_ greater than e.g., the 99<sup>th</sup> percentile of the _MtfScore_ distribution. Finally, each aptamer is ranked according to an aptamer score _AptScore_ that quantifies, in each sequence, the normalized abundance of top scoring motifs given the total number of retrieved motifs.
APTANI<sup>2</sup> comprises three additional post-analysis modules to inspect the evolution of aptamer enrichment along the SELEX process (_Evolution Analyzer_), to retrieve motifs in aptamer sequences (_Motif Fetcher_), and to visualize their predicted folding and motif structure (_grAPhTANI_). Outputs are in the form of tables and graphical representations to facilitate the downstream interpretation of results. In addition to the command line scripts, the tool is accessible through a graphic user interface (GUI) that simplifies the input selection and parameter setting through pre-compiled fields in a point-and-click environment.

### Installation
* APTANI<sup>2</sup> is fully coded in Python version 3.4. It is supported on any Unix system and can be run also in Windows environment
* The latest Python environment can be downloaded [here](https://www.python.org/downloads/)
* For the graphical representation of secondary structures, grAPhTANI exploits VARNA, available [here](https://varna.lisn.upsaclay.fr/index.php?lang=en&page=downloads&css=varna)
* **APTANI<sup>2</sup>** can be downloaded from the APTANI2\_v1 directory or as a zipped file (_APTANI2\_v1.zip_)

### Help and documentation
**APTANI<sup>2</sup>** builds on the evidence that sequence motifs with a significant binding potential are enriched at a given SELEX cycle and have an intrinsic high stability across the secondary structures of different aptamers. **APTANI<sup>2</sup>** is a single instruction software, which utilizes the Unix shell to function. **APTANI<sup>2</sup>** provides a lot of different parameters in order to be a flexible software and to adapt to a large specter of analysis.


#### Input parameters
The complete list and discussion of the available parameters of APTANI<sup>2</sup> can be found here:

* [APTANI2_v1 help](https://github.com/bicciatolab/APTANI/blob/main/APTANI2_v1_help.pdf) file
* Supplementary Information material of Caroli J, Forcato M, Bicciato S. APTANI2: update of aptamer selection through sequence-structure analysis. _Bioinformatics_ 2020 Apr 1;36(7):2266-2268.  doi: [10.1093/bioinformatics/btz897](https://doi.org/10.1093/bioinformatics/btz897)

#### Output
Results are summarized in several tables saved as tab delimited files. In particular, **APTANI<sup>2</sup>** returns:

1. a file with the values of _MtfScore_ for each type of motif
2. the _Aptamer_Scores.txt_ file containing the _AptScore_ values and a flag indicating the presence of G-quadruplexes for each aptamer sequence
3. the _counts.txt_ file, which contains the aptamer frequencies
4. the _Motifs.txt_ file listing the motifs, and their related scores, retrieved in each aptamer
5. the _Structures_Scores.txt_, which contains all the input information needed for the visualization of the secondary structure in the _grAPhTANI module_ (i.e., the aptamer sequence, the secondary structure in dots and brackets notation, and the sequence position in base pairs of each motif identified in any predicted secondary structure)

## Example
An example data set, with scripts and output files for both **APTANI**  and **APTANI<sup>2</sup>** is available at [bicciatolab_data](https://github.com/bicciatolab/bicciatolab_data/tree/main/APTANI_data_and_scripts) 
