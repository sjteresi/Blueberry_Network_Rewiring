# Blueberry Network Rewiring
Breed new blueberry cultivars that are resistant to the blueberry stem gall wasp. Examine standing genetic variation in the population of blueberry cultivars to identify resistant blueberry cultivars and the genomic loci responsible.

[Previously](https://github.com/EdgerLab/Blueberry_RNA_Seq_Expression_Analysis), differentially expressed genes were quantified from RNA-Seq data gathered from Liberty and Draper cultivars. Work now continues on comparing metabolic pathway flux and structure following gene expression perturbation following gall wasp infection.

# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| Undergraduate: | Alder Fulton  | [Personal GitHub](https://github.com/Alder-pixel) | <fultona5@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Genome Information:
| Regular CoGe ID                   | Masked ID                                    |
|-----------------------------------|----------------------------------------------|
| Arabidopsis thaliana Col-0 (id 1) | CNS PL.20 Masked repeats 50X (v10, id 16746) |
| Vaccinium corymbosum (id 39928)   | mask w/ RepeatMasker (v3, id 58746)          |

# Code:
The code is broken up into several different scripts inside the `src/` directory. The files are as follows:

## Arabidopsis Blueberry Orthology:
- `filter_orthologs.py`: Master code file for executing and filtering gene orthology data. Creates an ortholog table by merging a set of syntelogs (SynMap) and a set of homologs (BLAST).
	- Inputs: raw syntelog data from SynMap, see [SynMap Methods](#identifying-syntelogs), raw homolog data from BLAST, see [BLAST Methods](#identifying-homologs).
	- Outputs: An ortholog table. Created so that each blueberry gene can have only 1 Arabidopsis gene match, and Arabidopsis genes can be repeated in this table (non-unique to each match).
- `import_homologs.py`: Imports the homolog data from the raw file and manages data filtration. Helper file of `filter_orthologs.py`.
- `import_syntelogs.py`: Imports syntelog data from the raw file and manages data filtration. Helper file of `filter_orthologs.py`.
- `merge_homo_synt`: Merges the sets of homologs and syntelogs. Helper file of `filter_orthologs.py`. Prioritizes results from synteny over results from simple homology.
- `blastall.sb`: BASH file that runs the BLAST search on the computing cluster.

## FPKM and TPM Evaluation:
- `process_fpkm.py`: Generates an expression table in FPKM from a count matrix of blueberry genes.
	- Inputs: Genes input file, count matrix, output directory path
	- Outputs: Saves a table of FPKM values to the output directory.
- `fpkm.py`: Helper file, calculates the FPKM (fragments per kilobase trascript per million reads).
- `process_tpm.py`: Generates an expression table in TPM from a count matrix of blueberry genes.
	- Inputs: Gene input files, count matrix, output directory path (same arguments as `process_fpkm.py`)
	- Outputs: Saves a table of TPM values to the output directory.
- `tpm.py`: Helper file, calculates the TPM matrix.
- `gene_lengths.py`: Reads a gene file and sums the length of the gene. Helper file of `process_fpkm.py` and `process_tpm.py`. Sums by exon length.
- `count_matrix.py`: Reads in and cleans the count data. Helper file of `process_fpkm.py` and `process_tpm.py`.

## WGCNA:
- `run_wgcna.sb`: Used to run `Blueberry_WGCNA.R` on the computing cluster.
- `Blueberry_WGCNA.R`: Executes WGCNA
	- Inputs: TPM matrix of blueberry genes
	- Outputs: Modules of co-expressed genes (sorted by an arbitrary color identifier)
	- Version Control: `requirements/WGCNA_sessionInfo.txt`
- Tutorial Followed:
	- I followed the tutorials found at this [link](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/), where I specifically follow **part I**. The code found in `src/WGCNA` is an adapted version of the tutorials, there are no clinical traits that I am relating modules to.

## Modules TODO:
- `filter_modules.py`: Find union of differential expression / orthology set with the WGCNA output of 10 genes assigned to modules.
	- Inputs: WGCNA output file, Syntenny/homology output file, output folder
	- Outputs: `TODO`

##  TopGO:
- `generate_gene_w_GO_term.py`: Filters a GO term data table of Arabidopsis genes
	- Inputs: GO master file, downloaded from TAIR, output directory path
	- Outputs: `ArabidopsisGene_w_GO.tsv`, each Arabidopsis gene and its corresponding GO terms, filtered down so it can be used as an input in TopGO
- `topGO_Modules.R`: Runs TopGO on the modules from WGCNA
	- Inputs: Folder of modules in Arabidposis gene format, filtered GO output from `generate_gene_w_GO_term.py`, output directory, documentation directory
	- Outputs: Table of biological process GO terms and their level of overrepresentation per module ID, session info (described below)
	- Version Control: `doc/TopGO_sessionInfo.txt`
- `topGO_DEGs.R`: Runs TopGO on the unique set of DEGs Liberty, and separately Draper
	- Inputs: Files of blueberry DEGs in Arabidposis gene format, there is command in the Makefile to make these text files. A filtered GO output from `generate_gene_w_GO_term.py`, output directory, documentation directory
	- Outputs: Table of biological process GO terms and their level of overrepresentation per DEG set, session info (described below)
	- Version Control: `doc/TopGO_sessionInfo.txt`


## Log\_2FC\_Melanie TODO flesh out more at a later date after conferring with Melanie:
- `module_representation_log_2FC.py`: Filters a set of log 2 FC (fold-change) genes that have a signal for up or down regulation in a given RNA-seq library context. Input files derived from analysis generated by Melanie, filtration and subsequent analysis performed by Scott.
	- Inputs: Folder of tables of columns (blueberry genes, arabidopsis ortholog, E value, Point of Origin, RNA-Seq library context with float value for fold change)
	- Outputs: 1 file, a table of Module ID and columns of the percent AND number of genes in a given expression context, up or down regulation.

##  gene\_stats TODO:
- `operations.py`: Calculates the percentages of each gene belonging to each identification type. Note: this is not used anywhere, and whenever a gene is found by both, Syntenny was chosen.
- `summary_table.py`: Unifies the following dataframes - differentially expressed genes, Arabidopsis ortholog, Arabidopsis GO terms, and Blueberry gene network module identity.
	- Inputs: File of bluberry genes with their module colors, file of blueberry genes and their Arabidopsis orthologs, Directory containing the differentially expressed files, file of Arabidopsis genes and its GO term list, Output filename, Output directory.
	- Outputs: Creates a csv in the output directory of the unified dataframe.

## Proteins:
- `protein_table.py`: Filters the master TAIR Arabidopsis protein table to a more manageable format. Generates a table of Arabidopsis genes, their protein ID, and their protein names.
    - Inputs: Parent path of master protein/gene table from TAIR
    - Outputs: `Protein_and_Genes_Unfiltered.txt` (TEMP), `Filtered_Arabidopsis_Protein_Info.tsv`

## Module Overlap:
- `module_go_overlap.py`: Filters several datasets to show a set of interesting GO terms and their representation in the gene modules. 
	- Inputs:  See docstring
	- Outputs:  See docstring

## QTL:
- `deg_qtl.py`: Subsets the QTL output data of interesting genes to determine
  which genes are DEGs.
  - Inputs: Path of DEG results folder from EdgeR. Path of QTL output file (txt
    file of gene names).
  - Outputs: `all_qtl_genes.tsv` a file with DEG status columns appended to the
    gene names, is each gene up or down regulated in a given DEG comparison
    context? `at_least_one_deg_qtl_genes.tsv`, a subsetted `all_qtl_genes.tsv`
    file, a gene is included in this file if it is a DEG in at least one DEG
    comparison context.


## Miscellaneous:
- `exp_table_melanie.py`: Generates an FPKM table of a blueberry gene and its syntelog.
	- Inputs: Parent path of FPKM table, parent path of syntelog table, output directory
	- Outputs: Saves FPKM table as a tsv.

## Python Requirements:
We used Pip to manage our Python packages in a virtual environment. Python package version control information can be found at `requirements/common.txt`. Python 3.6.4

# Methods Descriptions:
## Ortholog Identification
Genome-wide analyses were performed to identify *Arabidopsis thaliana* – *Vaccinium corymbosum* orthologs using a combination of synteny- and BLAST- based approaches using SynMap and reciprocal BLASTp analyses (protein databases).
Ortholog pairs with E-values greater than or equal to 0.05 were excluded from both datasets.
In the event that SynMap or BLASTp returned multiple Arabidopsis genes for a single blueberry gene, the Arabidopsis ortholog with the best (lowest) E-value score was kept.
In the event that a blueberry gene had an ortholog identified through BLASTp and SynMap, we gave priority to the SynMap synteny results, preferentially keeping that ortholog pair.

### SynMap for Syntelogs:
[SynMap](https://genomevolution.org/CoGe/SynMap.pl) was run on the CoGe website. I ran SynMap with mostly [default options](https://genomevolution.org/wiki/index.php/SynMap), I did change one option: under *Merge Syntenic Blocks* I set it to `Quota Align Merge`. Here is the [link](https://genomevolution.org/r/1ejoj) to regenerate the analysis on CoGe.

### BLASTp for Homologs:
We are doing this step to identify homologs that may have been missed using a synteny-based approach.
Genes that could have been missed by the synteny search include single-gene transpositions (and others).
We are going to use a BLAST database of protein predictions.
First we generate a BLAST database and prepare it for protein indices.
Then we can run the BLAST algorithm on the two sequence files, this may take awhile.
For notes on the options for `blastall`, please refer to the [documentation](https://www.ncbi.nlm.nih.gov/Class/BLAST/blastallopts.txt).

TODO the blastp database was generated via...

## Project News:
[Combating the blueberry stem gall wasp](https://www.canr.msu.edu/news/combating-the-blueberry-stem-gall-wasp#:~:text=The%20blueberry%20stem%20gall%20wasp%20is%20a%20tiny%20insect%20that,shoot%20and%20decreases%20fruit%20production.)
