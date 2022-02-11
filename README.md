# Blueberry Network Rewiring
Breed new blueberry cultivars that are resistant to the blueberry stem gall wasp. Examine standing genetic variation in the population of blueberry cultivars to identify resistant blueberry cultivars and the genomic loci responsible.


[Previously](https://github.com/EdgerLab/Blueberry_RNA_Seq_Expression_Analysis), differentially expressed genes were quantified from RNA-Seq data gathered from Liberty and Draper cultivars. Work now continues on comparing metabolic pathway flux and structure following gene expression perturbation following gall wasp infection.

## Project Newspiece:
[Combating the blueberry stem gall wasp](https://www.canr.msu.edu/news/combating-the-blueberry-stem-gall-wasp#:~:text=The%20blueberry%20stem%20gall%20wasp%20is%20a%20tiny%20insect%20that,shoot%20and%20decreases%20fruit%20production.)


# Contact Information:
| Role          | Name          | GitHub                                                  | Email              |
|---------------|---------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi  | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| Undergraduate: | Alder Fulton  | [Personal GitHub](https://github.com/Alder-pixel) | <fultona5@msu.edu> |
| PI:           | Patrick Edger | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |

# Genome:
| Regular CoGe ID                   | Masked ID                                    |
|-----------------------------------|----------------------------------------------|
| Arabidopsis thaliana Col-0 (id 1) | CNS PL.20 Masked repeats 50X (v10, id 16746) |
| Vaccinium corymbosum (id 39928)   | mask w/ RepeatMasker (v3, id 58746)          |

# Code:
The code is broken up into several different scripts. The files are as follows:

**/src:**  
-  **/Arabidopsis_Blueberry_Orthology**  
	- `filter_orthologs.py`: Creates an ortholog table. `TODO`: Confirm and check this in more detail.
		- Inputs: Syntelogs, Orthologs, Output directory.
		- Outputs: An ortholog table. `TODO`: Confrim this. The output does not exist on my computer, not where the make file is supposed to send it.
	- `import_homologs.py`: Imports the homolog data and provides a calss for its access. Helper file of `filter_orthologs.py`.
	- `import_syntelogs.py`: Imports syntelogs from the raw file and manages data filtration. Helper file of `filter_orthologs.py`.
	- `merge_homo_synt`: Merges the homologs and syntelogs. Helper file of `filter_orthologs.py`
	- `blastall.sb`: Batch file that runs the BLAST search on the computing cluster. `TODO`: Not sure about this one.
-  **/FPKM_TPM**  
	- `process_fpkm.py`: Saves a FPKM table
		- Inputs: Genes input file, count matrix, Output directory 
		- Outputs: Saves a table of FPKM values to the output directory.
	- `gene_lengths.py`: Reads a gene file and sums the length of the gene. Helper file of `process_fpkm.py` and `process_tpm.py`.
	- `count_matrix.py`: Reads in and cleans the count data. Helper file of `process_fpkm.py` and `process_tpm.py`.
	- `fpkm.py`: Calculates the FPKM (fragments per kilabase trascript per million reads)
	- `process_tpm.py`: Generates and saves a TPM table from gene annotation and count matrix.
		- Inputs: Gene input files, count matrix, Output directory (same arguments as `process_fpkm.py`)
		- Outputs: Saves a TPM table.
	- `tpm.py`: Calculates the TPM matrix. Helper file of `process_tpm.py`
-  **/WGCNA**  
	- `run_wgcna.sb`: Used to run `Blueberry_WGCNA.R` on the computing cluster.
	- `Blueberry_WGCNA.R`: `TODO`
-  **/modules**  
	- `filter_modules.py`: Find union of differential expression / orthology set with the WGCNA output of 10 genes assigned to modules.
		- Inputs: WGCNA output file, Syntenny/homology output file, output folder
		- Outputs: `TODO`
-  **/TopGO**  
	- `topGo_blueberry.R`: Runs TopGo
		- Inputs: Folder of modules in Arabadposis gene fromat, `TODO`: Tsv containing?, Output directory, Documentation directory
		- Outputs: `TODO`: ???
	- `generate_gene_w_GO_term.py`: Filters a GO term database of Arabidopsis genes.
		- Inputs: Go master file - downloaded from TAIR, Output directory
		- Outputs: Creates a csv derived from the GO term database. 
-  **/gene_stats**  
	- `operations.py`: Calculates the percentages of each gene belonging to each identification type. Note: this is not used anywhere, and whenever a gene is found by both, Syntenny was chosen.
-  **/requirements**  
	- `common.txt`: A list of the common requirements. See next section.
- `summary_table.py`: Unifies the following dataframes - differentially expressed genes, Arabidopsis ortholog, Arabidopsis GO terms, and Blueberry gene network module identity.
	- Inputs: File of bluberry genes with their module colors, file of blueberry genes and their Arabidopsis orthologs, Directory containing the differentially expressed files, file of Arabidopsis genes and its GO term list, Output filename, Output directory.
	- Outputs: Creates a csv in the output directory of the unified dataframe.
- `exp_table_melanie.py`: Generates an FPKM table of a blueberry gene and its syntelog.
	- Inputs: Parent path of fpkm table, parent path of syntelog table, output directory
	- Outputs: Saves FPKM table as a csv.

## Requirements:
Please install **Pip** so that you may easily install Python packages. Then use Pip to go over our `src/requirements/common.txt` and install the needed Python packages: `pip install -r common.txt`. It is wise to create a virtual environment in case you have any conflicting package installations. Please refer to the [documentation](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) for notes on that.

## At_BB:
This block of work has both non-code tasks and code tasks. First, I ran the blueberry genome through CoGe to perform a syntelog search. Second, I supplement the results with a BLAST search.

### Running SynMap:
This section describes the methods to run [SynMap](https://genomevolution.org/CoGe/SynMap.pl) on CoGe. I ran SynMap with mostly [default options](https://genomevolution.org/wiki/index.php/SynMap), I did change one option: under *Merge Syntenic Blocks* I set it to `Quota Align Merge`. Here is the [link](https://genomevolution.org/r/1ejoj) to regenerate the analysis on CoGe.

### Running BLAST:
We are doing this step to identify homologs that may have been missed using a synteny-based approach. Genes that could have been missed by the synteny search include single-gene transpositions (and others). We are going to use a BLAST database of protein predictions.

Please refer to the script at `src/At_BB/blastall.sb` for more information. First we generate a BLAST database and prepare it for protein indices. Then we can run the BLAST algorithm on the two sequence files, this may take awhile. For notes on the options for `blastall`, please refer to the [documentation](https://www.ncbi.nlm.nih.gov/Class/BLAST/blastallopts.txt).

### Filtering and Combining our Data:
We have two data files, one from Synmap (synteny) and one from BLAST (homology). Now we must filter both outputs and merge the files. I orchestrate all of the running of the code from `generate_pairs.py` but compartmentalize each function/class object as necessary into their own files. Please refer to `generate_pairs.py` for the general control flow of the project from this point on out.

**Filter the SynMap Output File**:
The output file has a lot of extraneous information. We are going to distill it down into a 2-column tab-separated values (tsv) file. To do this we will use the Pandas module in Python to maniuplate the dataframe. If you have an aversion to Python you can always accomplish this task in R. Here I accomplish this with the use of the `import_syntelogs.py` file.

**Filter the BLAST Output File**:
This output file also has a lot of extraneous information for our purposes. We distill this file down as well. Please examine `import_homologs.py` for this portion.

**Merging Syntelogs and Homologs**:
Here I use a Python program to merge the files, making sure that I don't put the same gene twice, as the BLAST search could turn up some of the same gene pairs we found in the SynMap search. We defer to taking the SynMap result (Arabidopsis - Blueberry gene-pair) over the BLAST pair result. Please examine `merge_homo_synt.py` for this portion.


**Merging Genes From Expression Analyses with the Merged Data**:
Earlier in the year I worked on calculating differential expression data from galling wasp treatments on two blueberry cultivars (Draper and Liberty, shortened to DRA and LIB in the data files). There were several timepoints and a treatments, for each comparison group I create a file in the output data folder. 

We want to examine our merged data file, and take the subset of Arabidopsis-Blueberry gene pairs that match with the blueberry gene given from the differential expression data set. We perform this subsetting in two ways. First we take the Arabidopsis gene that belongs to every Arabidopsis-Blueberry pair with a matching differential expression blueberry gene. However some of these differential expression blueberry genes may match with the same Arabidopsis gene, so we run the analysis again separately and remove any matches in which the Arabidopsis would be repeated twice.

## fpkm:

## WGCNA:
[WGNCA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/index.html) is an R package for creating **Weighted Gene Co-Expression Networks**. We are using it to identify clusters (modules) of highly correlated genes; we will use that as a lens to interpret the functions and relations of genes involved in galling response.

In regards to the code, I follow the tutorials found at this [link](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/), where I follow **part I**. The code found in `src/WGCNA` is an adapted version of the tutorials, notably there are no clinical traits that I am relating modules to.

### WGCNA Code Version Control:
* R == 4.0.2
* WGCNA == 1.69
* tidyverse == 1.3.0
* docstring == 1.0.0
