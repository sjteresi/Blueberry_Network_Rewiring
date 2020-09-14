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
The code is dividied into 3 main directories. `At_BB`, `fpkm`, and `WGCNA`. Descriptions follow:

- `At_BB`: Use the output from a BLAST database search of blueberry proteins to Arabidopsis along with the output from **SynMap** on [CoGe](https://genomevolution.org/CoGe/SynMap.pl) to identify homologs and syntelogs. Join these two datasets together and then facilitate subsetting from a list of differentially expressed genes previously calculated with EdgeR at [Blueberry Expression Analysis](https://github.com/EdgerLab/Blueberry_RNA_Seq_Expression_Analysis). Output: Lists of *Arabidopsis* genes (that match with differentially expressed genes) for each RNA-seq library.
- `fpkm`: Use the count matrix previously calculated at [Blueberry Expression Analysis](https://github.com/EdgerLab/Blueberry_RNA_Seq_Expression_Analysis) with the gene annotation to calculate FPKM.
- `WGCNA`: `TODO` add more here. Use a TPM matrix in conjunction with the `WGCNA` package in R to generate a gene network.

## Requirements:
Please install **Pip** so that you may easily install Python packages. Then use Pip to go over our `scripts/common.txt` and install the needed Python packages: `pip install -r requirements.txt`. It is wise to create a virtual environment in case you have any conflicting package installations. Please refer to the [documentation](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) for notes on that.

## At_BB:
This block of work has both non-code tasks and code tasks. First, I ran the blueberry genome through CoGe to perform a syntelog search. Second, I supplement the results with a BLAST search.

### Running SynMap:
This section describes the methods to run [SynMap](https://genomevolution.org/CoGe/SynMap.pl) on CoGe. I ran SynMap with [default options](https://genomevolution.org/wiki/index.php/SynMap). Here is the [link](https://genomevolution.org/r/1ejoj) to regenerate the analysis on CoGe.

### Running BLAST:
We are doing this step to identify homologs that may have been missed using a synteny-based approach. Genes that could have been missed by the synteny search include single-gene transpositions (and others). We are going to use a BLAST database of protein predictions.

Please refer to the script at `scripts/At_BB/blastall.sb` for more information. First we generate a BLAST database and prepare it for protein indices. Then we can run the BLAST algorithm on the two sequence files, this may take awhile. For notes on the options for `blastall`, please refer to the [documentation](https://www.ncbi.nlm.nih.gov/Class/BLAST/blastallopts.txt).

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

