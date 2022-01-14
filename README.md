# metabarcodingRpipeline

## Clone the repository:

Change the directory where you would like to clone the repository.

	$ cd my_directory

Use ``git clone`` to clone on your computer the repository including the functions and test data.

	$ git clone https://github.com/fconstancias/metabarcodingRpipeline.git


## Configure a dedicated conda envirionment:

### Install conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html>
### Create a dedicated conda environment:
	$ conda create -n metabarcodingRpipeline -y
### Activate conda environment:
	$ conda activate metabarcodingRpipeline
### install R and atropos:
	(metabarcodingRpipeline)$ conda install -c bioconda atropos -y
	(metabarcodingRpipeline)$ conda install -c conda-forge r-base -y
	(metabarcodingRpipeline)$ conda install -c conda-forge r-devtools -y
	(metabarcodingRpipeline)$ conda install -c conda-forge r-optparse -y
	
### start the R terminal:
	(metabarcodingRpipeline)$ R

### install necessary R packages - within the R terminal:
	devtools::install_github("tidyverse/tidyverse");devtools::install_github("KlausVigo/phangorn");devtools::install_github("benjjneb/dada2")

	if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
	BiocManager::install("ShortRead");BiocManager::install("DECIPHER");BiocManager::install("phyloseq")

### quit R:
	quit(save = "no")
	
### for additional functionalities (i.e., post ASV clustering using vsearch/lulu & picrust2 functional potential estimation - please install the following tools:
### R packages - similarly as done before:

	devtools::install_github("tobiasgf/lulu");devtools::install_github("ycphs/openxlsx");devtools::install_github("mikemc/speedyseq")
### vsearch in the dedicated conda environment:
	(metabarcodingRpipeline)$ conda install -c bioconda vsearch -y
	
### picrust2 in the dedicated conda environment:
	(metabarcodingRpipeline)$ conda install picrust2 -y



## Run the pipeline:

activate the dedicated conda environment:

	$ conda activate metabarcodingRpipeline


Use ``Rscript`` to run the pipeline and specify some necessary parameters e.g.: *databases* 

- ``dada`` method: <https://benjjneb.github.io/dada2/training.html>


```bash

(metabarcodingRpipeline)$ Rscript scripts/dada2_metabarcoding_pipeline.Rscript \
-i test-data/ \
--preset V3V4 \
-T 8 \
--db ~/db/DADA2/silva_nr99_v138_train_set.fa.gz \
--db_species ~/db/DADA2/silva_species_assignment_v138.fa.gz \
--metadata test-data/metadata.xlsx \
--run_phylo FALSE \
--save_out test_pipe_Rscript.RDS \
-f scripts/functions_export_simplified.R > mylogs.txt 2>&1

```
		
The ``> mylogs.txt 2>&1`` trick will redirect what is printed on the screen to a file including potential errors and also parameters that you used.

## Add phylogenetic tree to a phyloseq object:

activate the dedicated conda environment:

	$ conda activate metabarcodingRpipeline

By default, based on <https://f1000research.com/articles/5-1492>. It might take quite some time depending on your configuration and the overall richness of your phyloseq object, that's why it is not computed by default running the pipeline.

	(metabarcodingRpipeline)$ Rscript scripts/add_phylogeny_to_phyloseq.R \
		-p dada2/physeq.RDS \
		-o dada2/physeq_phylo 


## Add/replace taxonomical information fron a phyloseq object:
### using dada2 assigntaxa/ assignspecies:

#### using Rscript
```bash
(metabarcodingRpipeline)$ Rscript scripts/run_phyloseq_dada2_tax.Rscript \
--phyloseq_path rscript-output/03_dada2_merged_runs_chimera_removed/physeq.rds \
--tax_threshold 60 \
--output rscript-output/04_dada2_taxonomy \
--db ~/db/DADA2/silva_nr99_v138_train_set.fa.gz \
--db_species ~/db/DADA2/silva_species_assignment_v138.fa.gz \
--reverse_comp TRUE \
-T 4 \
-f scripts/functions_export_simplified.R
```
#### within R
```r
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/dev/scripts/functions_export_simplified.R")

readRDS("/Users/physeq.RDS") %>%
  phyloseq_dada2_tax(physeq = ., 
                     threshold = 60, # 60 (very high),  50 (high), PB = 10
                     db ="~/db/DADA2/silva_nr99_v138_train_set.fa.gz"
                     db_species ="~/db/DADA2/silva_species_assignment_v138.fa.gz"
                     nthreads = 2,
                     tryRC = TRUE,
                     return = TRUE,
                     full_return = FALSE) -> physeq_new_tax

```
                     
### using DECIPHER IDtaxa::

You can also add or replace taxonomical information of your phyloseq object using DECIPHER IDtaxa function.


```bash
Rscript scripts/run_phyloseq_DECIPHER_tax.Rscript \
--phyloseq_path rscript-output/03_dada2_merged_runs_chimera_removed/physeq.rds \
--export rscript-output/04_dada2_taxonomy \
--reverse_comp TRUE \
--db ~/db/DADA2/SILVA_SSU_r132_March2018.RData \
--tax_threshold 50 \
-T 8 \
-f scripts/functions_export_simplified.R
```
#### within R

```r
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/dev/scripts/functions_export_simplified.R")

readRDS("/Users/physeq.RDS") %>%
  phyloseq_DECIPHER_tax(physeq = ., 
                        threshold = 60, # 60 (very high),  50 (high), PB = 10
                        db="~/db/DADA2/SILVA_SSU_r132_March2018.RData" 
                        nthreads = 6,
                        tryRC = TRUE,
                        return = TRUE) -> physeq_new_tax
```

## Add/update metadata a phyloseq object:


```r
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/dev/scripts/functions_export_simplified.R")

ps_tax_phylo %>%
  physeq_add_metadata(physeq = .,
                      metadata = "test-data/metadata.xlsx" %>%
  readxl::read_xlsx(),
                      sample_column = "sample_name") -> ps_tax_phylo_meta

ps_tax_phylo_meta
```


## TODO:

- <s>add phylogenetic tree to a phyloseq object</s>
- change name
- add https://zenodo.org/account/settings/github/ -> DOI
- Faster method for phylogenetic reconstruction (e.g., MAFFT + FastTree)
- add possibility to skip primer removal: skipping run_atropos() or changing atropos parameter?
- <s>replace taxonomic assignments of a phyloseq object using alternative approach/ database</s>
- <s>cluster ASV using DECIPHER</s>
- <s>cluster ASV using vsearch lulu</s>
- <s>run picrust2 from a phyloseq object</s>

## Help:


activate the dedicated conda environment:

	$ conda activate metabarcodingRpipeline

print help:
	
	(metabarcodingRpipeline)$ Rscript scripts/dada2_metabarcoding_pipeline.R --help

	Usage: scripts/dada2_metabarcoding_pipeline.R [options]


	Options:
	
	-i CHARACTER, --input_directory=CHARACTER
		Path of the input directory containing raw _R1_ and _R2_ raw reads in their respective run sub-directories 

       e.g., -i raw [contains raw/run1 and raw/run2]

       N.B.: sample name is extracted from .fastq.gz samples  before the first '_' e.g., XXX-XX 

       sample names must be unambigious and unique e.g., sample-1 = sample-11, sample-01 != sample-10

	-a CHARACTER, --atropos_binary=CHARACTER
		Path of atropos program [used for primer removal]

	-o CHARACTER, --output_directory=CHARACTER
		Name of the output directory

	-V CHARACTER, --pipeline=CHARACTER
		V4 or V3V4 will use default primers and parameters as used in the FBT lab [primers, trunc, maxee, overlap, expected length, ...]

	-t CHARACTER, --tax_method=CHARACTER
		User can specify using dada2 (=dada) or DECIPHER for taxonomic assignments [default dada]

	--tax_threshold=NUMERIC
		Thershold for taxonomic assignments [if --tax_metod dada: he minimum bootstrap confidence for assigning a taxonomic level. if --tax_method DECIPHER: Numeric specifying the confidence at which to truncate the output taxonomic classifications. ]

	--metadata=CHARACTER
		Path to excel document containing metadata [Sample identifier column should be sample_name]

	--database=CHARACTER
		Path to the taxonomic database

	--database_for_species_assignments=CHARACTER
		Path to the speies-level taxonomic database [only for --tax_metod  dada]

	--phylo=CHARACTER
		Compute phylogenetic tree from the ASV sequence ?

	--PRIMER_F=CHARACTER
		Sequence of the gene specific Fwd primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]

	--PRIMER_R=CHARACTER
		Sequence of the gene specific Rev primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]

	--minover=NUMERIC
		Minimum overlap for merginf R1 and R2 reads [if using -V V4 or V3V4, this parameter is already set]

	--trunclen=NUMERIC
		Nucleotide position to truncate the Fwd and Rev reads at [if using -V V4 or V3V4, this parameter is already set]

	--trim_length=NUMERIC
		ASV of length outside the range will be discarded [i.e., insilco size exclusion of ASV - if using -V V3 or V3V4, this parameter is already set]

	--maxee=NUMERIC
		Maximum expected error for Fwd and Rev reads [if using -V V4 or V3V4, this parameter is already set]

	--minLen=NUMERIC
		Minimul read length [if using -V V3 or V3V4, this parameter is already set]

	-T NUMERIC, --slots=NUMERIC
		Number of threads to perform the analyses

	-h, --help
		Show this help message and exit

