# metabaRpipe:

This pipeline is mainly based on `dada2` and the follwoing tutorials [https://f1000research.com/articles/5-1492](https://f1000research.com/articles/5-1492), [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html). 

Cite authors who deserve credits for their valuable work!

## Installation:
### Configure a dedicated conda environment:

1. [Install conda](<https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html>)
*e.g.*, for MAC users:`bash Miniconda3-latest-MacOSX-x86_64.sh`

2. Create a dedicated conda environment:


```bash
	conda create -n metabaRpipe -y
```

3. Activate the conda environment:

```bash
	conda activate metabaRpipe
```
4. Install `R` and `atropos` and `devtools` R package:

```bash
	conda install -c conda-forge r-base=4.1 -y
	conda install -c bioconda atropos=1.1.25 -y
	conda install -c conda-forge r-devtools -y
```
5. Confirm R was correctly installed within the conda environment:

```bash
 	which R
```	
which should result in something like this, indicating you will use R installed within your conda environment.

	/Users/xxx/miniconda3/envs/metabaRpipe/bin/R
	
6. Start `R` from the terminal an install the required `R packages`:

```R
	R
	
	install.packages("optparse", repos = "https://cloud.r-project.org")
	devtools::install_github('tidyverse/tidyverse')
	devtools::install_github("benjjneb/dada2")
	devtools::install_github("KlausVigo/phangorn")

	install.packages("BiocManager", repos = "https://cloud.r-project.org")
	BiocManager::install("ShortRead", update = FALSE)
	BiocManager::install("DECIPHER", update = FALSE)
	BiocManager::install("phyloseq", update = FALSE)

	quit(save = "no")
```

7. Clone the `metabaRpipe` repository

This step is required to get the `Rscripts` locally.

Change the directory where you would like to clone the repository.

```bash
	MY_DIR=/path_to_mydir/whereIwant/metabaRpipe/to/be/installed/
	cd ${MY_DIR}
```
Use ``git clone`` to clone on your computer the repository including the functions and test data.

```bash
	git clone https://github.com/fconstancias/metabaRpipe.git
```

<p align="right">(<a href="#top">back to top</a>)</p>

# Running the pipeline:

Everything is now ready to analyse your raw data. We will from the terminal using `Rscripts`, enabling you to run the pipeline and generate a phyloseq object from raw sequencing data using one single command.


* First, activate the dedicated conda environment:

```bash
	conda activate metabaRpipe
```

* Use ``Rscript`` to run the pipeline and specify some parameters e.g.: *databases* 

For instance, using the `test-data` available from the repository:

```bash
Rscript ${MY_DIR}/metabaRpipe/Rscripts/dada2_metabarcoding_pipeline.Rscript \
-i ${MY_DIR}/metabaRpipe/test-data/ \
--preset V3V4 \
--db ${MY_DIR}//metabaRpipe/databases/databases/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz \
--db_species ${MY_DIR}/metabaRpipe/databases//databases/GTDB_bac120_arc122_ssu_r202_Species.fa.gz \
--metadata ${MY_DIR}/metabaRpipe/metadata.xlsx \
--save_out test_pipe_Rscript.RDS \
-f ${MY_DIR}/metabaRpipe/Rscripts/functions.R  > mylogs.txt 2>&1
```		
The ``> mylogs.txt 2>&1`` trick will redirect what is printed on the screen to a file including potential errors and also parameters that you used.

* If you encounter a `Permission denied` error:

```bash
/metabaRpipe/Rscripts/dada2_metabarcoding_pipeline.Rscript: Permission denied

```
Mark the file as executable using `chmod`.

```bash
chmod +x ${MY_DIR}metabaRpipe/Rscripts/dada2_metabarcoding_pipeline.Rscript
```

The scripts can also be run within` R` using the R functions stored under `${MY_DIR}metabaRpipe/Rscripts/functions.R` - examples will come later.
<p align="right">(<a href="#top">back to top</a>)</p>

## Exploring the outputs:

By default the script generates several outputs under a ```dada2``` directory. Below the subfolders and outputs generated by the pipeline:

```bash
dada2/
├── 00_atropos_primer_removed
│   ├── 180914_M00842_0310_000000000-C3GBT
│   │   ├── R1F1-S66_primersout_R1_.fastq.gz
│   │   ├── R1F1-S66_primersout_R2_.fastq.gz
│   │   ├── R1F2-S300_primersout_R1_.fastq.gz
│   │   ├── R1F2-S300_primersout_R2_.fastq.gz
│   │   ├── R1F3-S90_primersout_R1_.fastq.gz
│   │   └── R1F3-S90_primersout_R2_.fastq.gz
│   └── 190719_M00842_0364_000000000-CKGHM
│       ├── Y2A15-2M-06-S78_primersout_R1_.fastq.gz
│       ├── Y2A15-2M-06-S78_primersout_R2_.fastq.gz
│       ├── Y2A15-2M-12-S77_primersout_R1_.fastq.gz
│       ├── Y2A15-2M-12-S77_primersout_R2_.fastq.gz
│       ├── Y3-R1F4-S136_primersout_R1_.fastq.gz
│       └── Y3-R1F4-S136_primersout_R2_.fastq.gz
├── 01_dada2_quality_profiles
│   ├── 180914_M00842_0310_000000000-C3GBT
│   │   ├── 180914_M00842_0310_000000000-C3GBT_forward.pdf
│   │   └── 180914_M00842_0310_000000000-C3GBT_reverse.pdf
│   └── 190719_M00842_0364_000000000-CKGHM
│       ├── 190719_M00842_0364_000000000-CKGHM_forward.pdf
│       └── 190719_M00842_0364_000000000-CKGHM_reverse.pdf
├── 02_dada2_filtered_denoised_merged
│   ├── 180914_M00842_0310_000000000-C3GBT
│   │   ├── 180914_M00842_0310_000000000-C3GBT.RData
│   │   ├── 180914_M00842_0310_000000000-C3GBT_seqtab.rds
│   │   ├── 180914_M00842_0310_000000000-C3GBT_track_analysis.tsv
│   │   ├── errors_180914_M00842_0310_000000000-C3GBT_fwd.pdf
│   │   ├── errors_180914_M00842_0310_000000000-C3GBT_rev.pdf
│   │   └── seq_distrib_180914_M00842_0310_000000000-C3GBT.pdf
│   └── 190719_M00842_0364_000000000-CKGHM
│       ├── 190719_M00842_0364_000000000-CKGHM.RData
│       ├── 190719_M00842_0364_000000000-CKGHM_seqtab.rds
│       ├── 190719_M00842_0364_000000000-CKGHM_track_analysis.tsv
│       ├── errors_190719_M00842_0364_000000000-CKGHM_fwd.pdf
│       ├── errors_190719_M00842_0364_000000000-CKGHM_rev.pdf
│       └── seq_distrib_190719_M00842_0364_000000000-CKGHM.pdf
├── 03_dada2_merged_runs_chimera_removed
│   ├── no-chim-seqtab.fasta
│   ├── no-chim-seqtab.rds
│   ├── physeq.rds
│   ├── seqtab_distrib.pdf
│   └── track_analysis.tsv
├── 04_dada2_taxonomy
│   ├── silva_nr99_v138_train_set_assignation.rds
│   └── silva_nr99_v138_train_set_table.tsv
└── phyloseq.RDS
```	

There are several outputs generated during the steps performed:

* the primers removal: `00_atropos_primer_removed` 
* the inspection of the quality profiles: `01_dada2_quality_profiles`
* the quality filtering, error learning, ASV inference and Fwd and Rev reads merging: `02_dada2_filtered_denoised_merged`
* the combinaiton of the data proceed from different runs and the removal of chimeras: `03_dada2_merged_runs_chimera_removed`
* the taxonomy assignments: `04_dada2_taxonomy`
* a `phyloseq` object has been generated and is waiting for you with the provided metadata:

```R
library(phyloseq)
readRDS("dada2/phyloseq.RDS")
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 322 taxa and 6 samples ]
sample_data() Sample Data:       [ 6 samples by 18 sample variables ]
tax_table()   Taxonomy Table:    [ 322 taxa by 7 taxonomic ranks ]
refseq()      DNAStringSet:      [ 322 reference sequences ]
```
This object contains the Amplicon Sequence Variants (ASV), their sequences `refseq()`, the ASV/sample count table `otu_table()`, the taxonomic path of the ASV `tax_table()` and the metadata `sample_data()`. This allows easy handling of the data.

More informations regarding `phyloseq` object can be found [here](https://joey711.github.io/phyloseq/).

If `--save_out test_pipe_Rscript.RDS` is specified, all the outputs are also saved within this `R` object.

```R
readRDS("test_pipe_Rscript.RDS") -> out
ls(out)

[1] "filtering_denoising" "merging"             "physeq"             
[4] "qplot"               "taxo" 
```

For more details, please check the [dada2 orignial tutorial](https://benjjneb.github.io/dada2/tutorial.html).
<p align="right">(<a href="#top">back to top</a>)</p>

## Raw data organisation:

Organisation of the raw sequencing data is crucial.

Fwd and Rev reads (*_R1_* and *_R2_*, respectively) are placed in run specific directory since error learning and ASV inference has to be perform on a run basis. If you are only analyzing one sequencing run, simply add only one subdirectory.


```bash
test-data/
├── 180914_M00842_0310_000000000-C3GBT
│   ├── R1F1-S66_L001_R1_001.fastq.gz
│   ├── R1F1-S66_L001_R2_001.fastq.gz
│   ├── R1F2-S300_L001_R1_001.fastq.gz
│   ├── R1F2-S300_L001_R2_001.fastq.gz
│   ├── R1F3-S90_L001_R1_001.fastq.gz
│   └── R1F3-S90_L001_R2_001.fastq.gz
├── 190719_M00842_0364_000000000-CKGHM
│   ├── Y2A15-2M-06-S78_L001_R1_001.fastq.gz
│   ├── Y2A15-2M-06-S78_L001_R2_001.fastq.gz
│   ├── Y2A15-2M-12-S77_L001_R1_001.fastq.gz
│   ├── Y2A15-2M-12-S77_L001_R2_001.fastq.gz
│   ├── Y3-R1F4-S136_L001_R1_001.fastq.gz
│   └── Y3-R1F4-S136_L001_R2_001.fastq.gz
└── metadata.xlsx
```

By default, sample names are retrieved from files names after the first ```_```, Fwd and Rev fastq files should be recognised by ```*_R1_*``` and ```*_R2_*```, respectively and a ```sample_name``` column in the metadata```.xlsx``` file links the sample names retrieved from the files together with the metadata.

<p align="right">(<a href="#top">back to top</a>)</p>

## Define your own presets:

The ```--preset``` option is an important parameter here since it will be used to remove primers used for PCR amplification - as recommended by dada2 authors -  and also adjust the length filtering and trimming parameters. For instance, ```V3V4``` is related to the following primers, dada2 options.

Depending on the library preparation protocol you could also skip the PCR primer removal step using `--rm_primers FALSE`


```bash
open ${MY_DIR}/metabaRpipe/Rscripts/functions.R 
...
  if(V == "V3V4"){
    
    PRIMER_F = "CCTAYGGGRBGCASCAG"
    PRIMER_R = "GGACTACNNGGGTATCTAAT"
    trim_length = c(240,600)
    trunclen =  c(260,250)
    maxee = c(4,5)
    minLen = 160
    minover = 10
    }
...
```

The sequences ```CCTAYGGGRBGCASCAG``` and ```GGACTACNNGGGTATCTAAT``` will be searched by ```atropos``` as a Fwd and Rev primers, respectively and removed. **Only** reads containing the expected primers will be kept.


Fwd and Rev reads will by truncated after ```260``` and ```250``` nucleotide positions, reads shorter then `160` nucleotides will be removed as well as the Fwd with a maximum expected error more then `4` and Rev of `5`. `10` nucleotides will be used to merged denoised Fwd and Rev reads and only ASV `>240 length <600` will be kept.

For more details, please check the dada2 orignial [tutorial](https://benjjneb.github.io/dada2/tutorial.html).

You can modify the `${MY_DIR}/metabaRpipe/Rscripts/functions.R` script by adding any `preset` of your choice adding another if statement and using the same variable names. Then  `--preset mypreset` can be used to call the parameters you defined in `${MY_DIR}/metabaRpipe/Rscripts/functions.R` when running the script ```Rscript ${MY_DIR}/metabaRpipe/Rscripts/dada2_metabarcoding_pipeline.Rscript ```. This allows automatisation of the process and the possiblity to run the pipeline from your HPC cluster.

<p align="right">(<a href="#top">back to top</a>)</p>

## ETH FBT users:

You can use the following commands to process the data targeting 16S V4 region using the PCR primers and parameters we use in the FBT group from the C18 bioinformatic workstation. Everything is installed and configured, you can follow the steps below to analyse your data.

* Get familiar with the pipeline and imporant information `p/Documentation/ILLUMINA_Sequencing/16S-bioinformatic-pipeline/FBT_bioinfo_pipeline.pdf`.

* Book the  C18 bioinformatic workstation using the [calendar](https://hest.sp.ethz.ch/bt-calendars/_layouts/15/start.aspx#/Lists/D42%20PharmaBiome%20office).

* The information to access the C18 bioinformatic workstation can be found here: `p/Documentation/ILLUMINA_Sequencing/16S-bioinformatic-pipeline/FBT_bioinfo_pipeline.pdf` on `slide 16`.  

* Create a directory were you are going to place the raw sequencing data and generate the outputs under `/Users/localadmin/WORKSHOP`. You could do that using Mac Finder or from the terminal..

* Organise the raw sequencing files in **Miseq-run-specific** **sub-directories**  as explained in the **Raw data structure** section above. It should look like `/Users/localadmin/WORKSHOP/My_dir/My_analysis/raw/my_run_1/` and if you are analysing more than 1 run you would other sub-directories, e.g., `/Users/localadmin/WORKSHOP/My_dir/My_analysis/raw/my_run_2/`

* [Open the terminal from the MAC dock](https://support.apple.com/en-gb/guide/terminal/apd5265185d-f365-44cb-8b09-71a064a42125/2.9/mac/10.14).

* Using the terminal, navigate to the directory you have created `cd` *i.e.*, change directory command.

```bash
cd /Users/localadmin/WORKSHOP/My_dir/My_analysis/
```

* Activate the conda environment. 

```bash
conda activate metabarcodingRpipeline
```

* Run the pipeline with the default V4 FBT parameters - press enter to start:

```bash
Rscript /Users/localadmin/ENGINEs/metabarcodingRpipeline/scripts/dada2_metabarcoding_pipeline.Rscript \
-i raw/ -T 4 \
--db  /Users/localadmin/ENGINEs/NEWPIPE/db/silva_nr99_v138.1_train_set.fa.gz \
--db_species /Users/localadmin/ENGINEs/NEWPIPE/db/silva_species_assignment_v138.1.fa.gz \
-f /Users/localadmin/ENGINEs/metabarcodingRpipeline/scripts/functions_export_simplified.R \
--metadata mapping_file.xlsx > run_pipe_logs.txt 2>&1
```

* Add a phylogenetic tree of the ASV directly to the R phyloseq object:

```bash
Rscript /Users/localadmin/ENGINEs/metabarcodingRpipeline/scripts/run_add_phylogeny_to_phyloseq.Rscript -p dada2/phyloseq.RDS -o dada2/phyloseq_phylo -f /Users/localadmin/ENGINEs/metabarcodingRpipeline/scripts/functions_export_simplified.R > add_phylo_logs.txt 2>&1

```

* Export qiime2 compatible files:

```bash
Rscript /Users/localadmin/ENGINEs/metabarcodingRpipeline/scripts/phyloseq_export_qiime.Rscript -i dada2/phyloseq_phylo/phyloseq_phylo.RDS -o dada2/qiime2 -f /Users/localadmin/ENGINEs/metabarcodingRpipeline/scripts/functions_export_simplified.R
```

You now have everything ready for analysis using your preferred platform !

<p align="right">(<a href="#top">back to top</a>)</p>

## Getting some help:


* Options from a specific Rscript can be access using `--help` argument.


```bash
${MY_DIR}/metabaRpipe/Rscripts/dada2_metabarcoding_pipeline.Rscript --help

Usage: /Users/test/Documents/GitHub/metabaRpipe/Rscripts/dada2_metabarcoding_pipeline.Rscript [options]


Options:
	-i CHARACTER, --input_directory=CHARACTER
		Path of the input directory containing raw _R1_ and _R2_ raw reads in their respective run sub-directories 

              e.g., -i raw [contains raw/run1 and raw/run2]

              N.B.: sample name is extracted from .fastq.gz samples  before the first '_' e.g., XXX-XX 

              sample names must be unambigious and unique e.g., sample-1 = sample-11, sample-01 != sample-10

	--raw_file_pattern=CHARACTER
		Patterns for R1 and R2 for raw files

	-a CHARACTER, --atropos_binary=CHARACTER
		Path of atropos program [used for primer removal]

	--cut_dir=CHARACTER
		Directory for atropos PCR primer removed reads

	--cut_pattern=CHARACTER
		Patterns for R1 and R2 of atropos files

	--rm_primers=CHARACTER
		Perform atropos PCR primer check & removal ?

	--filt_dir=CHARACTER
		Directory for filtered reads

	--merged_run_dir=CHARACTER
		Directory for run merging, chimera removal

	--taxa_dir=CHARACTER
		Directory for taxonomy step

	--sep=CHARACTER
		regex pattern to identify sample names  [default: after first _]

	--preset=CHARACTER
		Will use default primers and parameters as defined in the run_dada2_pipeline() [primers, trunc, maxee, overlap, expected length, ...]

	--PRIMER_F=CHARACTER
		Sequence of the gene specific Fwd primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]

	--PRIMER_R=CHARACTER
		Sequence of the gene specific Rev primer to be removed with atropos [if using -V V4 or V3V4, this parameter is already set]

	--tax_threshold=CHARACTER
		Threshold for taxonomic assignments [the minimum bootstrap confidence for assigning a taxonomic level.

	--nbases=CHARACTER
		Number of bases for error learning step

	--pool=CHARACTER
		Pooling strategy

	--minover=CHARACTER
		Minimum overlap for merginf R1 and R2 reads [if using -V V4 or V3V4, this parameter is already set]

	--trunclen=CHARACTER
		Nucleotide position to truncate the Fwd and Rev reads at [if using -V V4 or V3V4, this parameter is already set]

	--trim_length=CHARACTER
		ASV of length outside the range will be discarded [i.e., insilco size exclusion of ASV - if using -V V3 or V3V4, this parameter is already set]

	--maxee=CHARACTER
		Maximum expected error for Fwd and Rev reads [if using -V V4 or V3V4, this parameter is already set]

	--minLen=CHARACTER
		Minimul read length [if using -V V3 or V3V4, this parameter is already set]

	--chimera_method=CHARACTER
		 Chimera removal strategy

	--collapseNoMis=CHARACTER
		 Perform 100% similarity ASV clustering?

	--tryRC=CHARACTER
		Perform taxonomic assignments also on reversed complemented ASV sequence?

	--metadata=CHARACTER
		Path to excel document containing metadata [Sample identifier column should be sample_name]

	--db=CHARACTER
		Path to the taxonomic database

	--db_species=CHARACTER
		Path to the speies-level taxonomic database [only for --tax_metod  dada]

	--run_phylo=CHARACTER
		Compute phylogenetic tree from the ASV sequence ?

	--merge_phyloseq_export=CHARACTER
		Path fof the phyloseq object

	--output_phyloseq_phylo=CHARACTER
		Path fof the phyloseq object

	--save_out=CHARACTER
		Path fof the R object containing all output from the pipeline

	--export=CHARACTER
		Export pdf, table and intermediate RDS files?

	--remove_input_fastq=CHARACTER
		Remove intermediate fastq.gz files

	-T NUMERIC, --slots=NUMERIC
		Number of threads to perform the analyses

	--seed_value=NUMERIC
		Seed value for random number generator

	-f CHARACTER, --fun_dir=CHARACTER
		Directory containing the R functions

	-h, --help
		Show this help message and exit

```


<p align="right">(<a href="#top">back to top</a>)</p>


## Additional functionalities:

For additional functionalities (*i.e.*, post ASV clustering using vsearch/lulu & picrust2 functional potential estimation - please install the following tools:

### 1.  Computing an ASV-phylogenetic tree in a phyloseq object:

Based on the approach described[ here](<https://f1000research.com/articles/5-1492>). It can be intensive depending on your setup and the overall richness of your phyloseq object. This is therefore not computed by default running the pipeline.

```bash
	conda activate metabaRpipe
	Rscript ${MY_DIR}/metabaRpipe/Rscripts/run_add_phylogeny_to_phyloseq.Rscript \
	-p dada2/phyloseq.RDS \
	-o dada2/phyloseq_phylo \
	-f ${MY_DIR}/metabaRpipe/Rscripts/functions.R
```

N.B.: The phyloseq object should include the ASV sequences stored as `refseq()` :

```r
> readRDS("dada2/phyloseq.RDS")phyloseq-class experiment-level objectotu_table()   OTU Table:          [ 322 taxa and 6 samples ]:sample_data() Sample Data:        [ 6 samples by 18 sample variables ]:tax_table()   Taxonomy Table:     [ 322 taxa by 7 taxonomic ranks ]:refseq()      DNAStringSet     :      [ 322 reference sequences ]taxa are rows
```

<p align="right">(<a href="#top">back to top</a>)</p>

### 2. Adding/replacing taxonomical table from a phyloseq object:

* The `run_phyloseq_dada2_tax.Rscript` Rscript allow to update the taxonomy using `dada2::assignTaxonomy` and `dada2::assignSpecies` from the terminal. 

The `dada2 formated` databases can be downloaded [here](https://benjjneb.github.io/dada2/training.html).

Below the options we can specify to the script:


```
Rscript ${MY_DIR}/metabaRpipe/Rscripts/run_phyloseq_dada2_tax.Rscript --help

############################################################
Starting 
############################################################

Loading required package: optparse
Usage: /Users/test//metabaRpipe/Rscripts/run_phyloseq_dada2_tax.Rscript [options]


Options:
	-p CHARACTER, --phyloseq_path=CHARACTER
		Path of the input phyloseq object

	-o CHARACTER, --output=CHARACTER
		Name of the output directory

	--tax_threshold=NUMERIC
		Minimum bootstrap value / threshold

	--db=CHARACTER
		path of database

	--db_species=CHARACTER
		path of species level database

	--reverse_comp=CHARACTER
		Reverse complement sequences for taxonomic assignments?

	-T CHARACTER, --slots=CHARACTER
		Number of threads to perform the analyses

	-f CHARACTER, --fun_dir=CHARACTER
		Directory containing the R functions

	--seed=CHARACTER
		seed number for random generator

	-h, --help
		Show this help message and exit
```

To illustrate this we can update the phyloseq object we generated using the [Human Oral Microbiome Database](https://www.homd.org/).

```bash
conda activate metabaRpipe
Rscript ${MY_DIR}/metabaRpipe/Rscripts/run_phyloseq_dada2_tax.Rscript \
--phyloseq_path dada2/phyloseq.RDS \
--tax_threshold 60 \
--output 04_dada2_taxonomy \
--db ${MY_DIR}/metabaRpipe/databases/eHOMD_RefSeq_dada2_V15.22.fasta.gz \
--db_species ${MY_DIR}/metabaRpipe/databases/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz \
--reverse_comp TRUE \
-T 4 \
-f ${MY_DIR}/metabaRpipe/Rscripts/functions.R
```
We could also run the function directly within R/ Rstudio:

Below the orignial `tax_table()`

```r
require(tidyverse); require(phyloseq)

> readRDS("dada2/phyloseq.RDS") %>% tax_table() %>%  head()Taxonomy Table:     [ 6 taxa by 7 taxonomic ranks ]:       Kingdom  Phylum         Class    Order   Family   Genus   Species        <chr>    <chr>          <chr>    <chr>   <chr>    <chr>   <chr>   ASV001 Bacteria Proteobacteria Gammapr… Burkho… Burkhol… Ralsto… ""      ASV002 Bacteria Proteobacteria Alphapr… Rhizob… Rhizobi… Aureim… ""      ASV003 Bacteria Proteobacteria Gammapr… Burkho… Burkhol… Massil… ""      ASV004 Bacteria Proteobacteria Gammapr… Burkho… Burkhol… Parabu… "unknow…ASV005 Bacteria Proteobacteria Gammapr… Burkho… Burkhol… Massil… "sp0128…ASV006 Bacteria Proteobacteria Gammapr… Burkho… Burkhol… Massil… "" 
```

```r
source("https://raw.githubusercontent.com/fconstancias/metabaRpipe/master/Rscripts/functions.R")

readRDS("dada2/phyloseq.RDS") %>%
  phyloseq_dada2_tax(physeq = ., 
                     threshold = 60, 
                     db ="~/metabaRpipe/databases/eHOMD_RefSeq_dada2_V15.22.fasta.gz", #check the path
                     db_species ="~/metabaRpipe/databases/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz", #check the path
                     nthreads = 2,
                     tryRC = TRUE,
                     return = TRUE,
                     full_return = FALSE) -> physeq_new_tax
```

The `eHOMD_RefSeq_dada2_V15.22` updated  `tax_table()`

```r
> physeq_new_tax %>% tax_table() %>%  head()Taxonomy Table:     [ 6 taxa by 7 taxonomic ranks ]:        Kingdom  Phylum         Class    Order    Family   Genus  Species        <chr>    <chr>          <chr>    <chr>    <chr>    <chr>  <chr>  ASV0001 Bacteria Actinobacteria Actinob… Bifidob… Bifidob… Bifid… unknownASV0002 Bacteria Proteobacteria Gammapr… Enterob… Enterob… Esche… unknownASV0003 Bacteria Bacteroidetes  Bactero… Bactero… Prevote… Prevo… unknownASV0004 Bacteria Bacteroidetes  Bactero… Bactero… Prevote… Prevo… unknownASV0005 Bacteria Bacteroidetes  Bactero… Bactero… Bactero… Bacte… unknownASV0006 Bacteria Firmicutes     Bacilli  Lactoba… Strepto… Strep… unknown
```
          
* You can also perform taxonomic assignment of the ASV sequences using [DECIPHER IDtaxa function](https://github.com/benjjneb/dada2/issues/683). 



The `DECIPHER formated` databases (*i.e.*, Training sets for organismal classification (nucleotides) can be downloaded [here](http://www2.decipher.codes/Downloads.html).

As we have seen before, this can be performed from the terminal:

```bash
conda activate metabaRpipe

Rscript ${MY_DIR}/run_phyloseq_DECIPHER_tax.Rscript \
--phyloseq_path dada2/phyloseq.RDS \
--export dada2/04_dada2_taxonomy \
--reverse_comp TRUE \
--db ~/db/DADA2/SILVA_SSU_r132_March2018.RData \ # check where you downloaded your database
--tax_threshold 60 \
-T 4 \
-f ${MY_DIR}/metabaRpipe/Rscripts/functions.R
```

or from R/ Rstudio.

```r
require(tidyverse); require(phyloseq)

source("https://raw.githubusercontent.com/fconstancias/metabaRpipe/master/Rscripts/functions.R")


readRDS("dada2/physeq.RDS") %>%
  phyloseq_DECIPHER_tax(physeq = ., 
                        threshold = 60,
                        db="~/db/DADA2/SILVA_SSU_r132_March2018.RData"  # check where you downloaded your database
                        tryRC = TRUE,
                        return = TRUE) -> physeq_new_tax
```
 <p align="right">(<a href="#top">back to top</a>)</p>

### 3. Adding/updating the metadata information from a phyloseq object:

This can be done within R/Rstudio:

```r
require(tidyverse); require(phyloseq)

source("https://raw.githubusercontent.com/fconstancias/metabaRpipe/master/Rscripts/functions.R")

readRDS("dada2/phyloseq.RDS") %>%
  physeq_add_metadata(physeq = .,
                      metadata = "~/metabaRpipe/test-data/metadata.xlsx" %>%
                      readxl::read_xlsx()) -> ps_tax_phylo_meta

ps_tax_phylo_meta
```

 <p align="right">(<a href="#top">back to top</a>)</p>

#### 4. Perform post-clustering curation: <https://github.com/tobiasgf/lulu>

Install the following tools/ packages:

`R packages` - similarly as done before:

```bash
conda activate metabaRpipe
R
devtools::install_github("tobiasgf/lulu")
devtools::install_github("mikemc/speedyseq")
```
`vsearch` in the dedicated conda environment:

```bash
conda activate metabaRpipe 
conda install -c bioconda vsearch -y
```

```r
source("https://raw.githubusercontent.com/fconstancias/metabaRpipe/master/Rscripts/functions.R")

readRDS("dada2/phyloseq.RDS") %>%  
phyloseq_vsearch_lulu_cluster_ASV(vsearch ="/Users/test/miniconda3/envs/metabaRpipe4/bin/vsearch") -> phyloseq_lulu_clust
	
```

322 ASV were present in the original phyloseq object:

```r
> readRDS("dada2/phyloseq.RDS") phyloseq-class experiment-level objectotu_table()   OTU Table:          [ 322 taxa and 6 samples ]:sample_data() Sample Data:        [ 6 samples by 18 sample variables ]:tax_table()   Taxonomy Table:     [ 322 taxa by 7 taxonomic ranks ]:refseq()      DNAStringSet     :      [ 322 reference sequences ]taxa are rows
```

321 `vsearch-lulu` clustering :

```r
phyloseq_lulu_clust$physeq_curatedphyloseq-class experiment-level objectotu_table()   OTU Table:          [ 321 taxa and 6 samples ]:sample_data() Sample Data:        [ 6 samples by 18 sample variables ]:tax_table()   Taxonomy Table:     [ 321 taxa by 7 taxonomic ranks ]:refseq()      DNAStringSet     :      [ 321 reference sequences ]taxa are rows
```

More information are stored under:

```r
phyloseq_lulu_clust$curated_result$otu_map
```
Check the [lulu repo](https://github.com/tobiasgf/lulu#tutorial) for more information.

<p align="right">(<a href="#top">back to top</a>)</p>

#### 5. picrust2 functional potential estimation:

Install `picrust2` in the dedicated conda environment:

```bash
conda activate metabaRpipe
conda install  -c bioconda -c conda-forge picrust2=2.3.0-b -y
```

Check the official [picrust2 repository](https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.3.0-beta)) for more details regarding the parameters.

```bash
conda activate metabaRpipe

Rscript /Users/test/Documents/GitHub/metabaRpipe/Rscripts/run_phyloseq_picrust2.Rscript --help

############################################################
Starting 
############################################################

Loading required package: optparse
Usage: /Users/test/Documents/GitHub/metabaRpipe/Rscripts/run_phyloseq_picrust2.Rscript [options]


Options:
	-p CHARACTER, --phyloseq_path=CHARACTER
		Path of the input phyloseq object

	--picrust2_bin=CHARACTER
		picrust2 script path

	--working_dir=CHARACTER
		Name of the working directory

	--output=CHARACTER
		Name of the output directory - should not exist before running the analysis

	--traits=CHARACTER
		Traits c(COG,EC,KO,PFAM,TIGRFAM)

	--min_reads=CHARACTER
		 

	--min_samples=CHARACTER
		 

	-m CHARACTER, --hsp_method =CHARACTER
		"   SP method to use."mp": predict discrete traits using
  max parsimony. "emp_prob": predict discrete traits
  based on empirical state probabilities across tips.
  "subtree_average": predict continuous traits using
  subtree averaging. "pic": predict continuous traits
  with phylogentic independent contrast. "scp":
    reconstruct continuous traits using squared-change
  parsimony "

	--no_gap_fill=CHARACTER
		Disabel gap filling 

	--add_description=CHARACTER
		add description pf Pway

	--load_picrust2_data=CHARACTER
		import picrust2 data: could generate large files... 

	--return=CHARACTER
		return picrust2 data / write to a directory

	-T NUMERIC, --slots=NUMERIC
		Number of threads to perform the analyses

	-f CHARACTER, --fun_dir=CHARACTER
		Directory containing the R functions

	-h, --help
		Show this help message and exit

```

```bash
Rscript ${MY_DIR}/metabaRpipe/Rscripts/run_phyloseq_picrust2.Rscript --phyloseq_path ~/dada2/phyloseq.RDS \
--output ~/dada2/picrust2-out \
--traits COG,EC,KO,PFAM \
--add_description TRUE \
-f ${MY_DIR}/metabaRpipe/Rscripts/functions.R
```

Check the official [picrust2 repository](https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.3.0-beta)) for more details regarding the outputs.

<p align="right">(<a href="#top">back to top</a>)</p>
## To do:

- add https://zenodo.org/account/settings/github/ -> DOI
- <s>add phylogenetic tree to a phyloseq object</s>
- <s>change name</s>
- <s>add possibility to skip primer removal: skipping run_atropos() or changing atropos parameter?</s>
- <s>replace taxonomic assignments of a phyloseq object using alternative approach/ database</s>
- <s>cluster ASV using DECIPHER</s>
- <s>cluster ASV using vsearch lulu</s>
- <s>run picrust2 from a phyloseq object</s>

 <p align="right">(<a href="#top">back to top</a>)</p>
