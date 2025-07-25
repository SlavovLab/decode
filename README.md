# Decode: Computational pipeline for quantifying amino acid substitutions from alternate RNA decoding



* [Decode Website](https://decode.slavovlab.net) &nbsp; | &nbsp; [Download data](https://decode.slavovlab.net/mass-spec/data)
* [Preprint article](https://doi.org/10.1101/2024.08.26.609665)


The code in this repository can be used to identify, validate, and quantify amino acid substitutions in LC-MS proteomics data that arise from alternate decoding of RNA.

The [decode_pipeline](https://github.com/SlavovLab/decode/tree/main/decode_pipeline) has several steps and multiple sections of the pipeline require input from external analysis, such as database searches of proteomics data. As such, this is not a standalone software package that can be run automatically. Detailed instructions for running the pipeline are provided in the various README.md files in this repository.

Results and output files from running this pipeline on datasets as described by [Tsour *et al*](https://doi.org/10.1101/2024.08.26.609665) can be accessed here: [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs). These files can be downloaded and used to replicate figures using the code provided in [decode_figures](https://github.com/SlavovLab/decode/tree/main/decode_figures) and [decode_analysis](https://github.com/SlavovLab/decode/tree/main/decode_analysis).

## Overview of pipeline
1. for protein sequence database generation by *in-silico* translation of RNA-seq data ([custom_protein_database_pipeline](https://github.com/SlavovLab/decode/tree/main/custom_protein_database_pipeline))
2. Dependent peptide search with MaxQuant (external software analysis)
3. Identifying candidate peptides with amino acid substitutions  ([decode_pipeline](https://github.com/SlavovLab/decode/tree/main/decode_pipeline))
4. Validation database search (external software analysis)
5. Quantifying validated peptides with amino acid substitutions ([decode_pipeline](https://github.com/SlavovLab/decode/tree/main/decode_pipeline))
6. Downstream data analysis ([decode_figures](https://github.com/SlavovLab/decode/tree/main/decode_figures), [decode_analysis](https://github.com/SlavovLab/decode/tree/main/decode_analysis))

## Usage
### Step 1: Custom protein databases
Use RNA-seq data matched to LC-MS proteomics data to create sample-specific protein databases.

The code for this step is in [custom_protein_database_pipeline](https://github.com/SlavovLab/decode/tree/main/custom_protein_database_pipeline) and the [README.md](custom_protein_database_pipeline/README.md) in that directory contains detailed instructions for running the code.

If no matched RNA-seq data is available, this step can be skipped, but caution should be taken in interpreting quantified amino acid substitutions as there is lower confidence that they are not encoded in the genome.

### Step 2: Identifying modified peptides with MaxQuant
The dependent peptide search algorithm in MaxQuant is used to identify peptides with modifications in LC-MS proteomics data.

The LC-MS proteomics data is ideally searched against the sample-specific database generated in Step 1. If not available, species-specific [UniProt](https://www.uniprot.org/) fasta can be used.

A sample MaxQuant parameter file is provided in [MaxQuant_templates](https://github.com/SlavovLab/decode/tree/main/MaxQuant_templates), along with a script to create a new parameter file with user-defined parameters (raw files, fasta, etc.)

The output from this dependent peptide search is required to proceed with the next steps of the pipeline.

### Step 3: Identifying candidate alternate translation events
Search for modified peptides in dependent peptide search results that may represent amino acid substitutions. Add candidate peptides to custom protein sequence databases for validation search.

The code for this step can be found in [decode_pipeline/python_scripts](https://github.com/SlavovLab/decode/tree/main/decode_pipeline/python_scripts). [decode_pipeline/README.md](https://github.com/SlavovLab/decode/tree/main/decode_pipeline/README.md) contains detailed instructions for running this code.

### Step 4. Validation search with MaxQuant (or another proteomics data search engine)
Run a standard database search using the protein databases appended with candidate substituted peptides (step 3).

A sample MaxQuant parameter file is provided in [MaxQuant_templates](https://github.com/SlavovLab/decode/tree/main/MaxQuant_templates). The output from this validation search is required to proceed with the next steps of the pipeline.

### Step 5. Quantify alternate decoding events
The code for this step can be found in [decode_pipeline/python_scripts](https://github.com/SlavovLab/decode/tree/main/decode_pipeline/python_scripts). [decode_pipeline/README.md](decode_pipeline/README.md) contains detailed instructions for running this code.

### Step 6. Downstream data analysis
The data generated in the pipeline outlined above can be further analyzed in many ways. The code for the analysis described in [Preprint article](https://doi.org/10.1101/2024.08.26.609665) is provided in [decode_analysis](https://github.com/SlavovLab/decode/tree/main/decode_analysis) and [decode_figures](https://github.com/SlavovLab/decode/tree/main/decode_figures) and [gnomAD_analysis](https://github.com/SlavovLab/decode/tree/main/gnomAD_analysis).

These directories contain different subsets of the analyses described in the paper along with the code to reproduce the figures. Each contain a README.md describing the analyses contained.

The downstream analyses are dependent on the many data files generated in the pipeline described above. In order to keep this analysis as reproducible as possible, we have deposited all our relevant output files to [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs). If this data is downloaded and the "proj_dir" parameter set to the download location, all the figures in [decode_figures](https://github.com/SlavovLab/decode/tree/main/decode_figures) can be generated.

## License
The code is avilable under **[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)** license 
