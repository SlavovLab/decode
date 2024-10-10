---
layout: default
title: Data Analysis
nav_order: 3
permalink: mass-spec/decode_analysis
description: "Quantifying alternate RNA decoding: Data analysis pipeline"
nav_exclude: false
---
{% include social-media-links.html %}

# Quantifying alternate RNA decoding
{:.no_toc}

&nbsp;

{: .fs-5 .fw-300}
This section organizes instructions and links to software and pipelines for identifying, validating, and quantifying amino acid substitutions from alternate RNA decoding.

* Will be replaced with the ToC, excluding the section header
{:toc}

&nbsp;

## Creating sample-specific protein databases from RNA-seq data


### Pipeline for generating custom protein databases

This [pipeline](https://github.com/SlavovLab/decode/tree/main/custom_protein_database_pipeline) is used for generating sample-specific protein fasta databases from RNA-seq data.

-------


&nbsp;


## Identifying, validating, and quantifying amino acid substitutions from alternate RNA decoding

### Find modified peptides in LC-MS proteomics data with MaxQuant dependent peptide search
Search LC-MS proteomics data with [MaxQuant](https://www.maxquant.org/) dependent peptide search against sample-specific protein databases to identify peptides with modifications which may represent amino acid substitutions.

Sample [MaxQuant parameter file](https://github.com/SlavovLab/decode/tree/main/MaxQuant_templates) is provided in the [decode GitHub repository](https://github.com/SlavovLab/decode).

A [useful tutorial](https://atchen.me/research/2019/03/21/mq-linux.html) for running MaxQuant searches in Linux.

### Identify candidate <u>s</u>ubstituted <u>a</u>mino <u>a</u>cid <u>p</u>eptides (SAAP)
Apply [decode pipeline](https://github.com/SlavovLab/decode/decode_pipeline) steps 1 and 2 to identify peptides with potential amino acid substitutions from alternate RNA decoding. See [detailed pipeline instructions](https://github.com/SlavovLab/decode/tree/main/decode_pipeline#readme).

### Validate SAAP
[Decode pipeline](https://github.com/SlavovLab/decode/decode_pipeline) step 3.
Search LC-MS proteomics data with [MaxQuant](https://www.maxquant.org/) (or other search engine) standard database search against sample-specific protein databases appended with candidate SAAP.

Sample [MaxQuant parameter file](https://github.com/SlavovLab/decode/tree/main/MaxQuant_templates) is provided in the [decode GitHub repository](https://github.com/SlavovLab/decode).

### Quantify SAAP
Apply [decode pipeline](https://github.com/SlavovLab/decode/decode_pipeline) steps 4 and 5 to quantify validated amino acid substitutions from alternate RNA decoding. See [detailed pipeline instructions](https://github.com/SlavovLab/decode/tree/main/decode_pipeline#readme).


Results from running [decode pipeline](https://github.com/SlavovLab/decode/decode_pipeline) on TMT-labeled and label-free datasets, as described by [Tsour *et al*][Decode_article] are provided in [decode data directory](https://drive.google.com/drive/u/3/folders/15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb).

-------


&nbsp;


[Detailed methods](https://www.biorxiv.org/content/10.1101/2024.08.26.609665v1.full#:~:text=confirmed%20this%20result.-,Methods,-Sample%2Dspecific%20protein) are found in the article, [Tsour *et al*][Decode_article].


<!--- [plexDIA_Article]: https://doi.org/10.1101/2021.11.03.467007 "Multiplexed data-independent acquisition by plexDIA"
[plexDIA_Nature]: https://doi.org/10.1038/s41587-022-01389-w "Derks, J., Slavov, N. et al. Increasing the throughput of sensitive proteomics by plexDIA. Nat Biotechnol (2022)"--->
[decode_Code]: https://github.com/SlavovLab/decode "Decode data analysis pipeline, GitHub repository from the Slavov Laboratory"

-------


&nbsp;


## Data pipelines for reproducing decode data analysis
A [pipeline][decode_Code] for reproducing the analysis by [Tsour *et al*][Decode_article] is available at the [decode GitHub repository][decode_Code].  


* [Pipeline for processing RNA-seq data through custom protein database pipeline @ GitHub](https://github.com/SlavovLab/decode/custom_protein_database_pipeline)
* [Pipeline for processing LC-MS proteomics data through decode pipeline @ GitHub](https://github.com/SlavovLab/decode/decode_pipeline)

-------


[Decode_article]: https://www.biorxiv.org/content/10.1101/2024.08.26.609665v2.full "Alternate RNA decoding results in stable and abundant proteins in mammals"

&nbsp;  

&nbsp;

&nbsp;  

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;
