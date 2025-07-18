---
layout: default
title: Download data
nav_order: 2
permalink: mass-spec/data
description: "Data repository for results from analyzing alternate RNA decoding, intermediate pipeline outputs and results"
nav_exclude: false
---
{% include social-media-links.html %}

# Download alternate RNA decoding results

&nbsp;

&nbsp;

[decode data]({{site.baseurl}}#plexDIA-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
<!--- [plexDIA @ massIVE]({{site.baseurl}}#RAW-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 } --->
[decode code](https://github.com/SlavovLab/decode){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }

<h2 style="letter-spacing: 2px; font-size: 26px;" id="plexDIA-data" >Decode pipeline results organized by datasets</h2>
All results from processing LC-MS proteomics data through [decode pipeline](https://github.com/SlavovLab/decode) described in [Tsour *et al*](https://doi.org/10.1101/2024.08.26.609665) are organized in this [directory](https://drive.google.com/drive/folders/15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb?usp=sharing).

[Description](https://docs.google.com/document/d/17Bpu_kIfnSnGpETMQWQM7W9PUvMrpAe9/edit?usp=drive_link&ouid=109814487119977139380&rtpof=true&sd=true) of files in [directory](https://drive.google.com/drive/u/3/folders/15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb).

## Supplemental Data Files
* [Supplemental Data 1](https://drive.google.com/file/d/1klr3dtI5FcoV-eyTc696T5y0v51WA2n9/view?usp=sharing): A table of peptides identified through the dependent peptide search as having known post-translational or chemical modifications. This table includes peptides from all 6 CPTAC datasets and the label-free healthy human tissue dataset. Types and locations of modifications are specified for each peptide.
* [Supplemental Data 2](https://docs.google.com/spreadsheets/d/1_9EACBIJWJpJqy5eEyJ2bbrnfMHy5clg/edit?usp=sharing&ouid=100652957809040104900&rtpof=true&sd=true): A table containing 1 row per unique SAAP-BP pair per dataset (CPTAC and label-free). Each SAAP-BP pair is listed along with RAAS summary stats across the dataset and protein mapping information, including Pfam domains. The first tab of the file provides a description of the columns in the data table.
* [Supplemental Data 3](https://docs.google.com/spreadsheets/d/1r6vyessLaT9XHgyNYmJ9hTdWQ07di1wr/edit?usp=sharing&ouid=100652957809040104900&rtpof=true&sd=true): A table containing 1 row per unique SAAP-BP pair per TMT set (CPTAC) or healthy tissue type (label-free). Each SAAP-BP pair is listed along with their precursor ion level abundances and RAAS values computed from precursor ion intensities. The first tab of the file provides a description of the columns in the data table.
* [Supplemental Data 4](https://docs.google.com/spreadsheets/d/1qDePyzVOAq_JwVMMgqGaQkkbncee5VpQ/edit?usp=sharing&ouid=109814487119977139380&rtpof=true&sd=true): A table containing 1 row per unique SAAP-BP pair per patient sample (CPTAC only). Each SAAP-BP pair is listed along with their reporter ion level abundances and RAAS values computed from reporter ion intensities. The first tab of the file provides a description of the columns in the data table.
* [Supplemental Data 5](https://docs.google.com/spreadsheets/d/1vrVpevuzpacA73gyqZkBJSudse0z0Nv8/edit?usp=sharing&ouid=109814487119977139380&rtpof=true&sd=true): Table of SAAP-BP pairs identified in SILAC-labeled liver cells containing one row per unique SAAP-BP pair per sample. Precursor ion level abundances of the light and heavy peptides in each sample are listed, along with their sums, ratios and RAAS values. The first tab of the file provides a description of the columns in the data table.
* [Supplemental Data 6](https://docs.google.com/spreadsheets/d/18QHi-Jnnl26uiOPUx9HnJON49QMNRYAs/edit?usp=sharing&ouid=109814487119977139380&rtpof=true&sd=true): A subset of Supplemental_Data_2.SAAP_proteins.xlsx containing SAAP, BP and RAAS data pertaining to proteins that have been implicated in neurodegeneration and dementia (Uniprot).
* [Supplemental Data 7](https://drive.google.com/file/d/1XaiRGnfoTpMACKSmBPaFomcZ0lZOifbx/view?usp=sharing): Mapping of the amino acid substitution sites of unique pairs of BP/SAAP to proteins, their transcripts and genome coordinates, as defined in the Ensembl genome release *GRCh38.110*.
* [Supplemental Data 8](https://docs.google.com/spreadsheets/d/1uZPENtn3H-LIlfTeVPwnZvv56IJJ7zw3/edit?usp=sharing&ouid=100652957809040104900&rtpof=true&sd=true): Supplemental_Data_3 filtered for SAAP with RAAS > 0.1 and positional probability >0.9.
* [Supplemental Data 9](https://docs.google.com/spreadsheets/d/15Rrh-odEScRtvezZl7LbXlRUQJCwSqU6/edit?usp=sharing&ouid=100652957809040104900&rtpof=true&sd=true): Data tables with gene locus, allele frequency and constraint data for all observed substitutions.
* [Supplemental Data 10](https://drive.google.com/file/d/1bUkBZ3kJFBIv9OU7SMeyHGQEB8cuKllq/view?usp=sharing): Aligned reads to the transcripts coding for the base peptides (and corresponding SAAP) shown in Supplemental Figure 1 and Supplemental Figure 2.


&nbsp;

<!---
<h2 style="letter-spacing: 2px; font-size: 26px;" id="RAW-data" >plexDIA RAW data and search results from DIA-NN</h2>
The repositories below contain RAW mass-spectrometry data files generated by a first-generation Q-exactive instrument as well as the search results from analyzing the  RAW files by [DIA-NN](https://drive.google.com/file/d/1naoAhDX6VyvQ8Uc1ukfpcMcKzyTFbDCv/view?usp=sharing). Searching plexDIA data with DIA-NN is described in this [tutorial](https://youtu.be/0Wmg9LjDtgE).


* **MassIVE Repository for version 1 (Bulk plexDIA data):**
  - [**http:**  MSV000088302](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=8b0a2f5b2fc84964b4bd4ee64fc84d25)
  - [**ftp:** &nbsp; MSV000088302](ftp://massive.ucsd.edu/MSV000088302)

* **MassIVE Repository for version 2 (Bulk and single-cell plexDIA data):**
    - [**http:**  MSV000089093](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ae918c7ce5a94a4abd2c6b54a3806c9e)
    - [**ftp:** &nbsp; MSV000089093](ftp://massive.ucsd.edu/MSV000089093)


[plexDIA_Article]: https://doi.org/10.1101/2021.11.03.467007 "Increasing the throughput of sensitive proteomics by multiplexed data-independent acquisition using plexDIA"
--->

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

&nbsp;
