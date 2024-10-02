# Replication of figure generation from analysis of results from decode_pipeline

To replicate the figures in the [Tsour *et al*](https://doi.org/10.1101/2024.08.26.609665) data files must be downloaded from the [decode output](https://drive.google.com/open?id=15YoTBTZh4MdtAqHbibkYieEqyLyFi5hb&usp=drive_fs) google drive folder, specifically the supplemental data files and files in [pipeline_output/analysis_dependencies](https://drive.google.com/drive/folders/1fz9QbWfl5JxM6HdYedPqg4sc5_K4yx7j). 

Most figures can be replicated using the Jupyter notebooks in this directory. There are 5 notebooks, each corresponding to a main figure in [Tsour *et al*](https://doi.org/10.1101/2024.08.26.609665) and associated extended data and supplemental figures. The figure panels not replicated in these scripts are specified and can be replicated by running the code in [decode_analysis](https://github.com/SlavovLab/decode/tree/main/decode_analysis) or [gnomAD_analysis](https://github.com/SlavovLab/decode/tree/main/gnomAD_analysis). 

The data files and pipeline output files needed to replicate the panels in each notebook are specified in the notebook. User must edit the proj_dir variable in the second code cell to the local directory with the downloaded data. 

The script is organized by figure labels, e.g "Figure 1a" or "Extended Data Figure 2c", which can be found by searching the page. It is recommended to run specific desired chunks of the script rather than entire scripts at once. 

