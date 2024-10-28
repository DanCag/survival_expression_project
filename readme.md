# Survival outcome - expression project
This repository contains all scripts to reproduce analyses, tables and figures 
of survival-expression project, where we investigate the associations between 
gene expression and survival (see thesis for more details)

## Structure
Each folder refers to a figure or table of 3rd chapter of the thesis
where you find all necessary scripts to reproduce that figure or table.<br>

`functions` folder contains the functions called in different scripts.<br>

`01_data_processing` folder contains the script to extract TCGA expression
and survival info to build survival-expression tables. <br>

Repository contains only scripts, while data is stored in the ws.
If the script needs to be run in the workstation (old or new), you will find 
a comment at the beginning of the script.<br>

## Run
Open the `R.proj` file and run the script you are interested in. 
