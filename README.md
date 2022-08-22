# FMT_Meta_Machine_Learning

The code in this repository reproduces display items of the Machine Learning section of the Fecal Microbiota Transplantation Meta-Analysis conducted by Ianiro, Punchochar, Karcher et al. (2022, Nature Medicine). Specifically, it reproduces main text figure 4 and accompanying supplementary figures.

### Prepare

In order to set up your environment to be able to run this code, start by installing conda (as described under this link https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Then create a fresh conda environment using the `.yml` file in this repository like so

conda env create -f FMT.yml

This will create a new environment called `FMT` which you can activate by typing `conda activate FMT` after the environment was succesfully set up.

### Reproduce

Reproduce the display items by typing `Rscript FMT.r`. It will take around a day to finish.
