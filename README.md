# PaIRKAT (Shiny)
R Shiny app interface for PaIRKAT algorithm

A hosted instance of the application is available through Binder. 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Ghoshlab/PaIRKAT_Shiny/main?urlpath=shiny/app/)

PaIRKAT is model framework for assessing statistical relationships between networks of metabolites (pathways) and clinical outcome. PaIRKAT queries the KEGG database to determine interactions between metabolites from which network connectivity is constructed.

# Use
## Upload Data
Three datasets are required to use PaIRKAT

1. Phenotype Data: Contains outcomes of interest and any meaningful covariates to be adjusted for. Rows should be subjects and columns should be variables.
2. Metabolite Measurements: Contains measurments of metabolites for all subjects. Rows should be subjects and columns should be metabolite names. One column should have subject IDs matching subject IDs in clinical data.
3. Pathway Data: Contains linkage data and pathway information. Rows are metabolites and columns are variables. Should contain a column with KEGG IDs.

## Example Data
Example data to test PaIRKAT are available from `exampleData/Example_PaIRKAT_Data.zip`. This zipped folder contains 3 files which correspond to the 3 datasets detailed above. 

## Gather Pathways
The Gather Pathways tab will guide you though defining data linkage logic and querying the KEGG database to collect pathway information and form networks of metabolites.

## Run PaIRKAT
Define outcome of interest and clinical covaraites to control for. Output from this analysis is a list of significant pathways associated with the outcome of interest.

## Explore Results
Visual tools to explore the results of PaIRKAT analysis. There are two tools provided with PaIRKAT.

1. Network Graph: This tool allows you to view significant pathways and their connectivity both within and between pathways. Many options are available to color and size nodes to better emphasize features of the network.
2. Plot Builder: This flexible tool allows you to make plots using any of the information entered into the application or derived from its functions.
    
# Methods

PaIRKAT is a tool for improving testing power on high dimensional data by including graph topography in the kernel machine regression setting. Studies on high dimensional data can struggle to include the complex relationships between variables. The semi-parametric kernel machine regression model is a powerful tool for capturing these types of relationships. They provide a framework for testing for relationships between outcomes of interest and high dimensional data such as metabolomic, genomic, or proteomic pathways. We propose PaIRKAT, a method for including known biological connections between high dimensional variables into the kernel machine by representing them as edges of 'graphs' or 'networks.' It is common for nodes (e.g. metabolites) to be disconnected from all others within the graph, which leads to meaningful decreases in testing power whether or not the graph information is included. We include a graph regularization or 'smoothing' approach for managing this issue.

# Citation

PaIRKAT: A pathway integrated regression-based kernel association test with applications to metabolomics and COPD phenotypes 
Charlie M. Carpenter, Weiming Zhang, Lucas Gillenwater, Cameron Severn, Tusharkanti Ghosh, Russel Bowler, Katerina Kechris, Debashis Ghosh 
bioRxiv 2021.04.23.440821; doi: https://doi.org/10.1101/2021.04.23.440821

# Acknowledgements
The project described was supported by Award Number U01 HL089897 and Award Number U01 HL089856 from the National Heart, Lung, and Blood Institute, and U01 CA235488 from the National Cancer Institute. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Heart, Lung, and Blood Institute, National Cancer Institute, or the National Institutes of Health.

# License
PaIRKAT is released under the GNU General Public License version 3 (GPLv3)