The folders here contain R scripts and data files for the article "Host phylogeny and fungal foliar endophytes co-structure C4 grass metabolomes." written by Elizabeth A. Bowman*, Colin R. Morrison*, Caroline G. Chessher, Robert M. Plowes, Lawrence E. Gilbert

Elizabeth A. Bowman is the author for correspondence on this article. 
Address: 2907 Lake Austin Blvd.
Brackenridge Field Laboratory
University of Texas at Austin
Austin, TX 78703 USA. 
Email: eabowman@utexas.edu

Elizabeth A. Bowman and Colin R. Morrison wrote the scripts for data analysis. Authorship is indicated below.

------------------------------------------------------------------------------------------

Explanation of folders:
1. The data folder contains all data files used in the R script for analyses with
explanations of columns.
2. The figures folder is an output folder where figures generated in the R script will be
output.
3. The results folder is an output folder for results tables generated in the R script.
4. The script folder contains all scripts organized by type of analysis.

The data folder contains a 'Read me' file with column and file descriptions.

Explanation of each script is below.
------------------------------------------------------------------------------------------

Loadlibraries.R: Loads all libraries used in analyses. The installation of each library is commented out.

analysis.GrassPlantPhylogeny.R: Creates Supplementary Fig. S1. Written by E.A. Bowman

analysis.Q1_Metabolome.R: Code for assessing metabolomic variation between grass species; Panel A for Figure 1. Written by Colin R. Morrison. 

analysis.Q1_FE.R: Code for assessing culture-free (CF) and culture-based (CB) variation between grass species. Also assesses correlation between metabolome and endophytic communities; Panel B and C for Figure 1. Written by E.A. Bowman. 

analysis.Q2.R: Code for assessing correlation of the metabolome, CF and CB endophyte community, and host relatedness; Figure 2. Written by E.A. Bowman. 

analysis.Q3.R: Code for LASSO analysis of metabolites significantly associated with FE presence; Figure 3. Written by Colin R. Morrison.

CoverageBasedRarefaction.R: Code for doing coverage based rarefaction of CF Illumina data. Originally written by Komei Kadowaki, modified by E.A. Bowman.

varpart2.MEM.R: Function for distance based Moran's Eigenvector (dbMEM) for geographical distance. Written by Pierre Legendre. 

metabolomics_processing_pipeline.R: Code for pipeline to process UHPLC-MS/MS metabolomomic feature data, subset metabolomes for chemotaxonomic groups of interest. Written by Colin R. Morrison.
