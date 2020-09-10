# Analysis-of-Fate-Tracked-Montastraea-cavernosa-Colonies
 3D model generation & analysis protocols and statistical analysis pipeline
==========================================


### Ian Combs -- <icombs2017@fau.edu>
### version: May 21, 2020

------------------------------------------------------------------------
This repository contains protocols, scripts, and data associated with the 3D model generation and analysis of fate-tracked *Montastraea cavernosa* colonies within Broward County, FL.

------------------------------------------------------------------------
### Analysis and walkthroughs accompanying this repository:
[3D Model Analysis](file:///Users/iancombs/Documents/GitHub/analysisOf3dModels/3D%20Modeling/code/model.analysis.html)

------------------------------------------------------------------------

### Repository contents: 

- code/ 
  - *model.analysis.Rmd* -- analysis of 3D models Rmarkdown
  - *index.html* -- colony fate tracking statistical analysis webpage
  - *fateTracking.Rproj* -- R project file
  
- data/
  - *fateTrackingSurfaceArea.csv* -- surface Area datasheet for importing into R
  
- figures/
  - *Fig1.png* -- figure of coral colony, 3d model, and typical SCTLD manifestations
  - *Fig2.eps* -- map of study sites
  - *Fig4.png* -- box plot of tissue areas 
  - *FigS2.png* -- box plot and correlation plot of tissue loss metrics 
  - *FigS3.png* -- correlation plots of various tissue metrics
  
- protocols/ 
  - *3D_modeling_README* -- protocol for filming, rendering, and analyzing 3D Models

- scripts/
  - *ffmpeg.py* -- batch extract video files using ffmpeg (protocol found in *3D_modeling_README*)

- tables/
  - *Table2.docx* -- table of Friedman's test and pairwise comparison outputs for surface area measurements
  - *TableS1.docx* -- table of Kruskal-Wallis and pairwise comparisoon outputs for lesion count and site
  - *TableS2.docx* -- table of Kruskal-Wallis outputs for Rate of Tissue loss vs time 
  - *TableS3.docx* -- table showing Spearman's Rank Correlation of Rate of Tissue Loss and total colony size

- *gitignore* -- .gitignore
- *README.md* -- Repository readme document

