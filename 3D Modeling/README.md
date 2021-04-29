# Analysis-of-Fate-Tracked-Montastraea-cavernosa-Colonies
 3D model generation & analysis protocols and statistical analysis pipeline
==========================================


### Ian Combs -- <combsi9892@gmail.com>
### version: April 29, 2021

------------------------------------------------------------------------
This repository contains protocols, scripts, and data associated with the 3D model generation and analysis of fate-tracked *Montastraea cavernosa* colonies within Broward County, FL. Found in **Combs IR, Studivan MS, Eckert RJ, Voss JD. 2021. Quantifying the impacts of stony coral tissue loss disease on corals in Southeast Florida through surveys and 3D photogrammetry**

------------------------------------------------------------------------
### Analysis and walkthroughs accompanying this repository:
[3D Model Analysis](https://icombs2017.github.io/analysisOf3dModels/3D%20Modeling/code/)

------------------------------------------------------------------------

### Repository contents:

- code/
  - *model.analysis.Rmd* -- analysis of 3D models Rmarkdown
  - *index.html* -- colony fate tracking statistical analysis webpage
  - *fateTracking.Rproj* -- R project file

- data/
  - *S2_Dataset.csv* -- surface Area datasheet for importing into R

- figures/
  - *Fig1.png* -- figure of coral colony, 3d model, and typical SCTLD manifestations
  - *Fig2.eps* -- map of study sites
  - *Fig2.png* -- map of study sites
  - *Fig4.png* -- box plot of tissue areas through time
  - *Fig4.eps* -- box plot of tissue areas through time
  - *Fig5.png* -- box plot and correlation plot of tissue loss metrics
  - *Fig5.eps* -- box plot and correlation plot of tissue loss metrics
  - *Fig6.png* -- Correlation plot of Disease lesion area and total colony area
  -*Fig6.eps* -- Correlation plot of Disease lesion area and total colony area

- protocols/
  - *3D_modeling_README* -- protocol for filming, rendering, and analyzing 3D Models

- scripts/
  - *ffmpeg.py* -- batch extract video files using ffmpeg (protocol found in *3D_modeling_README*)

- tables/
  - *Table2.docx* -- table of Friedman's test and pairwise comparison outputs for surface area measurements
  - *TableS1.docx* -- table of Kruskal-Wallis and pairwise comparison outputs for lesion count and site


- *gitignore* -- .gitignore
- *README.md* -- Repository readme document
