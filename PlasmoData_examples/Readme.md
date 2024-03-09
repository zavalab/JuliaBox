# PlasmoData.jl Examples

This repository contains code for the examples given in the manuscript ["PlasmoData.jl -- A Julia Framework for Modeling and Analyzing Complex Data as Graphs"](https://arxiv.org/abs/2401.11404). The source code for the packages highlighted in this work are available here: [PlasmoData.jl](https://github.com/zavalab/PlasmoData.jl) and [PlasmoDataPlots.jl](https://github.com/zavalab/PlasmoDataPlots.jl/).

## Structure
The repository is structured as follows. The subdirectories include:
  * `tutorials/` - scripts for general PlasmoData.jl modeling tasks discussed in the paper. These generally follow the code snippets given in Sections 2 and 3.
 * `supporting_info_scripts/` - scripts for performing the analysis in the supplementary information
 * `CS1-Image_Analysis` - Scripts and data for performing the analysis in case study 1 of the manuscript on image analysis and feature extraction. The liquid crystal data comes from the work of [Jiang et al., 2023](https://doi.org/10.1021/acs.jpcc.3c03076) and [Bao et al., 2023](https://doi.org/10.1021/jacs.2c03424). The data has been transferred into a .jld file in this folder for convenience. When cloning the git repo, this data must be unpacked from zip file, `LC_data.zip`.
 * `CS2-Multivariate_Time_Series_Analysis` - Script and data for performing the analysis in case study 2 of the manuscript on representing multivariate time series with edge weighted graphs. The dengue data for Recife Brazil comes from the work of [de Souza et al., 2022](https://doi.org/10.1088/1742-5468/aca0e5).
 * `CS3-Connectivity_Analysis` - Data and scripts for performing the connectivity/pathway analysis in case study 3 of the manuscript.

## Versions

The analyses above were performed with Julia 1.7.3. Package versions used in this work can be accessed through the Project.toml and Manifest.toml files included in this directory. 

The environment can be activated by navigating to this directory, typing `]` in the Julia REPL, and then typing `activate .`.