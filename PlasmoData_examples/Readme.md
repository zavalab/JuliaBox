# PlasmoData.jl Examples

This repository contains code for the examples given in the manuscript "PlasmoData.jl -- A Julia Framework for Modeling and Analyzing Complex Data as Graphs". The source code for the packages highlighted in this work are available here: [PlasmoData.jl](https://github.com/zavalab/PlasmoData.jl) and [DataGraphPlots.jl](https://github.com/zavalab/DataGraphPlots.jl/).

## Structure
The repository is structured as follows. The subdirectories include:
 * `LC_data/` - Data for liquid crystal analysis from the work of [Jiang et al., 2023](https://doi.org/10.1021/acs.jpcc.3c03076) and [Bao et al., 2023](https://doi.org/10.1021/jacs.2c03424). The data has been transferred into a .jld file in this folder for convenience. Note that this data must be unpacked from the zip file `LC_data.zip`.
 * `pathway_example/` - data for a technology pathway example
 * `Recife_data/` - data for dengue outbreaks in Recife Brazil from [de Souza et al., 2022](https://doi.org/10.1088/1742-5468/aca0e5)
 * `supporting_info_scripts/` - scripts for performing the analysis in the supplementary information
 * `tutorials/` - scripts for general PlasmoData.jl modeling tasks discussed in the paper

The following scripts are also available in the main repository and correspond to the three examples.

 * `LC_SVMs.jl` - script for analyzing the liquid crystal data using PlasmoData.jl for feature extraction and classifying the problems using SVMs.
 * `LC_on_GNNs.jl` - script for analyzing the liquid crystal data using PlasmoData.jl for feature extraction and classifying the problems using GNNs.
 * `dengue_analysis.jl` - script for analyzing the dengue data.
 * `pathway_analysis_metrics.jl` - script for performing the pathway analysis.


## Versions

The analyses above were performed with Julia 1.7.3. Package versions used in this work can be accessed through the Project.toml and Manifest.toml files included in this directory. The specific commits used for PlasmoData.jl and DataGraphPlots.jl were 29122bc and ebe7fee, respectively, and these must be added separately from the Manifest.toml files. 