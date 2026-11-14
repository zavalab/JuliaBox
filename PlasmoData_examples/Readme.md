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

## References
This code accompanies the papers:

1. Nanqi Bao, Jake I. Gold, Jonathan K. Sheavly, James J. Schauer, Victor M. Zavala, Reid C. Van Lehn, Manos Mavrikakis, and Nicholas L. Abbott. "Ordering Transitions of Liquid Crystals Triggered by Metal Oxide-catalyzed Reactions of Sulfur Oxide Species." *Journal of the American Chemical Society* 144, no. 36 (2022): 16378-16388. https://doi.org/10.1021/jacs.2c03424

2. David L. Cole and Victor M. Zavala. "PlasmoData.jl — A Julia framework for modeling and analyzing complex data as graphs." *Computers & Chemical Engineering* 185 (2024): 108679. https://doi.org/10.1016/j.compchemeng.2024.108679

Preprint: https://arxiv.org/abs/2401.11404

3. Shengli Jiang, Nanqi Bao, Alexander D. Smith, Shraddha Byndoor, Reid C. Van Lehn, Manos Mavrikakis, Nicholas L. Abbott, and Victor M. Zavala. "Scalable Extraction of Information from Spatiotemporal Patterns of Chemoresponsive Liquid Crystals Using Topological Descriptors." *The Journal of Physical Chemistry C* 127, no. 32 (2023): 16081-16098. https://doi.org/10.1021/acs.jpcc.3c03076

Preprint: https://doi.org/10.26434/chemrxiv-2023-fjcxd

## BibTeX
```bibtex
@article{bao2022ordering,
  title   = {{Ordering Transitions of Liquid Crystals Triggered by Metal Oxide-catalyzed Reactions of Sulfur Oxide Species}},
  author  = {Bao, Nanqi and Gold, Jake I. and Sheavly, Jonathan K. and Schauer, James J. and Zavala, Victor M. and Van Lehn, Reid C. and Mavrikakis, Manos and Abbott, Nicholas L.},
  journal = {Journal of the American Chemical Society},
  volume  = {144},
  number  = {36},
  pages   = {16378--16388},
  year    = {2022},
  doi     = {10.1021/jacs.2c03424}
}
```

```bibtex
@article{cole2024plasmodata,
  title   = {{PlasmoData.jl — A Julia framework for modeling and analyzing complex data as graphs}},
  author  = {Cole, David L. and Zavala, Victor M.},
  journal = {Computers \& Chemical Engineering},
  volume  = {185},
  pages   = {108679},
  year    = {2024},
  doi     = {10.1016/j.compchemeng.2024.108679}
}
```

```bibtex
@article{jiang2023scalable,
  title   = {{Scalable Extraction of Information from Spatiotemporal Patterns of Chemoresponsive Liquid Crystals Using Topological Descriptors}},
  author  = {Jiang, Shengli and Bao, Nanqi and Smith, Alexander D. and Byndoor, Shraddha and Van Lehn, Reid C. and Mavrikakis, Manos and Abbott, Nicholas L. and Zavala, Victor M.},
  journal = {The Journal of Physical Chemistry C},
  volume  = {127},
  number  = {32},
  pages   = {16081--16098},
  year    = {2023},
  doi     = {10.1021/acs.jpcc.3c03076}
}
```
