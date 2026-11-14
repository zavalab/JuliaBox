# FERC2222Cases

Market model case studies of manure-to-biogas electricity generation on dairy farms before and after FERC Order 2222.

## Contents
- `MarketsModel.jl`: market model (requires JuMP and Clp)
- `Pre 174-83/`, `Post 171-69/`, `Post 174-83/`, `Post 363-71/`, `Post 861-71/`: input data for each case study

To run a case, set the `DataFolder` variable in `MarketsModel.jl` to the case folder and run the script (see `DIRECTIONS.rtf`). Results are saved to `_SolutionData.txt`.

## Reference
This code accompanies the paper:

Evan D. Erickson, Philip A. Tominac, and Victor M. Zavala. "Biogas production in United States dairy farms incentivized by electricity policy changes." *Nature Sustainability* 6, no. 4 (2023): 438-446. https://doi.org/10.1038/s41893-022-01038-9

## BibTeX
```bibtex
@article{erickson2023biogas,
  title   = {{Biogas production in United States dairy farms incentivized by electricity policy changes}},
  author  = {Erickson, Evan D. and Tominac, Philip A. and Zavala, Victor M.},
  journal = {Nature Sustainability},
  volume  = {6},
  number  = {4},
  pages   = {438--446},
  year    = {2023},
  doi     = {10.1038/s41893-022-01038-9}
}
```

## Authors
Evan D. Erickson, Philip A. Tominac, and Victor M. Zavala
