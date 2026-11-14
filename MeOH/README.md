# MeOH

Design and operation of flexible electrolysis for retrofitting methanol synthesis, using day-ahead and real-time electricity prices for two locations in Texas (Panhandle and Houston).

## Contents
- `DAM_modeling.jl`, `RTM_modeling.jl`: design and operation of the facility under day-ahead (DAM) and real-time (RTM) market prices
- `H2_levelized_cost.jl`: levelized cost of hydrogen with flexible electrolysis
- `DAM_*_price.csv`, `RTM_*_price.csv`: electricity prices for each location
- `data_summary.xlsx`: physical and economic parameters

We suggest Julia 1.7.3. See also `readme.dm`.

## Reference
This code accompanies the paper:

Jiaze Ma, Michael Rebarchik, Saurabh Bhandari, Manos Mavrikakis, George W. Huber, and Victor M. Zavala. "Exploiting electricity market dynamics using flexible electrolysis units for retrofitting methanol synthesis." *Energy & Environmental Science* 16, no. 5 (2023): 2346-2357. https://doi.org/10.1039/d3ee00568b

Preprint: https://doi.org/10.26434/chemrxiv-2023-fx947

## BibTeX
```bibtex
@article{ma2023exploiting,
  title   = {{Exploiting electricity market dynamics using flexible electrolysis units for retrofitting methanol synthesis}},
  author  = {Ma, Jiaze and Rebarchik, Michael and Bhandari, Saurabh and Mavrikakis, Manos and Huber, George W. and Zavala, Victor M.},
  journal = {Energy \& Environmental Science},
  volume  = {16},
  number  = {5},
  pages   = {2346--2357},
  year    = {2023},
  doi     = {10.1039/d3ee00568b}
}
```

## Authors
Jiaze Ma, Michael Rebarchik, Saurabh Bhandari, Manos Mavrikakis, George W. Huber, and Victor M. Zavala
