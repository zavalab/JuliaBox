# ElectrolysisProtocols

Flexible versus structured operating protocols for scheduling electrolyzer stacks under real-time electricity prices (ERCOT, July 2024).

## Contents
- `data.jl`: loads the data; `ercot_rtm_July2024.xlsx`: ERCOT real-time market prices for July 2024
- `func_flex_scheduling.jl`: flexible scheduling of the stacks
- `func_strc_scheduling.jl`: scheduling with structured operating protocols
- `run_flex_strc_day.ipynb`, `run_flex_strc_week.ipynb`: compare the two over a day and over a week

## Reference
This code accompanies the paper:

Jukbin Kim, Pen-Yeh Yen, Abiela M. Bradley, Rachel N. Gaines, Paul J. A. Kenis, and Victor M. Zavala. "Achieving Flexible Operation by Aggregating Structured Operating Protocols: An Electrolysis Case Study Using Real Market Data." *Industrial & Engineering Chemistry Research* 65, no. 40 (2026): 21675-21696. https://doi.org/10.1021/acs.iecr.6c01628

Preprint: https://doi.org/10.26434/chemrxiv.15001860/v1

## BibTeX
```bibtex
@article{kim2026achieving,
  title   = {{Achieving Flexible Operation by Aggregating Structured Operating Protocols: An Electrolysis Case Study Using Real Market Data}},
  author  = {Kim, Jukbin and Yen, Pen-Yeh and Bradley, Abiela M. and Gaines, Rachel N. and Kenis, Paul J. A. and Zavala, Victor M.},
  journal = {Industrial \& Engineering Chemistry Research},
  volume  = {65},
  number  = {40},
  pages   = {21675--21696},
  year    = {2026},
  doi     = {10.1021/acs.iecr.6c01628}
}
```

## Authors
Jukbin Kim, Pen-Yeh Yen, Abiela M. Bradley, Rachel N. Gaines, Paul J. A. Kenis, and Victor M. Zavala
