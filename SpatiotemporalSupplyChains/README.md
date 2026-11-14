# SpatiotemporalSupplyChains

Spatio-temporal economic properties of multi-product supply chains with storage, with Wisconsin case studies.

## How to run
Run `STMarketsModelV3.jl` in any directory. Under the `### SET DIRECTORY ###` heading in that file, choose the case study to run (uncomment the desired case study and comment out the rest).

`RecordMaker.jl` saves problem data to .csv files for easy access and plotting. Records are saved within each case study folder inside `TestingFiles/`.

## Contents
- `STMarketsModelV3.jl`: spatio-temporal supply chain model
- `RecordMaker.jl`: saves results to .csv files
- `TestingFiles/`: case studies (e.g., `WisconsinV3/`, `WisconsinV3NoStorage/`)

## Reference
This code accompanies the paper:

Philip A. Tominac, Weiqi Zhang, and Victor M. Zavala. "Spatio-temporal economic properties of multi-product supply chains." *Computers & Chemical Engineering* 159 (2022): 107666. https://doi.org/10.1016/j.compchemeng.2022.107666

Preprint: https://arxiv.org/abs/2106.13836

## BibTeX
```bibtex
@article{tominac2022spatio,
  title   = {{Spatio-temporal economic properties of multi-product supply chains}},
  author  = {Tominac, Philip A. and Zhang, Weiqi and Zavala, Victor M.},
  journal = {Computers \& Chemical Engineering},
  volume  = {159},
  pages   = {107666},
  year    = {2022},
  doi     = {10.1016/j.compchemeng.2022.107666}
}
```

## Authors
Philip A. Tominac, Weiqi Zhang, and Victor M. Zavala
