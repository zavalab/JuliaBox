# IPCluster_examples

Julia interface and stochastic farmer examples for the clustering-based interior-point solver ClusterIPM.

## Contents
- `ClusterIPM.jl`: Julia interface to ClusterIPM (Nai-Yuan Chiang and Kibaek Kim, 2015)
- `farm_IPCluster.jl`, `farm_IPCluster.ipynb`, `farm_IPCluster_parallel.ipynb`: stochastic farmer problem, serial and parallel
- `farmlarge_IPCluster.jl`, `farmlarge_IPCluster.ipynb`, `farmlarge_IPCluster_parallel.ipynb`: a larger version of the farmer problem

## Reference
This code accompanies the paper:

Yankai Cao, Carl D. Laird, and Victor M. Zavala. "Clustering-based preconditioning for stochastic programs." *Computational Optimization and Applications* 64, no. 2 (2016): 379-406. https://doi.org/10.1007/s10589-015-9813-x

## BibTeX
```bibtex
@article{cao2016clustering,
  title   = {{Clustering-based preconditioning for stochastic programs}},
  author  = {Cao, Yankai and Laird, Carl D. and Zavala, Victor M.},
  journal = {Computational Optimization and Applications},
  volume  = {64},
  number  = {2},
  pages   = {379--406},
  year    = {2016},
  doi     = {10.1007/s10589-015-9813-x}
}
```

## Authors
Yankai Cao, Carl D. Laird, and Victor M. Zavala
