# PIPS_IPM_example

Examples of solving structured stochastic programs with PIPS-NLP and ClusterIPM through an early version of Plasmo.jl: stochastic farmer problems and a stochastic natural gas network.

## Contents
- `farmLink.jl`, `farmlarge.jl`: stochastic farmer problems
- `stochpdegas_PIPS.jl`, `set_data_paper.jl`: stochastic natural gas network model and its data setup
- `DAT_PAPER/`: gas network data and scenario generator
- `Plasmo.jl/`: early version of Plasmo.jl with interfaces to PIPS-NLP, ClusterIPM, and Ipopt
- `pipsnlp.parameter`: PIPS-NLP options
