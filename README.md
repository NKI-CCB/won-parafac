# WON-PARAFAC

### Weighted orthogonal non-negative (WON) parallel factor analsyis (PARAFAC)

WON-PARAFAC is a variant of parallel factor analysis (PARAFAC), a tensor factorization method.
WON-PARAFAC impose the following three constraints on the standard PARAFAC:
1. Weighting scheme
- For balanced integration of the multiple data types
2. Orthogonality constraint
- To reduce overlapping between a factor (originally used on gene mode). This also introduces extra sparcity on the mode.
3. Non-negativity
- To induce sparse and parts-based representation.

### Implementation / Dependency

A multiplicative update rule was used to derive the algorithm, as in the original NMF implementation from [Lee & Seung (Nature, 1999)](https://www.nature.com/articles/44565).
The code requires [tensor toolbox version 2.6 (by Tamara Kolda)](https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
), freely available for non-commercial use upon registration.

For running the code, tenstor toobox must be avilable on the path environment, using `addpath` command in MATLAB.

### Demo code and data

You can load demo data, which contains pan-cancer multiomics data produced in [GDSC1000 project (Sanger)](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html).
You can load the data by:

```matlab
load Demo.mat
```
The command will load a varialbe `X`, a 3-way tensor (1815 gene by 935 cell lines and 5 data types).
Note that the 5 data types corresponds to below:
- positive gene expression levels (non-negative continuous; GE(+))
- absolute value of negative gene expression levels (non-negative continuous; GE(-))
- mutation (binary; MT)
- copy number gain (binary; CN(+))
- copy number loss (binary; CN(-))

The list of genes names in `X` is indicated in `genenames`, which will also be loaded together with `X`.

`Demo.m` will perform WON-PARAFAC analysis using random 100 genes by default, and varying number of factors and strength of orthogonal constraint on gene factor matrix.

- Number of basis: 10, 20, 30, ..., 200
- Strength of orthogonal constraint: 0 (no constraint), 0.2, 0.5, 1

Finally, a plot will be generated to show the performance of WON-PARAFAC for reconstructing input tensor (see below for an example).

![alt text](https://github.com/NKI-CCB/won-parafac/blob/master/Demo_plot.png "Demo plot")
