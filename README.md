### Dependencies 

Sagemath (any sufficiently recent version should be fine)

### Description of files:

- $\mathtt{test}\textunderscore\mathtt{algorithm.ipynb}$: executes the sampling algorithm for a single code; does also some sanity-checks to verify rank, row weight, self-orthogonality, etc.
- $\mathtt{verify}\textunderscore\mathtt{distribution.ipynb}$: samples random codes from $\mathcal H_{n,r,v}$, computes experimentally the Expected Weight Distribution and compares it with the theoretical one
- $\mathtt{test}\textunderscore\mathtt{ISD.sage}$: executes Lee and Brickell ISD on codes sampled from $\mathcal H_{n,r,v}$ and measures relevant quantities (success of calls, number of found codewords, etc.)
- $\mathtt{benchmark}\textunderscore\mathtt{algorithm.sage}$: executes the sampling algorithm for several codes and measures success rate, timings, etc.

