# Installation

The first step is to install the `HEAAN` and `FullRNSHEAAN` libraries, which can be found in the parent directory.

To install HEAAN, cd to `CKKS_noise/HEAAN/HEAAN/lib` and run `make all`.
To install FullRNSHEAAN, cd to `CKKS_noise/FullRNS-HEAAN/lib` and run `make all`.

# Experiments

To run the textbook experiments, do:

`make textbook`

the output should look like:

```
The parameters are: log(N) = 13, log(Q) = 109, log(Delta) = 40, h = 64.
over 10 trials, the results are:
REAL ADDITION, average = 2.56985e-08, maximum = 6.03574e-08
REAL MULTIPLICATION, average = 3.60548e-08, maximum = 6.70139e-08
COMPLEX ADDITION, average = 3.56668e-08, maximum = 6.06288e-08
COMPLEX MULTIPLICATION, average = 4.53978e-08, maximum = 8.98513e-08
RING ADDITION, average = 24.3, maximum = 31
RING MULTIPLICATION, average = 36.1, maximum = 46

[...]
```

To run the RNS experiments, do:

`make rns`

the output should look like:

```
The parameters are: log(N) = 12, log(Q) = 100, log(Delta) = 40.
over 10 trials, the results are:
REAL ADDITION, average = 4.35007e-08, maximum = 4.90518e-08
REAL MULTIPLICATION, average = 2.36729e-07, maximum = 2.4406e-07
COMPLEX ADDITION, average = 5.5694e-08, maximum = 6.05713e-08
COMPLEX MULTIPLICATION, average = 4.35002e-07, maximum = 4.52752e-07

[...]
```

The `loop` parameter, which determines the number of iterations, is set to 10 in both `textbook-ckks.cpp` and `rns-ckks.cpp` for speed purposes.
To increase the number of iterations, manually increase this parameter in both files and re-run the code.
