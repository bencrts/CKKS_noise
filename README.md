# CKKS Precision Loss
A repository for code used in the paper "On the precision loss in approximate homomorphic encryption"

## File structure

```
+-- HEAAN
|   +-- An adapted version of the C++ HEAAN v1.0 library. Experimental code is in /run.
+-- noise-heuristics
|   +-- Python files for noise heuristics
```

Changes to the v1 HEAAN library can be viewed [here](https://github.com/bencrts/CKKS_noise/commit/1630c5e6694d7772c4eedd1b1aa9c403e6ae1fb1).

## Using the code

To re-run the experiments used in the paper, the following steps should be taken. <br>

1. For the noise heuristics, navigate to `/noise-heuristics` and run the python file main.py via e.g. `python3 main.py`. <br>
2. For the (ring-to-ring) experimental results, navigate to `/HEAAN_v1.0/HEAAN/run` and use the command `make encryption`. <br>
3. For the (complex-to-complex) experimental results, navigate to `/HEAAN_v1.0/HEAAN/run` and use the command `make complex`. <br>

In each case, results will print to the terminal. Note that we do not fix seeds in this public implementation, so running the experiments multiple times will yield different results. 

## Example

As an example, we consider running the Python heuristics.

```
user@PC % python3 main.py
------- CE Results -------
{'type': 'Rq -> Rq CE 6√V', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 52.35, 'fresh-noise': 17.5, '+': 18.5, 'x': 31.8}
{'type': 'Rq -> Rq CE', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 52.35, '+': 18.04, 'x': 31.39}
{'type': 'RR->RR CE', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 40.0, '+': -21.93, 'x': -20.91}
{'type': 'CC->CC CE', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 40.5, '+': -21.93, 'x': -20.42}
{'type': 'CC->CC CE 6√V', 'log(N)': 13.0, 'log(Δ)': 40.0, 'fresh-prec-loss': -22.51, 'α': 1e-05, 'log(MB)': 40.5, '+': -21.51, 'x': -19.99}



------- WCR Results -------
{'type': 'Rq->Rq WCR 6rootV', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': None, 'fresh noise': 10.97, 'log(MB)': 40.0, '+': 11.97, 'x': 25.97}
{'type': 'Rq->Rq WCR', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 40.0, '+': 11.99, 'x': 25.99}
{'type': 'RR->RR WCR', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 40.0, '+': -15.66, 'x': -1.66}
{'type': 'CC->CC WCR', 'log(N)': 13.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 40.5, '+': -15.66, 'x': -1.16}



------- CLT Results -------
{'type': 'No-α CLT', 'log(N)': 13.0, 'log(q)': 109.0, 'log(Δ)': 40.0, 'α': None, 'fresh noise': 10.97, 'log(MB)': 93.0, '+': 11.47, 'x': 18.76}
{'type': 'Rq->Rq CLT', 'log(N)': 13.0, 'log(q)': 109.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 93.0, '+': 11.49, 'x': 18.78}
{'type': 'RR->RR CLT', 'log(N)': 13.0, 'log(q)': 109.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 80.0, 'c+': -22.46, 'cx': -21.67, 'r+': -22.54, 'rx': -21.74}
{'type': 'CC->CC CLT', 'log(N)': 13.0, 'log(q)': 109.0, 'log(Δ)': 40.0, 'α': 1e-05, 'log(MB)': 81.0, 'c+': -22.46, 'cx': -21.17, 'r+': -22.54, 'rx': -21.25}
```
