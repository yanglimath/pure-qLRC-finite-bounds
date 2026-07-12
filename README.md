# Finite-length comparison of nine qLRC bounds

`qlrc_finite_bounds.py` is the self-contained Python program used to obtain the numerical data for Table II. It compares nine upper bounds on the quantum dimension `kappa` for 60 parameter tuples. The program uses only the Python standard library, prints the results to the console, and neither accesses the Internet nor creates output files.

## Run the program

With Python:

```bash
python3 qlrc_finite_bounds.py
```

## Parameter tuples

The program uses

```text
q = 2, 3,
delta = 3, 5, 7, 9, 11,
```

and the following six explicitly listed `(n,r)` pairs for every `(q,delta)`:

```text
delta = 3:  (38,2), (52,3), (67,4), (86,5), (100,6), (114,7)
delta = 5:  (45,2), (54,3), (69,4), (90,5), (102,6), (113,7)
delta = 7:  (42,2), (60,3), (68,4), (84,5),  (99,6), (120,7)
delta = 9:  (38,2), (53,3), (68,4), (86,5),  (99,6), (115,7)
delta = 11: (45,2), (60,3), (73,4), (85,5),  (97,6), (111,7)
```

The six positions cover the code-length strata `30-45`, `46-60`, `61-75`, `76-90`, `91-105`, and `106-120`, with locality `r=2,3,4,5,6,7`, respectively. Hence the output contains

```text
2 × 5 × 6 = 60
```

parameter blocks. Each block begins with `(q, delta, n, r)` and then lists the nine bounds in manuscript order.

## CM-like bound and Grassl table data

The CM-like bound requires the function `k_opt^(Q)(N,delta)`, which denotes the maximum possible dimension of a `Q`-ary linear code of length `N` and minimum distance at least `delta`. For a candidate `kappa`, the program sets

```text
Q = q^2,
m = (n+kappa)/2,
N = m - ell*(r+1),
```

and tests

```text
kappa <= min_ell {ell*r + k_opt^(Q)(N,delta)}.
```

The minimization is taken over every integer `ell` from `0` to `floor((m-1)/(r+1))`. Thus `N` is calculated separately for every tested pair `(kappa,ell)` and is not fixed in advance. For the 60 parameter tuples used in Table II, all required values satisfy `1 <= N <= 120`.

For quaternary and nonary linear codes, the [Grassl table](https://www.codetables.de) records the parameters of the best known linear codes. The program therefore uses the dimension in the lower-bound entry of the Grassl table as the optimal value required by the CM-like bound:

```text
k_Grassl^(Q)(N,delta) = k_opt^(Q)(N,delta),   Q = 4, 9.
```

The required Grassl table data for `delta = 3, 5, 7, 9, 11` and `0 <= N <= 120` are included directly in the Python program. The variable `grassl_same` stores these data in compressed form. For each `(Q,delta)`,

```text
N in grassl_same[(Q,delta)]
```

means that the recorded dimension does not increase from length `N-1` to length `N`:

```text
k_Grassl^(Q)(N,delta) = k_Grassl^(Q)(N-1,delta).
```

At every length not contained in `grassl_same[(Q,delta)]`, the dimension increases by one. Starting from `k_Grassl^(Q)(0,delta)=0`, `expand_grassl` applies this rule successively for `N=1,...,120` and reconstructs the complete dimension sequence. The reconstructed value is accessed in the CM-like bound as

```python
k_Grassl[(Q, delta)][N]
```

For example, consider `(q,delta,n,r,kappa)=(2,3,30,2,14)`. Then `Q=4`, `m=22`, and

```python
grassl_same[(4, 3)] = (1, 2, 6, 22, 86)
```

For `ell=5`, the residual length is `N=22-5*(2+1)=7`. Up to length 7, the dimension remains unchanged at `N=1,2,6` and increases at `N=3,4,5,7`. Hence `expand_grassl` gives

```python
k_Grassl[(4, 3)][7] = 4
```

and the corresponding CM term is

```text
ell*r + k_Grassl[(4,3)][N] = 5*2 + 4 = 14.
```

Testing every allowed `ell=0,...,7` gives CM terms `18, 18, 17, 16, 15, 14, 14, 14`. Their minimum is 14, so the candidate `kappa=14` satisfies the CM-like bound. The function `max_kappa` repeats this test for every admissible `kappa` and reports the largest passing value.

## Variables

| Name | Meaning |
| --- | --- |
| `q_values` | Quantum alphabet sizes used in Table II. |
| `delta_n_r` | The tested distances and the six `(n,r)` pairs assigned to each distance. |
| `grassl_same` | Compressed positions where the Grassl dimension remains unchanged as `N` increases. |
| `k_Grassl` | Complete reconstructed Grassl dimension sequences indexed by `(Q,delta)` and `N`. |
| `bounds` | Names, formula numbers, functions, and output order of the nine bounds. |
| `q` | Quantum alphabet size. |
| `Q` | Classical alphabet size `q^2` used by the Hermitian construction. |
| `delta` | Minimum distance. |
| `n` | Quantum code length. |
| `r` | Locality. |
| `kappa` | Candidate quantum dimension. |
| `m` | Auxiliary classical dimension `(n+kappa)/2`. |
| `ell` | Integer optimized over in the displayed bound. |
| `N` | Residual classical code length `m-ell*(r+1)`. |
| `V` | Hamming-ball volume used in the pure sphere-packing-like test. |

## Functions

| Function | Purpose |
| --- | --- |
| `expand_grassl` | Expands one compressed `grassl_same` tuple into values indexed by `N=0,...,120`. |
| `GG` | Tests the GG Singleton-like bound in (1). |
| `pure_S` | Tests the pure Singleton-like bound in (7). |
| `pure_G` | Tests the pure Griesmer-like bound in (8). |
| `pure_P` | Tests the pure Plotkin-like bound in (9). |
| `pure_SP` | Tests the pure sphere-packing-like bound in (10). |
| `ref25_G` | Tests the Griesmer-like bound in [25, Theorem 6]. |
| `ref25_CM` | Tests the CM-like bound in [25, Theorem 6] by enumerating every allowed `ell` and using `k=k_Grassl^(q^2)(N,delta)`. |
| `ref25_S` | Tests the Singleton-like bound in [25, Theorem 6]. |
| `ref25_P` | Tests the Plotkin-like bound in [25, Theorem 6]. |
| `max_kappa` | Enumerates all admissible `kappa` and returns the largest candidate allowed by one bound. |
| `main` | Prints the 60 parameter blocks and nine results per block. |

## Output order

1. GG Singleton-like bound in (1)
2. Pure Singleton-like bound in (7)
3. Pure Griesmer-like bound in (8)
4. Pure Plotkin-like bound in (9)
5. Pure sphere-packing-like bound in (10)
6. [25, Theorem 6] Griesmer-like bound
7. [25, Theorem 6] CM-like bound
8. [25, Theorem 6] Singleton-like bound
9. [25, Theorem 6] Plotkin-like bound

The console output contains only the parameter tuple and the nine final numerical results. The CM algorithm and the Grassl table data used for item 7 are documented in the preceding section and are not repeated in the output. A smaller reported value gives a tighter upper restriction on `kappa` for that parameter tuple. The computation establishes comparisons only over the explicitly listed finite parameter set.

## License

The source code is released under the MIT License.
