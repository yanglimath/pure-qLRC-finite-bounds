# Finite-length comparison of nine qLRC bounds

`qlrc_finite_bounds.py` is the self-contained Python program used to generate the numerical data in Tables II and III. It compares nine upper bounds on the quantum dimension `kappa` for 80 parameter tuples. The program uses only the Python standard library, prints results to the console, and neither accesses the Internet nor creates files.

## Run

```bash
python3 qlrc_finite_bounds.py
```

## Tested parameters and sampling

The candidate ranges are

```text
q = 2, 3
delta = 3, 5, 7, 9, 11
30 <= n <= 130
2 <= r <= 6
```

For each fixed `(q,delta)`, the program forms all `(n,r)` pairs for which every bound has at least one admissible `kappa`, then samples eight pairs uniformly without replacement. A pair already selected for one `delta` is not reused for another `delta` with the same `q`.

One pseudorandom generator is initialized with seed `2026` before the loops over `q = 2,3` and the five increasing values of `delta`; its state is not reset. Thus the program produces

```text
2 * 5 * 8 = 80
```

distinct parameter tuples, exactly reproducing the two tables.

## Conversion to upper bounds on kappa

For each tuple `(q,delta,n,r)` and each bound, `max_kappa` enumerates every integer

```text
0 <= kappa <= n,    n + kappa = 0 (mod 2),
```

tests the bound and its applicability conditions, and returns the largest passing value. A smaller result is a tighter upper bound on `kappa`.

## CM-like bound and Grassl data

The CM-like bound requires `k_opt^(Q)(N,delta)`, the maximum possible dimension of a `Q`-ary linear code of length `N` and minimum distance at least `delta`. For every candidate `kappa`, the program sets

```text
Q = q^2
m = (n + kappa)/2
N = m - ell*(r + 1)
```

and tests

```text
kappa <= min_ell {ell*r + k_opt^(Q)(N,delta)},
```

where `ell = 0,...,floor((m-1)/(r+1))`. Hence `N` is determined separately for every tested `(kappa,ell)` and is not fixed beforehand.

For `Q = 4,9`, the Grassl table records the dimensions of the best-known linear codes. In this finite comparison, the program uses

```text
k_Grassl^(Q)(N,delta) = k_opt^(Q)(N,delta)
```

for all required values. This is the assumption underlying the displayed CM values; when the Grassl entry is not certified optimal, the resulting comparison is conditional.

The necessary Grassl data for `delta = 3,5,7,9,11` and `0 <= N <= 130` are included in compressed form. For each `(Q,delta)`, membership

```text
N in grassl_same[(Q,delta)]
```

means

```text
k_Grassl^(Q)(N,delta) = k_Grassl^(Q)(N-1,delta).
```

At every other length the dimension increases by one. Starting from `k_Grassl^(Q)(0,delta)=0`, `expand_grassl` restores the complete sequence. For example,

```python
grassl_same[(4, 3)] = (1, 2, 6, 22, 86)
k_Grassl[(4, 3)][7] = 4
```

because the value stays unchanged at `N=1,2,6` and increases at `N=3,4,5,7`.

## Main names

| Name | Meaning |
| --- | --- |
| `q_values`, `delta_values`, `n_values`, `r_values` | Candidate parameter ranges. |
| `sample_size`, `seed` | Eight samples per `(q,delta)` and random seed 2026. |
| `grassl_same` | Positions where the compressed Grassl dimension does not increase. |
| `k_Grassl` | Grassl dimensions restored for `N=0,...,130`. |
| `bounds` | Bound names, functions, and output order. |
| `expand_grassl` | Restores one compressed Grassl sequence. |
| `GG` | Tests the GG Singleton-like bound in (1). |
| `pure_S`, `pure_G`, `pure_P`, `pure_SP` | Test the four pure bounds in (7)--(10). |
| `ref25_G`, `ref25_CM`, `ref25_S`, `ref25_P` | Test the four bounds in [25, Theorem 6]. |
| `max_kappa` | Returns the largest parity-compatible `kappa` allowed by one bound. |
| `admissible` | Checks whether all nine bounds have a numerical result. |
| `parameter_tuples` | Reconstructs the 80 seeded random tuples. |
| `main` | Prints each tuple followed by its nine results. |

## Output order

1. GG Singleton-like bound in (1)
2. Pure Singleton-like bound in (7)
3. Pure Griesmer-like bound in (8)
4. Pure Plotkin-like bound in (9)
5. Pure Sphere-packing-like bound in (10)
6. [25, Theorem 6] Griesmer-like bound
7. [25, Theorem 6] CM-like bound
8. [25, Theorem 6] Singleton-like bound
9. [25, Theorem 6] Plotkin-like bound

The console output contains only `(q, delta, n, r)` and the nine final numerical results. The computation establishes comparisons only for the reproducible finite parameter set described above.
