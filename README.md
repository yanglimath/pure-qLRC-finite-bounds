# Finite-length comparison of nine qLRC bounds

`qlrc_finite_bounds_sagemath.py` is the self-contained program used to obtain the numerical data for Table II. It compares nine upper bounds on the quantum dimension `kappa` for 60 parameter tuples. The program uses only the Python standard library, prints the results to the console, and neither accesses the Internet nor creates output files.

## Run the program

With Python:

```bash
python3 qlrc_finite_bounds_sagemath.py
```

## Parameter tuples

The program uses

```text
q = 2, 3,
delta = 3, 5, 8, 10, 12,
```

and six explicitly listed `(n,r)` pairs for every `(q,delta)`. Hence the output contains

```text
2 × 5 × 6 = 60
```

parameter blocks. Each block begins with `(q, delta, n, r)` and then lists the nine bounds in manuscript order.

## How the residual length N is determined

The CM-like bound does not require the optimal value of `kappa` to be known in advance. For fixed `(q,delta,n,r)`, `max_kappa` enumerates every integer

```text
0 <= kappa <= n
```

that satisfies `n+kappa` even. For each candidate `kappa`, `ref25_CM` sets

```text
m = (n+kappa)/2
ell_max = floor((m-1)/(r+1))
```

and enumerates every integer `ell` from `0` through `ell_max`. The corresponding residual length is calculated at that moment as

```text
N = m - ell*(r+1).
```

Thus `N` is generated separately for every tested pair `(kappa,ell)`. It is not selected before the search and it is not tied to a preselected `ell`.

Because `0 <= kappa <= n`, one has `m <= n`; because of the definition of `ell_max`, one has `N >= 1`. All Table II cases satisfy `n <= 130`, so every required residual length satisfies

```text
1 <= N <= 130.
```

The embedded Grassl snapshot therefore stores values for all indices `0 <= N <= 130`. Index `N=0` is included only to simplify reconstruction of the stored sequences.

For each candidate `kappa`, the program evaluates

```text
min_ell {ell*r + k_Grassl^(Q)(N,delta)},   Q = q^2,
```

over the complete allowed range of `ell`. The candidate passes the CM-like test when this minimum is at least `kappa`. The largest passing candidate is the reported upper value of `kappa`.

## Grassl snapshot and compression

The fixed snapshot was obtained from [codetables.de](https://www.codetables.de) on 2026-07-13 for

```text
Q = 4, 9,
delta = 3, 5, 8, 10, 12,
0 <= N <= 130.
```

`k_Grassl[(Q,delta)][N]` is the largest dimension recorded by the snapshot for a known linear code of length `N` and minimum distance at least `delta`. The finite-length comparison adopts the explicit numerical assumption

```text
k_Grassl^(Q)(N,delta) = k_opt^(Q)(N,delta).
```

For a Grassl entry whose lower and upper bounds do not coincide, this is a computational convention rather than an independently proved optimality statement. The assumption is documented here instead of being repeated in every console-output block.

The complete dimension sequences are stored losslessly through `grassl_same`. For each `(Q,delta)`,

```text
N in grassl_same[(Q,delta)]
```

means

```text
k_Grassl^(Q)(N,delta) = k_Grassl^(Q)(N-1,delta).
```

At every length not listed in `grassl_same`, the stored dimension increases by one. `expand_grassl` starts from the value zero at `N=0` and reconstructs the complete tuple of 131 values.

For example, for `(Q,delta)=(4,8)`, the number `19` occurs in `grassl_same[(4,8)]`, whereas `20` does not. Consequently,

```text
k_Grassl^(4)(18,8) = 9,
k_Grassl^(4)(19,8) = 9,
k_Grassl^(4)(20,8) = 10.
```

The snapshot is embedded data, not a crawler. `expand_grassl` does not read a file and does not access the Internet.

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
| `expand_grassl` | Expands one compressed `grassl_same` tuple into values indexed by `N=0,...,130`. |
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

The console output contains only the parameter tuple and the nine final numerical results. The CM algorithm, Grassl snapshot, and assumption used for item 7 are documented in the preceding sections and are not repeated in the output. A smaller reported value gives a tighter upper restriction on `kappa` for that parameter tuple. The computation establishes comparisons only over the explicitly listed finite parameter set.

## License

The source code is released under the MIT License.
