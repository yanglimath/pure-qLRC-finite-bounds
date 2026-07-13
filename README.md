# Finite-length comparison of nine qLRC bounds

`qlrc_finite_bounds.py` is the self-contained Python program used for the
finite-length comparison tables in the manuscript [**LLLL2026**]. It evaluates nine upper
bounds on the quantum dimension `kappa` using only the Python standard
library. The program prints results to the console and creates no files.

[**LELLS2026**] Y. Li, S., G. Luo, and S. Ling, *Improved bounds and optimal constructions of pure
quantum locally recoverable codes*, preprint, 2026. 


## Parameters and random sampling

The candidate ranges are

```text
q = 2, 3
delta = 3, 5, 7, 9, 11
30 <= n <= 130
2 <= r <= 6
```

For each fixed `(q,delta)`, the program uses
`random.Random(202607).sample` to select eight `(n,r)` pairs uniformly without
replacement from the full candidate set. For a fixed `q`, a pair selected
for one `delta` is excluded from later distances. Bound applicability does
not affect parameter selection. Hence the program produces
2 * 5 * 8 = 80 distinct and exactly reproducible parameter tuples.

## Some notes

-For every fixed `(q,delta,n,r)`, `max_kappa` checks every integer
`0 <= kappa <= n` and returns the largest candidate accepted by the selected bound. 
All these bounds hold under the same pure assumption and their own condititons. 


-If a bound accepts no candidate, `max_kappa` returns `None`, which the console
displays as `not applicable`. Such a result does not remove the sampled
parameter tuple. A smaller numerical value is a tighter upper restriction on
`kappa`.

- The CM-like bound in [26, Theorem 6] enumerates every allowed `ell` and evaluates

  ```text
  Q = q^2
  m = (n+kappa)/2
  N = m - ell*(r+1)
  kappa <= min_ell {ell*r + k_opt^(Q)(N,delta)}.
  ```
The variable `grassl_same` is used to reconstruct the largest known dimensions of linear codes recorded in Grassl’s code table [14] for 
\(Q\in\{4,9\}\), \(\delta\in\{3,5,7,9,11\}\), and \(0\leq N\leq130\). For the numerical evaluation, we treat these best-known dimensions as the exact values of \(k_{\mathrm{opt}}^{(Q)}(N,\delta)\). Therefore, whenever the optimality of a recorded dimension has not been established, the resulting CM-like bound and the associated numerical comparisons are conditional on this assumption.

## Output order

1. GG Singleton-like bound in (1)
2. Pure Singleton-like bound in (7)
3. Pure Griesmer-like bound in (8)
4. Pure Plotkin-like bound in (9)
5. Pure Sphere-packing-like bound in (10)
6. [26, Theorem 6] Griesmer-like bound
7. [26, Theorem 6] CM-like bound
8. [26, Theorem 6] Singleton-like bound
9. [26, Theorem 6] Plotkin-like bound
