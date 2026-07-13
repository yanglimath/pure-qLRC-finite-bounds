# Finite-length comparison of nine qLRC bounds

`qlrc_finite_bounds.py` is the self-contained Python program used for the
finite-length comparison tables in the manuscript. It evaluates nine upper
bounds on the quantum dimension `kappa` using only the Python standard
library. The program prints results to the console and creates no files.

## Run

```bash
python3 qlrc_finite_bounds.py
```

## Parameters and sampling

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

```text
2 * 5 * 8 = 80
```

distinct and exactly reproducible parameter tuples.

## Enumeration and applicability

For every fixed `(q,delta,n,r)`, `max_kappa` checks every integer

```text
0 <= kappa <= n
```

and returns the largest candidate accepted by the selected bound. Each of the
nine predicates directly rejects odd `n+kappa` and candidates with
`2r > n+kappa`, because the associated Hermitian classical dimension
`k=(n+kappa)/2` must be integral and the underlying classical LRC satisfies
`r <= k`. In addition, `ref26_G` and `ref26_P` require `r < kappa`; the other
seven bounds have no further applicability condition.

If a bound accepts no candidate, `max_kappa` returns `None`, which the console
displays as `not applicable`. Such a result does not remove the sampled
parameter tuple. A smaller numerical value is a tighter upper restriction on
`kappa`.

## Exact arithmetic and CM data

- Ceilings and floors use `math.ceil` and `math.floor` applied to exact
  `fractions.Fraction` values.
- Both Plotkin-like bounds use exact rational arithmetic.
- The sphere-packing-like bound is checked with arbitrary-precision integer
  powers, without logarithms or floating-point approximations.
- The CM-like bound enumerates every allowed `ell` and evaluates

  ```text
  Q = q^2
  m = (n+kappa)/2
  N = m - ell*(r+1)
  kappa <= min_ell {ell*r + k_opt^(Q)(N,delta)}.
  ```

The embedded `grassl_same` data reconstruct the required Grassl-table
dimensions for `Q=4,9`, `delta=3,5,7,9,11`, and `0 <= N <= 130`. The CM
comparison assumes that each recorded best-known dimension equals the exact
value `k_opt^(Q)(N,delta)`; wherever optimality is not certified, the CM
comparison is conditional on this assumption.

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
