# Finite-length comparison of nine qLRC bounds

`qlrc_finite_bounds.py` is the self-contained Python program used for the
finite-length comparison tables in the manuscript [**LLLL2026**]. It evaluates nine upper
bounds on the quantum dimension `kappa` using only the Python standard
library. The program prints results to the console and creates no files.

[**LLLL2026**] Y. Li, S., G. Luo, and S. Ling, *Improved bounds and optimal constructions of pure
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

- For every fixed `(q,delta,n,r)`, `max_kappa` checks every integer
`0 <= kappa <= n` and returns the largest candidate accepted by the selected bound. 
All these bounds hold under the same pure assumption and their own conditions. 


- If a bound accepts no candidate, `max_kappa` returns `None`, which the console
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
The variable `grassl_same` is used to reconstruct the largest known dimensions of linear codes recorded in Grassl’s code table [14] for `Q in {4, 9}`, `delta in {3, 5, 7, 9, 11}`, and `0 <= N <= 130`. For the numerical evaluation, we treat these best-known dimensions as the exact values of `k_opt^(Q)(N, delta)`. Therefore, whenever the optimality of a recorded dimension has not been established, the resulting CM-like bound and the associated numerical comparisons are conditional on this assumption. 


### A simple example

Consider `Q = 4` and `delta = 3`. The program stores
`grassl_same[(4, 3)] = (1, 2, 6, 22, 86)`. 
The reconstruction starts with 
`values[0] = 0`. 
The function then checks each length `N` in order.

When `N = 1`, the number `1` appears in `grassl_same[(4, 3)]`. Therefore, the dimension remains `0`.

When `N = 2`, the number `2` also appears in the tuple. Therefore, the dimension remains `0`.

When `N = 3`, the number `3` does not appear in the tuple. Therefore, the dimension increases from `0` to `1`.

The same rule gives

```text
N:          0  1  2  3  4  5  6  7
Best-known dimension:  0  0  0  1  2  3  3  4
```

For example, the dimension remains `3` when `N = 6` because `6` appears in the tuple. It then increases to `4` when `N = 7` because `7` does not appear in the tuple. This rule is implemented by `values.append(values[-1] + (N not in same))`
In Python, `N not in same` has value `0` when `N` appears in the tuple and value `1` when it does not. In this way, the function reconstructs the complete list of largest known dimensions for `0 <= N <= 130`.
