# This Python program compares the tightness of the nine bounds mentioned in the manuscript for finite code lengths; the resulting data are used to prepare Table II.

from fractions import Fraction
from math import comb
from random import Random


q_values = (2, 3)
delta_values = (3, 5, 7, 9, 11)
n_values = range(30, 131)
r_values = range(2, 7)
sample_size = 8
seed = 2026


grassl_same = {
    (4, 3): (1, 2, 6, 22, 86),
    (4, 5): (1, 2, 3, 4, 6, 12, 22, 44, 86),
    (4, 7): (
        1, 2, 3, 4, 5, 6, 8, 10, 18, 22, 27, 43, 47, 71, 114, 123,
    ),
    (4, 9): (
        1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 15, 19, 28, 29, 32, 43,
        52, 66, 69, 88, 112, 113,
    ),
    (4, 11): (
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 17, 18, 30, 31,
        33, 36, 40, 43, 65, 66, 70, 86, 89, 93, 94, 113,
    ),
    (9, 3): (1, 2, 11, 92),
    (9, 5): (1, 2, 3, 4, 11, 21, 73, 97),
    (9, 7): (
        1, 2, 3, 4, 5, 6, 11, 18, 23, 42, 53, 88,
    ),
    (9, 9): (
        1, 2, 3, 4, 5, 6, 7, 8, 11, 20, 21, 29, 42, 43, 61, 90,
    ),
    (9, 11): (
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 17, 21, 29, 31, 41,
        42, 82, 85, 92,
    ),
}


def expand_grassl(same):
    values = [0]
    for N in range(1, 131):
        values.append(values[-1] + (N not in same))
    return tuple(values)


k_Grassl = {key: expand_grassl(same) for key, same in grassl_same.items()}


def GG(q, delta, n, r, kappa):
    return kappa <= n - 2 * (delta - 1) - (n - delta + 1) // (r + 1) - (n - 2 * (delta - 1) - (n - delta + 1) // (r + 1)) // (r + 1)


def pure_S(q, delta, n, r, kappa):
    if 2 * r > n + kappa:
        return False
    return 2 * delta <= n - kappa - 2 * ((n + kappa + 2 * r - 1) // (2 * r)) + 4


def pure_G(q, delta, n, r, kappa):
    if 2 * r > n + kappa:
        return False
    return n >= max(ell * (r + 1) + sum((delta + q**t - 1) // q**t for t in range(0, n + kappa - 2 * ell * r - 1, 2)) for ell in range((n + kappa + 2 * r - 1) // (2 * r)))


def pure_P(q, delta, n, r, kappa):
    if 2 * r > n + kappa:
        return False
    return delta <= min(Fraction(q ** (n + kappa - 2 * ell * r - 2) * (q * q - 1) * (n - ell * (r + 1)), q ** (n + kappa - 2 * ell * r) - 1) for ell in range((n + kappa + 2 * r - 1) // (2 * r)))


def pure_SP(q, delta, n, r, kappa):
    if 2 * r > n + kappa:
        return False
    for ell in range((n - 1) // (r + 1) + 1):
        V = sum(comb(n - ell * (r + 1), i) * (q * q - 1) ** i for i in range(min((delta - 1) // 2, n - ell * (r + 1)) + 1))
        if n - kappa - 2 * ell < 0:
            return False
        if V * V > (q * q) ** (n - kappa - 2 * ell):
            return False
    return True


def ref25_G(q, delta, n, r, kappa):
    if kappa <= r:
        return True
    return n + kappa >= 2 * max(ell * (r + 1) + sum((delta + q ** (2 * i) - 1) // q ** (2 * i) for i in range(kappa - ell * r)) for ell in range(1, (kappa + r - 1) // r))


def ref25_CM(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0:
        return False
    m = (n + kappa) // 2
    return kappa <= min(ell * r + k_Grassl[(q * q, delta)][m - ell * (r + 1)] for ell in range((m - 1) // (r + 1) + 1))


def ref25_S(q, delta, n, r, kappa):
    return 2 * delta <= n - kappa - 2 * ((kappa + r - 1) // r) + 4


def ref25_P(q, delta, n, r, kappa):
    if kappa <= r:
        return True
    return 2 * delta <= min(Fraction(q ** (2 * (kappa - ell * r) - 2) * (q * q - 1) * (n + kappa - 2 * ell * (r + 1)), q ** (2 * (kappa - ell * r)) - 1) for ell in range(1, (kappa + r - 1) // r))


bounds = (
    ("GG Singleton-like bound in (1)", GG),
    ("Pure Singleton-like bound in (7)", pure_S),
    ("Pure Griesmer-like bound in (8)", pure_G),
    ("Pure Plotkin-like bound in (9)", pure_P),
    ("Pure Sphere-packing-like bound in (10)", pure_SP),
    ("[25, Theorem 6] Griesmer-like bound", ref25_G),
    ("[25, Theorem 6] CM-like bound", ref25_CM),
    ("[25, Theorem 6] Singleton-like bound", ref25_S),
    ("[25, Theorem 6] Plotkin-like bound", ref25_P),
)


def max_kappa(bound, q, delta, n, r):
    return max((kappa for kappa in range(n + 1) if (n + kappa) % 2 == 0 and bound(q, delta, n, r, kappa)), default=None)


def admissible(q, delta, n, r):
    return all(any(bound(q, delta, n, r, kappa) for kappa in range(n % 2, n + 1, 2)) for _, bound in bounds)


def parameter_tuples():
    rng = Random(seed)
    for q in q_values:
        used = set()
        for delta in delta_values:
            pairs = [(n, r) for n in n_values for r in r_values if (n, r) not in used and admissible(q, delta, n, r)]
            selected = rng.sample(pairs, sample_size)
            used.update(selected)
            for n, r in selected:
                yield q, delta, n, r


def main():
    for q, delta, n, r in parameter_tuples():
        print(f"\n(q, delta, n, r) = ({q}, {delta}, {n}, {r})")
        for index, (name, bound) in enumerate(bounds, start=1):
            print(f"{index}. {name}: {max_kappa(bound, q, delta, n, r)}")


if __name__ == "__main__":
    main()
