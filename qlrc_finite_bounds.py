# This Python program compares the tightness of nine finite-length qLRC bounds used in the manuscript tables.

from fractions import Fraction
from math import ceil, comb, floor
from random import Random


q_values = (2, 3)
delta_values = (3, 5, 7, 9, 11)
n_values = range(30, 131)
r_values = range(2, 7)
sample_size = 8
seed = 202607


grassl_same = {
    (4, 3): (1, 2, 6, 22, 86),
    (4, 5): (1, 2, 3, 4, 6, 12, 22, 44, 86),
    (4, 7): (1, 2, 3, 4, 5, 6, 8, 10, 18, 22, 27, 43, 47, 71, 114, 123),
    (4, 9): (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 15, 19, 28, 29, 32, 43, 52, 66, 69, 88, 112, 113),
    (4, 11): (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 17, 18, 30, 31, 33, 36, 40, 43, 65, 66, 70, 86, 89, 93, 94, 113),
    (9, 3): (1, 2, 11, 92),
    (9, 5): (1, 2, 3, 4, 11, 21, 73, 97),
    (9, 7): (1, 2, 3, 4, 5, 6, 11, 18, 23, 42, 53, 88),
    (9, 9): (1, 2, 3, 4, 5, 6, 7, 8, 11, 20, 21, 29, 42, 43, 61, 90),
    (9, 11): (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 17, 21, 29, 31, 41, 42, 82, 85, 92),
}


def expand_grassl(same):
    values = [0]
    for N in range(1, 131):
        values.append(values[-1] + (N not in same))
    return tuple(values)


k_opt = {key: expand_grassl(same) for key, same in grassl_same.items()}


def GG(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa:
        return False
    first_floor = floor(Fraction(n - (delta - 1), r + 1))
    second_floor = floor(Fraction(n - 2 * (delta - 1) - first_floor, r + 1))
    return kappa <= (n - 2 * (delta - 1) - first_floor - second_floor)


def pure_S(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa:
        return False
    return 2 * delta <= (n - kappa - 2 * ceil(Fraction(n + kappa, 2 * r)) + 4)


def pure_G(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa:
        return False
    ell_upper = ceil(Fraction(n + kappa, 2 * r)) - 1
    return n >= max(ell * (r + 1) + sum(ceil(Fraction(delta, q**t)) for t in range(0, n + kappa - 2 * ell * r - 1, 2)) for ell in range(0, ell_upper + 1))


def pure_P(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa:
        return False
    ell_upper = ceil(Fraction(n + kappa, 2 * r)) - 1
    return delta <= min(Fraction(q ** (n + kappa - 2 * ell * r - 2) * (q * q - 1) * (n - ell * (r + 1)), q ** (n + kappa - 2 * ell * r) - 1) for ell in range(0, ell_upper + 1))


def pure_SP(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa:
        return False
    ell_upper = floor(Fraction(n - 1, r + 1))
    radius = floor(Fraction(delta - 1, 2))
    for ell in range(0, ell_upper + 1):
        residual_length = n - ell * (r + 1)
        volume = sum(comb(residual_length, i) * (q * q - 1) ** i for i in range(0, radius + 1))
        exponent = n - kappa - 2 * ell
        if exponent < 0 or volume * volume > (q * q) ** exponent:
            return False
    return True


def ref26_G(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa or r >= kappa:
        return False
    ell_upper = ceil(Fraction(kappa, r)) - 1
    return n + kappa >= 2 * max(ell * (r + 1) + sum(ceil(Fraction(delta, q ** (2 * i))) for i in range(0, kappa - ell * r)) for ell in range(1, ell_upper + 1))


def ref26_CM(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa:
        return False
    ell_upper = floor(Fraction(n + kappa - 1, 2 * (r + 1)))
    return kappa <= min(ell * r + k_opt[(q * q, delta)][(n + kappa) // 2 - ell * (r + 1)] for ell in range(0, ell_upper + 1))


def ref26_S(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa:
        return False
    return 2 * delta <= (n - kappa - 2 * ceil(Fraction(kappa, r)) + 4)


def ref26_P(q, delta, n, r, kappa):
    if (n + kappa) % 2 != 0 or 2 * r > n + kappa or r >= kappa:
        return False
    ell_upper = ceil(Fraction(kappa, r)) - 1
    return 2 * delta <= min(Fraction(q ** (2 * (kappa - ell * r) - 2) * (q * q - 1) * (n + kappa - 2 * ell * (r + 1)), q ** (2 * (kappa - ell * r)) - 1) for ell in range(1, ell_upper + 1))


bounds = (
    ("GG Singleton-like bound in (1)", GG),
    ("Pure Singleton-like bound in (7)", pure_S),
    ("Pure Griesmer-like bound in (8)", pure_G),
    ("Pure Plotkin-like bound in (9)", pure_P),
    ("Pure Sphere-packing-like bound in (10)", pure_SP),
    ("[26, Theorem 6] Griesmer-like bound", ref26_G),
    ("[26, Theorem 6] CM-like bound", ref26_CM),
    ("[26, Theorem 6] Singleton-like bound", ref26_S),
    ("[26, Theorem 6] Plotkin-like bound", ref26_P),
)


def max_kappa(bound, q, delta, n, r):
    return max((kappa for kappa in range(0, n + 1) if bound(q, delta, n, r, kappa)), default=None)


def parameter_tuples():
    rng = Random(seed)
    for q in q_values:
        used = set()
        for delta in delta_values:
            pairs = [(n, r) for n in n_values for r in r_values if (n, r) not in used]
            selected = rng.sample(pairs, sample_size)
            used.update(selected)
            for n, r in selected:
                yield q, delta, n, r


def main():
    for q, delta, n, r in parameter_tuples():
        print(f"\n(q, delta, n, r) = ({q}, {delta}, {n}, {r})")
        for index, (name, bound) in enumerate(bounds, start=1):
            value = max_kappa(bound, q, delta, n, r)
            result = "not applicable" if value is None else value
            print(f"{index}. {name}: {result}")


if __name__ == "__main__":
    main()
