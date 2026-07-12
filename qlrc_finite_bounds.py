# This Python program compares the tightness of the nine bounds mentioned in the manuscript for finite code lengths; the resulting data are used to prepare Table II.

from fractions import Fraction
from math import comb


q_values = (2, 3)
delta_n_r = {
    3: ((30, 2), (50, 3), (70, 5), (90, 7), (110, 9), (130, 12)),
    5: ((30, 2), (50, 3), (70, 5), (90, 7), (110, 9), (130, 12)),
    8: ((30, 2), (50, 3), (70, 5), (90, 7), (110, 9), (130, 12)),
    10: ((40, 2), (60, 3), (80, 5), (100, 7), (120, 9), (130, 12)),
    12: ((50, 2), (70, 3), (90, 5), (110, 7), (120, 9), (130, 12)),
}


grassl_same = {
    (4, 3): (1, 2, 6, 22, 86),
    (4, 5): (1, 2, 3, 4, 6, 12, 22, 44, 86),
    (4, 8): (
        1, 2, 3, 4, 5, 6, 7, 9, 11, 19, 21, 27, 37, 40, 66, 87, 91, 93,
    ),
    (4, 10): (
        1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 16, 19, 29, 30, 32,
        37, 41, 54, 67, 69, 87, 93, 97,
    ),
    (4, 12): (
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 18, 19, 31,
        32, 33, 37, 39, 43, 66, 67, 69, 72, 87, 91, 94, 114, 115, 121, 130,
    ),
    (9, 3): (1, 2, 11, 92),
    (9, 5): (1, 2, 3, 4, 11, 21, 73, 97),
    (9, 8): (
        1, 2, 3, 4, 5, 6, 7, 11, 19, 21, 29, 42, 83, 89,
    ),
    (9, 10): (
        1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 21, 22, 29, 31, 40, 62, 81, 91,
    ),
    (9, 12): (
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 17, 21, 29, 31,
        42, 43, 83, 84, 87, 93,
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


def main():
    for q in q_values:
        for delta, pairs in delta_n_r.items():
            for n, r in pairs:
                print(f"\n(q, delta, n, r) = ({q}, {delta}, {n}, {r})")
                for index, (name, bound) in enumerate(bounds, start=1):
                    print(f"{index}. {name}: {max_kappa(bound, q, delta, n, r)}")


if __name__ == "__main__":
    main()
