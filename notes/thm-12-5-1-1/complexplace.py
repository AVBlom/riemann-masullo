"""Complex place (F = Q(i)): compare candidate readings of Masullo's term with the classical one.

phi(g) = f(|g|_K), |z|_v = |z|^2 at the complex place, F(u) = f(e^u), F(0) = 0.
Classical archimedean term for zeta_{Q(i)} = zeta * L(., chi_-4) (conductor terms vanish as F(0)=0):
    A_C(F) = (1/2pi) int h(r) [Re psi(1/4 + ir/2) + Re psi(3/4 + ir/2)] dr
Masullo's printed term: (1/2pi) int_{C^x} f(|z|^2) |z|_v^{1/2} / |1 - z|^e dmu_C(z),
with e = 2 (normalized absolute value |1-z|_C) or e = 1 (ordinary modulus), and
dmu_C = dx dy / |z|^2  (reported ratios scale linearly for other constants).
   e = 2:  int_{C^x} ... = pi * int F(u) e^{u/2} / |1 - e^u| du
   e = 1:  int_{C^x} ... = int_0^inf f(r^2) * [4 K(k)/(1+r)] dr,   k^2 = 4r/(1+r)^2
The identity needs (printed term) = -A_C(F).
"""
import math
import os
import sys

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classical import TestFunction, arch  # noqa: E402
from realplace import TESTS  # noqa: E402

xg, wg = np.polynomial.legendre.leggauss(200)


def panels(lo, hi, n=200):
    e = np.linspace(lo, hi, n + 1)
    u = ((e[:-1, None] + e[1:, None]) / 2 + (e[1:, None] - e[:-1, None]) / 2 * xg[None, :]).ravel()
    w = ((e[1:, None] - e[:-1, None]) / 2 * wg[None, :]).ravel()
    return u, w


def split_int(T, kernel):
    tot = 0.0
    for lo, hi in ((T.lo, min(0.0, T.hi)), (max(0.0, T.lo), T.hi)):
        if hi > lo:
            u, w = panels(lo, hi)
            tot += float(np.sum(T.F(u) * kernel(u) * w))
    return tot


def e2(T):
    return math.pi * split_int(T, lambda u: np.exp(u / 2) / np.abs(1 - np.exp(u)))


def e1(T):
    # variable u = log(r^2) = 2 log r, r = e^{u/2}, dr = r du / 2
    def kern(u):
        r = np.exp(u / 2)
        k2 = 4 * r / (1 + r) ** 2
        K = np.array([float(mp.ellipk(x)) for x in k2])  # mpmath ellipk takes m = k^2
        return 4 * K / (1 + r) * r / 2
    return split_int(T, kern)


if __name__ == "__main__":
    for T in TESTS:
        A = arch(T, 0) + arch(T, 1)
        a2, a1 = e2(T) / (2 * math.pi), e1(T) / (2 * math.pi)
        print(f"{T.name:40s} -A_C={-A:+.10f}  printed(e=2)={a2:+.10f} ratio={a2/(-A):.8f}   "
              f"printed(e=1)={a1:+.10f} ratio={a1/(-A):.6f}")
