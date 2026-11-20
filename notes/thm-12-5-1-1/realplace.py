"""Compare Masullo's real-place integral with the classical archimedean term (F = Q).

  R(F)   := int_{R^x} f(|x|) |x|^{1/2} / |1 - x| dx/|x|     (f(1) = 0, so no principal value needed)
          = int F(u) e^{u/2} (1/|1 - e^u| + 1/(1 + e^u)) du,   F(u) = f(e^u)
  A(F)   := classical archimedean term for zeta (see classical.py), i.e.
            sum_rho h(gamma) = h(i/2) + h(-i/2) - primes + F(0) log(1/pi) + (1/2pi) int h(r) Re psi(1/4 + ir/2) dr
            and A(F) := F(0) log(1/pi) + (1/2pi) int h(r) Re psi(1/4 + ir/2) dr.
  The printed Thm 12.5.1.1 for F = Q (if mu_infty = dx/|x| and vol(hat Z^x) = 1) asserts
            (1/2) R(F) + primes = h(-i/2) + h(i/2) - sum_rho h(gamma),   i.e.  (1/2) R(F) = -A(F).
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classical import TestFunction, arch, bump

TESTS = [
    TestFunction("S1: (x-0.4)*bump((x-0.4)/2) * x", 0.4, 2.0, lambda x: x * (x - 0.4)),
    TestFunction("S2: bump(x-1.5)", 1.5, 1.0),
    TestFunction("S3: x^2*bump(x/2.5)", 0.0, 2.5, lambda x: x * x),
    TestFunction("S4: x*bump((x+0.7)/1.5)", -0.7, 1.5, lambda x: x),
    TestFunction("S5: x*exp(x)*bump((x-0.2)/1.8)", 0.2, 1.8, lambda x: x * np.exp(x)),
]


def R(T):
    # F(u) ~ F'(0) u near 0 makes F/|1-e^u| bounded but discontinuous (sign flip) at u = 0:
    # integrate on each side separately; Gauss-Legendre avoids evaluating at u = 0.
    tot = 0.0
    xg, wg = np.polynomial.legendre.leggauss(400)
    for lo, hi in ((max(T.lo, -50), min(0.0, T.hi)), (max(0.0, T.lo), T.hi)):
        if hi <= lo:
            continue
        # composite Gauss-Legendre on 200 panels
        edges = np.linspace(lo, hi, 201)
        uu = ((edges[:-1, None] + edges[1:, None]) / 2 + (edges[1:, None] - edges[:-1, None]) / 2 * xg[None, :]).ravel()
        ww = ((edges[1:, None] - edges[:-1, None]) / 2 * wg[None, :]).ravel()
        FF = T.F(uu)
        kk = np.exp(uu / 2) * (1 / np.abs(1 - np.exp(uu)) + 1 / (1 + np.exp(uu)))
        tot += float(np.sum(FF * kk * ww))
    return float(tot)


if __name__ == "__main__":
    for T in TESTS:
        r, a = R(T), arch(T, 0)
        a += float(T.F(0.0)) * math.log(1 / math.pi)
        print(f"{T.name:40s} F(0)={float(T.F(0.0)):+.1e}  R={r:+.10f}  -A={-a:+.10f}  R/(-A)={r/(-a):.8f}  R/2+A={r/2+a:+.3e}")
