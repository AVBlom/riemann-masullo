"""Classical explicit formula harness (control): validates zeros, quadrature and test functions.

Variables: x = log t.  A test function is F(x) = f(e^x), smooth with compact support.
  h(r) := int F(x) e^{i r x} dx          (so h(gamma) = int f(t) t^{i gamma} dt/t)
Weil explicit formula for a primitive Dirichlet L-function L(s, chi) of conductor q, parity a
(a = 0 even, a = 1 odd), real chi (zeros symmetric), rho = 1/2 + i gamma:
  sum_gamma h(gamma) = [chi trivial] (h(i/2) + h(-i/2))
                       - sum_n Lambda(n) n^{-1/2} (chi(n) F(log n) + conj chi(n) F(-log n))
                       + F(0) log(q/pi)
                       + (1/2pi) int h(r) Re psi((1/2 + a + i r)/2) dr
For zeta_{Q(i)} = zeta * L(., chi_-4) add the two identities.
"""
import math
import os

import numpy as np
from scipy.special import digamma

HERE = os.path.dirname(os.path.abspath(__file__))


def bump(y):
    out = np.zeros_like(y, dtype=float)
    m = np.abs(y) < 1
    out[m] = np.exp(-1.0 / (1.0 - y[m] ** 2))
    return out


class TestFunction:
    """F(x) = poly(x) * bump((x - c)/a); support [c-a, c+a]."""

    def __init__(self, name, c, a, poly=lambda x: np.ones_like(x)):
        self.name, self.c, self.a, self.poly = name, c, a, poly
        self.lo, self.hi = c - a, c + a
        n = 4001
        self.x = np.linspace(self.lo, self.hi, n)
        self.dx = self.x[1] - self.x[0]
        self.Fx = self.F(self.x)

    def F(self, x):
        x = np.asarray(x, dtype=float)
        return self.poly(x) * bump((x - self.c) / self.a)

    def h(self, r):
        """h(r) = int F(x) e^{i r x} dx for complex r (trapezoid; spectrally accurate for C_c^inf)."""
        r = np.atleast_1d(np.asarray(r, dtype=complex))
        out = np.empty(r.shape, dtype=complex)
        for i in range(0, len(r), 200):
            rr = r[i:i + 200, None]
            out[i:i + 200] = (self.Fx[None, :] * np.exp(1j * rr * self.x[None, :])).sum(axis=1) * self.dx
        return out


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return np.array([float(l) for l in fh if l.strip()])


def primes_upto(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


def prime_sum(T, chi):
    """sum_n Lambda(n) n^{-1/2} (chi(n) F(log n) + chi(n) F(-log n)) for real chi."""
    N = int(math.exp(max(abs(T.lo), abs(T.hi)))) + 2
    tot = 0.0
    for p in primes_upto(N):
        pk = p
        while pk <= N:
            lx = math.log(pk)
            tot += math.log(p) * pk ** -0.5 * chi(pk) * (T.F(lx) + T.F(-lx))
            pk *= p
    return float(tot)


def zero_sum(T, gammas):
    g = np.concatenate([gammas, -gammas])
    return T.h(g).real.sum(), abs(T.h(np.array([gammas[-1]]))[0])


def arch(T, a, R=1000.0, n=100001):
    r = np.linspace(-R, R, n)
    hr = T.h(r).real
    w = digamma((0.5 + a + 1j * r) / 2).real
    return float(np.trapezoid(hr * w, r) / (2 * math.pi))


def explicit(T, gammas, q, a, chi, trivial):
    zs, tail = zero_sum(T, gammas)
    pole = (T.h(np.array([0.5j]))[0] + T.h(np.array([-0.5j]))[0]).real if trivial else 0.0
    ps = prime_sum(T, chi)
    cond = float(T.F(0.0)) * math.log(q / math.pi)
    ar = arch(T, a)
    rhs = pole - ps + cond + ar
    return dict(zeros=zs, pole=pole, primes=ps, cond=cond, arch=ar, rhs=rhs, diff=zs - rhs, last_h=tail)


TESTS = [
    TestFunction("T1: x*bump(x/2)", 0.0, 2.0, lambda x: x),
    TestFunction("T2: bump(x-1.5)", 1.5, 1.0),
    TestFunction("T3: x^2*bump(x/2.5)", 0.0, 2.5, lambda x: x * x),
    TestFunction("T4: bump((x+0.3)/1.2)  [F(0)!=0]", -0.3, 1.2),
]


def chi4(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


if __name__ == "__main__":
    zz, cz = load("zeta_zeros.txt"), load("chi4_zeros.txt")
    for T in TESTS:
        e1 = explicit(T, zz, 1, 0, lambda n: 1, True)
        e2 = explicit(T, cz, 4, 1, chi4, False)
        print(T.name)
        print(f"   zeta : zeros={e1['zeros']:+.12f} rhs={e1['rhs']:+.12f} diff={e1['diff']:+.2e} |h(gamma_max)|={e1['last_h']:.1e}")
        print(f"   chi-4: zeros={e2['zeros']:+.12f} rhs={e2['rhs']:+.12f} diff={e2['diff']:+.2e} |h(gamma_max)|={e2['last_h']:.1e}")
