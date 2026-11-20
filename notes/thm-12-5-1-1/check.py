"""Numerical check of Masullo, *The Music of the Primes*, Theorem 12.5.1.1 (pp. 519-520),
for F = Q, Q(i) and Q(sqrt-3), with test functions phi(g) = f(|g|_F) depending only on the idele norm.

Such phi lie in D(C_F) (they are invariant under all local units, smooth, compactly supported because
C_F^1 is compact) and satisfy phi(C_F^1) = 0 iff f(1) = 0.  Write F(x) := f(e^x), h(r) := int F(x) e^{irx} dx.

Printed LHS, with the paper's normalizations (Def. A.3.2.2, p. 560: mu_v = dx/|x| at real v,
2 dx dy/|z|^2 at complex v; |.|_v normalized, so |z|_v = |z|^2; Def. B.4.4.1 / Lemma B.4.4.3: with
phi(1) = 0 the principal values are ordinary integrals):
  real v     : (1/2)    int_{R^x} f(|x|) |x|^{1/2} / |1 - x|        dx/|x|
  complex v  : (1/2pi)  int_{C^x} f(|z|^2) |z| / |1 - z|^2        2 dx dy/|z|^2
  finite v   : log(q_v)/mu_v(O_v^x) int f(|g|_v)|g|_v^{1/2}/|1-g|_v dmu_v
               = log q_v * sum_{k>=1} q_v^{-k/2} (F(k log q_v) + F(-k log q_v))     (exact; f(1) = 0)
Printed RHS: M(phi)(|.|^{1/2}) + M(phi)(|.|^{-1/2}) - sum_{zeros c in interior of 1/2 + I} N(c) M(phi)(c|.|^{-1/2})
  with M(phi)(c) = int_{C_F} phi c d(mu-bar), mu-bar Tate's measure (Def. A.3.2.1), so for norm-only phi
  M(phi)(|.|^s) = vol * int f(t) t^s dt/t,  vol = mu-bar_1(C_F^1) = 2^{r1}(2pi)^{r2} h R/(sqrt|d| w)
  (= 1 for Q, pi/4 for Q(i), pi/(3 sqrt 3) for Q(sqrt-3)).  Only zeros of zeta_F contribute (other Hecke
  characters integrate to 0 over C_F^1).  A zero c = |.|^{1/2 + i gamma} on the critical line contributes
  vol * h(gamma).  zeta_{Q(i)} = zeta L(s, chi_-4), zeta_{Q(sqrt-3)} = zeta L(s, chi_-3).
Paper's own explicit formula, Thm B.4.4.7 (p. 623): the same with every global term divided by vol,
and the zero sum over all zeros.

Zeros: zeta, L(s, chi_-4), L(s, chi_-3) zeros on the critical line with 0 < gamma <= 1000 (zeros.py), counts
checked against N(T); |h(gamma)| < 1e-14 beyond that height for these test functions.
Columns: LHS/RHS printed (any alpha > 0: every zero up to height 1000 is on the line), RHS/LHS, the printed
identity with I = {0} (empty interior, no zeros), and the residual with every global term divided by vol
as in the paper's own Theorem B.4.4.7 (p. 623).
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classical import TestFunction, arch, load, primes_upto  # noqa: E402
from complexplace import e2  # noqa: E402
from realplace import R  # noqa: E402

TESTS = [
    TestFunction("B: bump(x-1.5)          [|g|>1]", 1.5, 1.0),
    TestFunction("G: (x-1.3) bump((x-1)/0.9) [|g|>1]", 1.0, 0.9, lambda x: x - 1.3),
    TestFunction("A: x(x-0.4) bump((x-0.4)/2)", 0.4, 2.0, lambda x: x * (x - 0.4)),
    TestFunction("C: x^2 bump(x/2.5)", 0.0, 2.5, lambda x: x * x),
    TestFunction("D: x bump((x+0.7)/1.5)", -0.7, 1.5, lambda x: x),
    TestFunction("E: x e^x bump((x-0.2)/1.8)", 0.2, 1.8, lambda x: x * np.exp(x)),
]


def finite_places(field, N):
    """List of (q_v) for all finite places of the field with q_v <= N."""
    out = []
    for p in primes_upto(N):
        p = int(p)
        if field == "Q":
            out.append(p)
            continue
        # imaginary quadratic field of discriminant -D (D = 4 or 3): chi_{-D}(p) decides the splitting
        D = 4 if field == "Q(i)" else 3
        chi = 0 if p % D == 0 or (D == 4 and p == 2) else (1 if (p % D == 1) else -1)
        if chi == 0:
            out.append(p)                     # ramified: residue field F_p
        elif chi == 1:
            out += [p, p]                     # split: two places, q_v = p
        elif p * p <= N:
            out.append(p * p)                 # inert: q_v = p^2
    return out


def lhs_finite(T, field):
    N = int(math.exp(max(abs(T.lo), abs(T.hi)))) + 2
    tot = 0.0
    for q in finite_places(field, N):
        k = 1
        while q ** k <= N:
            x = k * math.log(q)
            tot += math.log(q) * q ** (-k / 2) * float(T.F(x) + T.F(-x))
            k += 1
    return tot


def lhs_arch(T, field):
    if field == "Q":
        return 0.5 * R(T)
    return (1 / (2 * math.pi)) * 2 * e2(T)   # mu_C = 2 dx dy/|z|^2; e2 uses dx dy/|z|^2


def mellin_terms(T, zeros):
    poles = float((T.h(np.array([-0.5j]))[0] + T.h(np.array([0.5j]))[0]).real)
    g = np.concatenate([zeros, -zeros])
    zsum = float(T.h(g).real.sum())
    return poles, zsum


FIELDS = {
    "Q": dict(vol=1.0, zeros=lambda: load("zeta_zeros.txt")),
    "Q(i)": dict(vol=math.pi / 4, zeros=lambda: np.concatenate([load("zeta_zeros.txt"), load("chi4_zeros.txt")])),
    # Q(sqrt-3): r2 = 1, h = 1, R = 1, |d| = 3, w = 6  ->  2 pi / (sqrt 3 * 6)
    "Q(sqrt-3)": dict(vol=math.pi / (3 * math.sqrt(3)),
                      zeros=lambda: np.concatenate([load("zeta_zeros.txt"), load("chi3_zeros.txt")])),
}


def main():
    print(__doc__.split("\n")[0])
    for field, d in FIELDS.items():
        zeros, vol = d["zeros"](), d["vol"]
        print(f"\n=== F = {field}   vol = mu-bar_1(C_F^1) = {vol:.12f}   zeros used: {len(zeros)} ordinates in (0, 1000]")
        print(f"{'test function':36s} {'LHS printed':>16s} {'RHS printed':>16s} {'LHS-RHS':>11s} "
              f"{'RHS/LHS':>10s} {'LHS-RHS(I={0})':>14s} {'B.4.4.7 resid.':>14s}")
        for T in TESTS:
            lhs = lhs_arch(T, field) + lhs_finite(T, field)
            poles, zsum = mellin_terms(T, zeros)
            rhs = vol * (poles - zsum)          # any alpha > 0: every known zero is on the line
            rhs0 = vol * poles                  # alpha = 0: the interior of the strip is empty
            b447 = lhs - (poles - zsum)          # paper's own B.4.4.7 (global terms divided by vol)
            print(f"{T.name:36s} {lhs:+16.12f} {rhs:+16.12f} {lhs-rhs:+11.2e} {rhs/lhs:10.8f} {lhs-rhs0:+14.2e} {b447:+14.2e}")
    # Control: classical explicit formula in its standard form (independent of the paper's conventions)
    print("\nControl (classical Weil explicit formula in standard form, incl. a test with F(0) != 0):")
    from classical import TESTS as CT, chi4, explicit  # noqa: E402

    def chi3(n):
        return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)

    zz, cz4, cz3 = load("zeta_zeros.txt"), load("chi4_zeros.txt"), load("chi3_zeros.txt")
    for T in CT[1:]:
        r1 = explicit(T, zz, 1, 0, lambda n: 1, True)["diff"]
        r4 = explicit(T, cz4, 4, 1, chi4, False)["diff"]
        r3 = explicit(T, cz3, 3, 1, chi3, False)["diff"]
        print(f"  {T.name:36s} residuals: zeta {r1:+.1e}   L(chi_-4) {r4:+.1e}   L(chi_-3) {r3:+.1e}")


if __name__ == "__main__":
    main()
