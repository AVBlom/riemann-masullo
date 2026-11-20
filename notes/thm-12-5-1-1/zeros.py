"""Compute and cache ordinates of nontrivial zeros on the critical line.

  zeta:          gamma_n = Im of the n-th zero of zeta(s), via mpmath.zetazero (n = 1..N).
  L(s, chi_-q):  for the odd real primitive characters chi_-4 (mod 4) and chi_-3 (mod 3):
                 sign changes of the real Hardy function
                   Z(t) = exp(i*theta(t)) * L(1/2+it, chi),
                   theta(t) = (t/2) log(q/pi) + Im loggamma((1/2 + 1 + it)/2)   (a = 1, root number 1),
                 on a grid of step 0.04, refined with mpmath.findroot.

Usage: python -I zeros.py T_MAX [zeta] [chi4] [chi3]   (default: all; writes <name>_zeros.txt here)
"""
import math
import os
import sys
import time

import mpmath as mp

mp.mp.dps = 20
HERE = os.path.dirname(os.path.abspath(__file__))

CHARS = {"chi4": (4, [0, 1, 0, -1]), "chi3": (3, [0, 1, -1])}


def zeta_zeros(T):
    out, n = [], 1
    while True:
        g = float(mp.im(mp.zetazero(n)))
        if g > T:
            return out
        out.append(g)
        n += 1


def hardy_Z(t, q, chi):
    t = mp.mpf(t)
    s = mp.mpc(0.5, t)
    theta = (t / 2) * mp.log(q / mp.pi) + mp.im(mp.loggamma((s + 1) / 2))
    val = mp.exp(1j * theta) * mp.dirichlet(s, chi)
    return val.real, abs(val.imag)


def l_zeros(T, q, chi, step=0.04):
    out, max_imag = [], 0.0
    t0, z0 = 0.0, hardy_Z(0.0, q, chi)[0]
    t = step
    while t <= T + step:
        z, im = hardy_Z(t, q, chi)
        max_imag = max(max_imag, im)
        if z == 0 or (z > 0) != (z0 > 0):
            out.append(float(mp.findroot(lambda u: hardy_Z(u, q, chi)[0], (t0, t), solver="anderson")))
        t0, z0 = t, z
        t += step
    return out, max_imag


def write(name, zs):
    with open(os.path.join(HERE, f"{name}_zeros.txt"), "w") as fh:
        fh.write("\n".join(repr(g) for g in zs) + "\n")


def main():
    T = float(sys.argv[1]) if len(sys.argv) > 1 else 1000.0
    which = sys.argv[2:] or ["zeta", *CHARS]
    if "zeta" in which:
        t0 = time.time()
        zz = zeta_zeros(T)
        write("zeta", zz)
        print(f"zeta: {len(zz)} zeros up to {T} in {time.time()-t0:.0f}s; "
              f"Riemann-von Mangoldt main term {T/(2*math.pi)*math.log(T/(2*math.pi*math.e))+7/8:.1f}", flush=True)
    for name, (q, chi) in CHARS.items():
        if name not in which:
            continue
        t0 = time.time()
        zs, mi = l_zeros(T, q, chi)
        write(name, zs)
        main_term = T / (2 * math.pi) * math.log(q * T / (2 * math.pi * math.e))
        print(f"{name}: {len(zs)} zeros in (0, {T}] in {time.time()-t0:.0f}s; main term (T/2pi)log(qT/2pi e) = "
              f"{main_term:.1f}; max |Im Z| on grid = {mi:.2e}; first zeros {zs[:3]}", flush=True)


if __name__ == "__main__":
    main()
