# Theorem 12.5.1.1 (Trace Formula): a numerical check against the paper's own definitions

Paper: A. M. Masullo, *The Music of the Primes*, draft of 20 November 2026, pp. 519–520.

## Summary

We specialized both sides of the printed identity to test functions φ(g) = f(|g|_F) that depend only
on the idele norm, and evaluated them numerically. We used the paper's definitions for every measure
and constant, and read every formula from the page images.

1. **F = ℚ: consistent with Weil's explicit formula.** For any I = [−α, α] with α > 0, the two sides
   agree to 10⁻¹³. That is the classical explicit formula for ζ. The constant ½ at the real place,
   1/(2π) at the complex place and log q_v/μ_v(O_v^×) at the finite places are right for the paper's
   μ_v. This depends on how ∫′ is read at the archimedean places; see (b) under Side findings.

2. **F = ℚ(i) and F = ℚ(√−3): the printed identity is false, by a constant factor.**
   - The right side uses M(φ)(c) = ∫_{C_F} φ c dμ̄ with the Tamagawa measure μ̄ (p. 519 and Def. A.3.2.1).
     Under μ̄, C_F¹ has volume μ̄₁(C_F¹) = 2^{r₁}(2π)^{r₂}hR/(|d|^{1/2}w).
   - The local terms on the left side are all taken against canonical measures (Remark B.4.4.2).
   - So RHS = μ̄₁(C_F¹)·LHS for every test function. The ratio is π/4 for ℚ(i) and π/(3√3) for ℚ(√−3).
   - The failure persists:
     - for φ supported in |g| > 1, which satisfies the extra support hypothesis of Thm 11.3.2.1;
     - under every reading of the archimedean denominator and of ∫′ that we tried.
   - With every global term divided by μ̄₁(C_F¹), as in the paper's own Theorem B.4.4.7 (p. 623), the
     identity holds to 10⁻¹³.
   - The proof of 12.5.1.1 (p. 520) quotes Theorem 12.4.1.1 without that theorem's factor
     1/μ̄_{S,1}(J¹_{F,S}). Under the proof's hypotheses (Pic(O_{F,S}) = 1, S containing the ramified
     places), μ̄_{S,1}(J¹_{F,S}) = μ̄₁(C_F¹).
   - μ̄₁(C_ℚ¹) = 1, which is why the error does not show over ℚ. We expect the same failure for every
     F with μ̄₁(C_F¹) ≠ 1, but we tested only these two fields.

3. **I is not quantified in the statement.**
   - The statement on p. 519 does not say what I is. The proof (p. 520) takes I = [−α, α].
   - The results it cites (Def. 10.4.3.7, Thms 10.4.2.8 and 11.3.2.1), and those they rest on
     (Thms 10.2.2.1 and 10.3.2.1), all allow α ≥ 0.
   - For α = 0 the interior of the strip 1/2 + I is empty. The identity then fails already for ℚ, by
     about 10⁻².
   - The introduction (pp. 13, 27, 53) takes α > 0 ("a > 0 arbitrary"). So this is a missing
     quantifier in the statement, not a failure of the intended theorem.

4. **What remains of the theorem.** Divide by μ̄₁(C_F¹) and take α ≥ 1/2. The strip then contains all
   zeros of the completed L-functions, and the theorem is the paper's Theorem B.4.4.7: Weil's explicit
   formula in Connes' form, which is classical.
   - **Fixed α with 0 < α < 1/2.** The theorem adds that zeros with |σ − 1/2| ≥ α contribute nothing,
     i.e. a zero-free band. Since the introduction takes α arbitrarily small, the deduction on pp. 29–30
     gets GRH from this.
   - **Where that content would have to be proved.** The left side does not depend on I, and Thm 11.3.2.1
     gives an I-independent value. So the I-dependence has to disappear in the steps of the proof on
     p. 520:
     - (11.3.1.4), taken from the overview §11.3.1, which says it omits several items;
     - Thm 11.3.2.1, whose hypothesis supp φ ⊂ C_F^{>1} the theorem does not assume;
     - Thm 12.4.1.1, which is about a different, I-free kernel on J_{F,S}, and has no
       M(φ)(|·|^{±1/2}) terms;
     - (12.5.1.4), justified "by design".
   - The section titled for this bridge, §12.3 ("Putting the pieces together: global to S-local
     deconstruction of the geometric kernel", p. 514), has only its heading and the footnote
     "Include section." The proof does not cite it.
   - None of this can be tested numerically: for ℚ, every computed zero is on the critical line.

Points 2 and 3 show that the statement *as printed* is false. Both are repairable. The substance of
the dispute is point 4, which this check neither confirms nor refutes. We do not know whether point 2
is the counterexample reported on Zulip. It is consistent with the README's description, though:
"correcting one of the errors reduces the theorem to a classical formula".

## 1. Definitions used

All read from the page images.

| Object | Definition | Where |
|---|---|---|
| μ_v, real v | dx/\|x\| (= \|·\|_v^{−1}ν_v with ν_v = dx) | Def. A.3.2.2, p. 560; Def. A.3.1.2, p. 556 |
| μ_v, complex v | 2 dx dy/\|z\|² (ν_v = \|dz∧dz̄\|), with \|z\|_v = \|z\|² | Def. A.3.2.2; Def. A.3.1.2; p. 7 |
| μ_v, finite v | Haar, with μ_v(O_v^×) = \|disc(O_v/ℤ_p)\|_p^{1/2} | (A.3.2.2), p. 560 |
| μ̄ on C_F | Tamagawa: dμ̄ = dt/t ⊗ dμ̄₁, μ̄₁(C_F¹) = 2^{r₁}(2π)^{r₂}h_F R_F/(\|d_F\|^{1/2}w_F) | Def. A.3.2.1, (A.3.2.1), pp. 559–560; §A.3.3.1, p. 561 |
| M(φ)(c) | ∫_{C_F} φ c dμ̄ | Thm 12.5.1.1, p. 519 |
| canonical measures | ½μ_v (real), μ_v/(2π) (complex), (log q_v/μ_v(O^×))μ_v (finite); μ̄/μ̄₁(C_F¹) on C_F | Remark B.4.4.2, p. 620; (D.2.1.6), p. 794 |
| ∫′ | Connes–Weil principal value for the standard character (Def. B.4.4.1, p. 619). With φ(1) = 0: at finite places the ordinary integral (Lemma B.4.4.3); at archimedean places the ordinary integral if Prop. B.4.3.2 is used as printed, see (b) under Side findings | pp. 619–621 |
| L(c) | completed Hecke L-function, L(χ\|·\|^s) = L(s, χ) | p. 152 (from the text extraction) |

At complex places the theorem prints |1 − g_v| with no subscript. Only the normalized
|1 − g_v|_v = |1 − g_v|² agrees with the definition (B.4.4.2) and gives the classical term. With the
ordinary modulus, the ratio to the classical term varies with the test function, and the identity
still fails.

## 2. Specialization

Take φ(g) = f(|g|_F) with f ∈ C_c^∞(ℝ_{>0}) and f(1) = 0. Then φ ∈ 𝒟(C_F): it is invariant under every
local unit group, and it is compactly supported because C_F¹ is compact. Also φ(C_F¹) = 0. The author
uses such norm-composed test functions himself (§1.4.2, p. 28). Write F(x) = f(eˣ) and
h(r) = ∫ F(x)e^{irx} dx.

- **Finite place v.** Exact, since f(1) = 0. On the shell |g|_v = q_v^{−m} one has
  |1 − g|_v = max(1, q_v^{−m}), and each shell has measure μ_v(O_v^×). So the term is
  log q_v Σ_{k≥1} q_v^{−k/2}(F(k log q_v) + F(−k log q_v)). Summed over v, this is the von Mangoldt sum
  for ζ_F.
- **Real place.** ½ ∫ F(u) e^{u/2} (1/|1 − eᵘ| + 1/(1 + eᵘ)) du.
- **Complex place.** ∫ F(u) e^{u/2}/|1 − eᵘ| du. This is (1/2π) of the paper's integral; it uses
  that the angular average of |1 − re^{iθ}|^{−2} is 1/|1 − r²|.
- **Global terms.** M(φ)(|·|^s) = μ̄₁(C_F¹) ∫ f(t) t^s dt/t. Every character that is nontrivial on
  C_F¹ integrates φ to 0, so only the zeros of ζ_F contribute. A zero ½ + iγ contributes μ̄₁(C_F¹)·h(γ).
  Here ζ_{ℚ(i)} = ζ·L(·, χ₋₄) and ζ_{ℚ(√−3)} = ζ·L(·, χ₋₃).

## 3. Numbers

Test functions: F(x) = p(x)·bump((x − c)/a), with bump(y) = exp(−1/(1 − y²)) and F(0) = 0. B and G
are supported in |g| > 1.

Zeros: every critical-line zero with 0 < γ ≤ 1000. For ζ that is 649, which is complete
(N(1000) = 649). For L(s, χ₋₄) it is 868, and the Riemann–von Mangoldt main term is 868.4. For
L(s, χ₋₃) it is 823, against a main term of 822.6. Beyond γ = 1000, |h(γ)| < 10⁻¹⁴.

Columns:
- "LHS − RHS": the printed identity for any α > 0.
- "RHS/LHS": the ratio of the two sides.
- "I = {0}": the printed identity with α = 0.
- "B.4.4.7": the residual after dividing every global term by μ̄₁(C_F¹).

```
F = Q          vol = 1                  LHS printed   RHS printed    LHS-RHS    RHS/LHS   LHS-RHS(I={0})  B.4.4.7
B: bump(x-1.5)               [|g|>1]   +1.162721242  +1.162721242   +1.2e-13  1.00000000     -9.8e-03     +1.2e-13
G: (x-1.3) bump((x-1)/0.9)   [|g|>1]   -0.233718229  -0.233718229   +3.1e-15  1.00000000     +1.4e-02     +3.1e-15
A: x(x-0.4) bump((x-0.4)/2)            +1.400281978  +1.400281978   +1.4e-13  1.00000000     +7.7e-03     +1.4e-13
C: x^2 bump(x/2.5)                     +2.812349007  +2.812349007   -4.3e-13  1.00000000     +1.1e-02     -4.3e-13
D: x bump((x+0.7)/1.5)                 -1.118537703  -1.118537703   +3.4e-14  1.00000000     +3.6e-03     +3.4e-14
E: x e^x bump((x-0.2)/1.8)             +2.050440970  +2.050440970   +2.2e-13  1.00000000     +3.2e-02     +2.2e-13

F = Q(i)       vol = pi/4 = 0.785398
B: bump(x-1.5)               [|g|>1]   +1.130904795  +0.888210549   +2.4e-01  0.78539816     +2.1e-01     +1.3e-13
G: (x-1.3) bump((x-1)/0.9)   [|g|>1]   -0.229733310  -0.180432120   -4.9e-02  0.78539816     -3.5e-02     +7.6e-15
A: x(x-0.4) bump((x-0.4)/2)            +1.389972423  +1.091681788   +3.0e-01  0.78539816     +3.0e-01     +1.6e-13
C: x^2 bump(x/2.5)                     +2.864556454  +2.249817378   +6.1e-01  0.78539816     +6.6e-01     -5.2e-13
D: x bump((x+0.7)/1.5)                 -1.191121530  -0.935504662   -2.6e-01  0.78539816     -3.1e-01     +3.4e-14
E: x e^x bump((x-0.2)/1.8)             +2.184777464  +1.715920208   +4.7e-01  0.78539816     +6.0e-01     +2.2e-13

F = Q(sqrt-3)  vol = pi/(3 sqrt 3) = 0.604600
B: bump(x-1.5)               [|g|>1]   +1.213713262  +0.733810781   +4.8e-01  0.60459979     +5.0e-01     +1.3e-13
G: (x-1.3) bump((x-1)/0.9)   [|g|>1]   -0.241706511  -0.146135705   -9.6e-02  0.60459979     -9.2e-02     +9.4e-15
A: x(x-0.4) bump((x-0.4)/2)            +1.413069042  +0.854341243   +5.6e-01  0.60459979     +5.7e-01     +1.4e-13
C: x^2 bump(x/2.5)                     +2.865931347  +1.732741485   +1.1e+00  0.60459979     +1.2e+00     -5.0e-13
D: x bump((x+0.7)/1.5)                 -1.076029540  -0.650567232   -4.3e-01  0.60459979     -4.0e-01     +2.6e-14
E: x e^x bump((x-0.2)/1.8)             +1.975422297  +1.194339902   +7.8e-01  0.60459979     +7.6e-01     +2.1e-13
```

**Control.** We also evaluated the classical explicit formula for ζ, L(s, χ₋₄) and L(s, χ₋₃) in its
standard form, with the same zeros and quadrature but independently of the paper's conventions. It
balances to between 10⁻¹³ and 10⁻¹⁵, including a test function with F(0) ≠ 0.

## 4. Side findings

These are separate from Thm 12.5.1.1 and could go in ERRATA.

- **(a) Proof of B.4.4.7, p. 625.** The proof defines D(φ) with dμ̄ and no 1/μ̄₁(C_F¹), which
  contradicts the statement (B.4.4.16) on p. 623. The statement is the correct version: it balances to
  10⁻¹³ above.
- **(b) Prop. B.4.3.2, (B.4.3.13)–(B.4.3.14), p. 619.** As printed, these are exactly M_v times the true
  pairing ⟨Ĵ_v, f⟩ = ∫ log|x|_v f̂ dν_v. The factor is M_v = 2 (real) and 2π (complex), checked with the
  self-dual Gaussians.
  - If Def. B.4.4.1 is read literally with the true Fourier transform, the archimedean ∫′ is (1/M_v)
    times the ordinary integral. The identity then fails even for ℚ.
  - The proof of B.4.4.7 ((B.4.4.20)) uses the printed formulas, which give the ordinary integral.
    That is the reading used above.
  - Footnote 16 on p. 619 flags "factors of 2 and 2π".
- **(c) Note after Thm 12.4.1.1, p. 514.** It calls μ_v canonical at real places. By Remark B.4.4.2
  the canonical measure is ½μ_v. The formula itself has the ½ and is fine.
- **(d) Introduction, p. 28.** The strip is printed as σ(c) ∈ I, not 1/2 + I. The real-place term there
  uses |1 − t| for both signs.

## 5. Not checked

- **Test functions and fields.** Only norm-only φ. We did not test characters nontrivial on C_F¹,
  ramified characters or fields other than ℚ, ℚ(i) and ℚ(√−3). The failure for every field with
  μ̄₁(C_F¹) ≠ 1 is a prediction.
- **Proofs.** The proofs of 10.4.2.8, 11.3.2.1 and 12.4.1.1, and the split (11.3.1.4), were not
  followed line by line. Point 4 lists only the gaps that are visible on the surface.
- **Unquantified parameters.** The Sobolev parameter s is not quantified in 12.5.1.1. The results it
  depends on require s > 1, s ≥ 3 or s ≥ 0, and the introduction uses s = 0. Neither is the character χ
  in the proof's choice of S. Dependence on s was not examined.
- **Text-only readings.** Three things were read from the text extraction only: that L(c) is the
  completed L-function (pp. 12, 152), the definition of σ(c) (3.1.1.2), and the closed-strip form of the
  spectral trace (10.3.2.3)–(10.3.2.4) on p. 454.
- **The range 0 < α < 1/2.** This is GRH-type and not testable numerically.
- **Zero data.** The numbers assume that all zeros up to height 1000 are simple and on the line. That
  is known for ζ, and the counts and the B.4.4.7 residuals support it for the two L-functions.
- **Zulip.** The thread was not readable from here.

## Reproducing

Python 3.11 with `mpmath`, `numpy` and `scipy`. Run in this folder:

```sh
python -I -B zeros.py 1000   # zeros of zeta, L(s,chi_-4), L(s,chi_-3) up to height 1000 (~40 min)
python -I -B check.py        # the table and the control (~3 min)
```

`classical.py` is the classical explicit formula used as the control. `realplace.py` and
`complexplace.py` compare the archimedean terms with the classical ones. The ratio is exactly 1 for
the paper's μ_v, and it varies with the test function if |1 − z| is read as the ordinary modulus.

This check was carried out with an AI assistant (Claude): five agents extracted definitions from the
page images, two re-derived the result blind with their own code, and three tried to refute the
claims and audited the code. No specialist in the area has reviewed it.
