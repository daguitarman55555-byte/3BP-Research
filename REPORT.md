# Fifth-stage research: the hidden finite cover for signed-area histories

**Setting.** Planar Newtonian three-body problem, known G (=1 here) and known pairwise-distinct
positive masses, observable the signed area A(t) = ½ det(r₂−r₁, r₃−r₁), states taken modulo
translation, Galilean boost and proper rotation (not reflection).

**Target.** The degree d = [L : K(h₇)], where L is the state field, K = ℝ(A,…,A⁽⁶⁾), h₇ = A⁽⁷⁾,
N = [L:K], n = [K(h₇):K], N = nd. Equivalently: does the complete area-history relation
𝓡 = {(X,Y): A_X⁽ᵏ⁾ = A_Y⁽ᵏ⁾ ∀k} have a non-diagonal 7-dimensional component?

No web, literature, saved memory or other conversation was used. The three supplied packages
were treated only as claims to audit.

Status labels: **PROVED**, **EXACT COMPUTATIONAL CERTIFICATE**, **STRONGLY SUPPORTED NUMERICALLY**,
**CONJECTURE**, **FALSE**, **UNRESOLVED**.

---

## DECISIVE RESULT

**C. Partial result only. It is the strongest result here, and it points strongly toward d = 1.**

1. **EXACT COMPUTATIONAL CERTIFICATE (new).** For masses (1,2,3) there is an explicit complex
   sixth-jet target c₀ (stored as exact binary64 numbers) with **24 148 distinct, nondegenerate
   complex points of the lifted state variety V** in F₆⁻¹(c₀). Each point lies in a 212-bit
   ball-arithmetic Krawczyk box of radius < 1.6·10⁻⁵¹ that contains exactly one solution. The 24 148
   enclosures of A⁽⁷⁾ are **pairwise disjoint**: the smallest gap is 1.79, while every
   enclosure radius is below 10⁻²¹. Consequences, which are **PROVED** from the certificate:

   * N ≥ 24 148 and **n ≥ 24 148**. Every polynomial seventh-order area relation has degree
     ≥ 24 148 in A⁽⁷⁾. The fourth stage's bound was 37.
   * **d ≤ N / 24 148.** So d ≥ 2 forces N ≥ 48 296: a nontrivial cover would need at least as
     many uncertified sheets as certified ones.
   * Both bounds hold for generic masses too, because the certified roots persist under small
     changes of mass.
2. **STRONGLY SUPPORTED NUMERICALLY.** The certified set is about 98.8 % of the whole fibre.
   Random points of V, tracked back to c₀ one path at a time, landed in the known set
   2 862 times out of 2 898 (and 1 893 of 1 918 in an earlier test). The misses are genuine
   new sheets, and they repeat. This estimates **N ≈ 24.4·10³ ± a few hundred**. Since d ≥ 2
   needs at least 50 % of the sheets to be missing, the observed ≈1.2 % gives
   **n = N and d = 1** (h₇ is a primitive element: L = K(A⁽⁷⁾)) with overwhelming numerical
   weight.
3. **PROVED (reduction).** If c₀ is not an asymptotic critical value of F₆ (F₆ is proper over
   a neighbourhood of c₀) and F₆⁻¹(c₀) has exactly as many points as are certified, then d = 1
   exactly, and generic complete-area-history uniqueness follows. The only missing ingredient for
   option A is therefore a proof that a *single* sixth-jet fibre is complete. Monodromy cannot
   supply that proof. An exact 𝔽ₚ Gröbner computation (msolve) of the same fibre was attempted
   and is infeasible at this size; §6 gives the details.

Item 2 holds for two unrelated mass triples, (1,2,3) and (10,17,29). The same pipeline does
detect the genuine d = 2 cover at equal masses (1,2,2), which serves as a positive control.

Option A is not claimed, because fibre completeness is not proved. Option B is
**excluded within the certified part of the fibre**: no two of the 24 148 certified analytic
sheets are glued even at order seven.

---

## 1. Phase 1 — audit of the foundations actually used

| Claim | Check performed here | Verdict |
|---|---|---|
| Reduced state has dimension 7 | 12 − 2 (translation) − 2 (boost) − 1 (rotation) | **PROVED** (correct) |
| Generic full rank of (A,…,A⁽⁶⁾) | New exact engine (`jets.py`, FLINT polynomial derivation). Exact rational Jacobian determinant on V at Heron 13-14-15 states, masses (1,2,3) and (3,5,7): det = 7.58·10⁻⁷ and 1.03·10⁻¹¹, exactly nonzero | **EXACT CERTIFICATE**, confirms prior claim |
| Engine validity | Taylor jets match direct high-accuracy ODE integration; remainder scales as h⁹ (2.6e-10 → 4.8e-13 when h halves) | confirmed |
| Differential-field closure at order 7 (h_j ∈ ℝ(h₀..h₇)) | Proof re-derived: D(K) ⊂ K(h₇); differentiating p(h₇)=0 gives h₈ ∈ K(h₇) | **PROVED** (correct). It does not imply d=1 (fourth-stage toy z₇² example is valid) |
| Multi-branch sixth-jet fibres | Our full fibre, moved to the third-stage target state (3,0,4,1,2,3,−1), see §9 | see §9 |
| Local 10th/11th-order equilateral classification | Independent mpmath Taylor engine (60 digits), turning pair λ=0, ω=±1, total mass 9, masses ∝(1,2,3): jets of the partners agree through order 11 (differences ~1e-52); pair Jacobian singular values for rows 0..10 are all ≥1.1e-3 (rank 11), rows 0..9 rank 10, rows 0..11 has 12th singular value 1.5e-38 (rank stays 11) | numerically reproduced at this mass; the universal-mass claim was not re-certified |

No foundational claim failed.

## 2. The algebraic set-up used for exact work

V ⊂ ℂ¹⁰ has coordinates (l,x,y,u,v,w,z,h_a,h_b,h_c), with a = (l,0), b = (x,y), inertial relative
velocities (u,v),(w,z), and the equations h_a l = 1, h_b²(x²+y²) = 1, h_c²((x−l)²+y²) = 1.
V is an irreducible smooth affine variety, and its function field is the state field L of the
fourth stage. Fixing h_a = 1/l, rather than h_a² l² = 1, removes the rotation-by-π gauge
component. The area jets A⁽ᵏ⁾ are explicit polynomials on V (sizes for (1,2,3): 1, 3, 5, 24,
99, 437, 1508, 5116, 15017 terms for k = 0..8). The sixth-jet fibre over c is the square
system {A⁽ᵏ⁾ = c_k (k ≤ 6)} ∪ {the three equations of V}: 10 equations in 10 unknowns, total
degree up to 21.

N counts **all** complex points of V, including complex states and points with "wrong-sign"
inverse distances. Physical states are a Zariski-dense real subset, so d = 1 algebraically
implies generic physical uniqueness. Conversely, d > 1 need not give a *physical* twin.

## 3. Phase 3 — the differential-Galois reformulation (PROVED)

**Lemma 3.1.** Let 𝒜 = K(h₇). It is D-stable and equals the whole area differential field.
Let M be a normal closure of L/𝒜. D extends uniquely to M, and every τ ∈ Gal(M/𝒜) commutes
with D, because τDτ⁻¹ is another extension of D|_𝒜. The d embeddings L → M over 𝒜 are
therefore **differential** embeddings. Geometrically, d > 1 holds exactly when there is a
non-identity algebraic correspondence Ψ on V that commutes with the reduced flow and preserves
every A⁽ᵏ⁾. Equivalently, 𝓡 has a non-diagonal 7-dimensional component, and through every
F₆-regular point X that component meets the diagonal only off (X,X).

**Lemma 3.2 (constants).** Let C(·) denote field constants. A Wronskian argument gives
linear disjointness: C(L) and 𝒜 are linearly disjoint over C(𝒜). Hence
[C(L)·𝒜 : 𝒜] = [C(L):C(𝒜)] ≤ d. Moreover C(M) is algebraic over C(𝒜), because the minimal
polynomial of a constant has constant coefficients. Consequences:

* if d = 1, then H, J (and every algebraic first integral) are rational functions of
  A,…,A⁽⁷⁾;
* any hidden branch Ψ acts on (H,J) by an algebraic correspondence H∘Ψ = φ(H,J),
  J∘Ψ = ψ(H,J);
* if H or J were not determined by the area history, then d > 1.

**Lemma 3.3 (block counting).** If d ≥ 2, the F₇-classes are blocks of size d for the monodromy
group of F₆. A set of sheets with pairwise distinct A⁽⁷⁾ meets each block at most once, so it
contains at most N/d sheets. Hence N ≥ d·24 148. This proves the inequality in the Decisive
Result.

**Primitive-element attack.** The structure above reduces "L = K(h₇)" to the statement that h₇
separates one complete generic fibre of F₆. The certificate (§5) proves separation on 24 148
sheets. What is missing is completeness of the fibre, not any failure of primitivity. No hidden
element of L that is invariant under all area derivatives was found; such an element would have
to take equal values on two sheets with equal A⁽⁷⁾, and none exist among the certified sheets.

## 4. Phase 2 — computing the fibre (methods)

* `monodromy.py`, `orbit.py`, `orbit_extend*.py`, `orbit_grow.py`: batched, numba-compiled
  predictor–corrector (RK4 plus Newton) parameter homotopy on the 10×10 system, run as monodromy
  orbit closure under random triangle loops c₀ → p → q → c₀ in ℂ⁷.
* **Important negative lesson.** A first orbit closure under 6 fixed loops *stopped* at
  19 893 sheets: the queue emptied and no new points appeared. That plateau was **false**. Fresh
  loops produced more sheets: 22 542, then 23 796, then 24 123. Closure under one finite loop set
  is not completeness. Path failures (≈2–8 %) leave holes in the permutation data.
* `membership.py`: an independent completeness estimator. Sample a random point Y′ ∈ V (random
  complex state at several scales, both sign branches of h_b and h_c), track the single path
  from F₆(Y′) to c₀ through a random intermediate target, and test whether its endpoint is
  already known. Misses are new sheets.
* `cert_arb.py` / `run_arb*.py`: rigorous certification in python-flint `acb` ball arithmetic.
  Each point is polished by Newton's method at 212 bits, then given the Krawczyk test
  K(X) = x̃ − Y F(x̃) + (I − Y J(X))(X − x̃) ⊂ int X on a box X. J(X) is enclosed by interval
  evaluation of the exact Jacobian polynomials. Inclusion proves that X holds exactly one root
  and that J is invertible there. A⁽⁷⁾ is then enclosed on X.

## 5. Phase 2 — results for masses (1,2,3)

| Quantity | Value | Status |
|---|---|---|
| Certified distinct regular points of F₆⁻¹(c₀) | **24 148** (all at 212 bits, no box overlaps) | **EXACT CERTIFICATE** |
| Pairwise-disjoint A⁽⁷⁾ enclosures | 24 148 / 24 148; min gap 1.79, max radius 9.3·10⁻²²; min relative gap 1.3·10⁻⁴ | **EXACT CERTIFICATE** |
| N ≥ 24 148, n ≥ 24 148, d ≤ N/24 148 | consequences of the two rows above | **PROVED** |
| Membership tests on the final set | 2 862 / 2 898 tracked paths land in the set (1.24 % new); earlier 1 893 / 1 918 (1.3 %) and 193 / 196 | numerical |
| Estimated N | ≈ 24 450 (capture–recapture with repeated misses) | **STRONGLY SUPPORTED NUMERICALLY** |
| d = 1, n = N | the miss fraction is ≈1.2 %, while d ≥ 2 requires ≥ 50 % | **STRONGLY SUPPORTED NUMERICALLY** |
| Real points at the complex target | 0 (as expected) | — |

Files: `certificates/fibre_123_final_summary.json`, `certificates/arb_fibre_123_final.pkl`
(all centres, radii, A⁽⁷⁾ midpoints and radii, and precisions), `data_orbit_123_ext4.npz`,
`certificates/membership_*.json`.

**Why the miss statistics bear on d.** Lemma 3.3 shows that d ≥ 2 forces at least half of all
sheets to lie outside the certified set. Random points of V land on sheets with a non-uniform,
but not adversarial, distribution. Observing ≲ 1.3 % misses in about 5 000 trials is
incompatible with a missing mass of ≥ 50 % unless the landing measure almost entirely avoids an
entire half of the sheets. No mechanism for that is known. The landing distribution has not been
analysed rigorously, which is why this is numerical support and not proof.

## 6. Attempted exact counting (Gröbner)

`make_msolve.py` writes the same square system over 𝔽_p (p = 1073741827) at a uniformly random
target. msolve 0.6.5 (F4, one thread) reached step degree 14, with Macaulay matrices of
658 237 × 741 113, and did not finish in 25 minutes (`logs/msolve_p1.log`). With N ≈ 2.4·10⁴
solutions and degree-21 equations in 10 variables, a direct modular Gröbner basis is beyond
this session's budget. Even a finished modular count would give only N(c) ≤ N_generic at that
target, plus a lucky-prime caveat. Upper bounds on N need a properness or compactification
argument, which §10 discusses.

## 7. Phases 4 & 5 — pair relation and hidden symmetries

* **Pair relation (PROVED, mostly prior).** The equilateral opposite-spin manifold is an exact
  3-dimensional component of 𝓡 near the classified pairs, and it is lower-dimensional. Any
  7-dimensional non-diagonal component Z misses (X,X) for every F₆-regular X: near such a point
  𝓡 ⊂ 𝓡₆ is locally the diagonal. Combined with §5, **Z contains no pair (Yᵢ,Yⱼ) of the 24 148
  certified sheets over c₀, nor any pair of their analytic continuations near c₀.** If Z exists,
  its fibre over each certified sheet consists of uncertified sheets or escapes to the boundary
  of V.
* **Direct substitution of candidate symmetries** (`symmetry_tests.py`, three random states).
  The first jet order at which A_X and A_ΨX differ:
  * spin flip (shape velocities kept, J→−J): order 2;
  * time reversal: order 1;
  * shape-velocity flip: order 1;
  * radial (dilation) flip: order 1;
  * reflection, and reflection with time reversal: order 0 (A changes sign).

  None is a hidden symmetry for unequal masses.
* **Nonlinear, algebraic, multi-valued and open-dense correspondences** are all subsumed by
  Lemma 3.1. Any of them would show up as two sheets of a generic sixth-jet fibre with equal
  A⁽⁷⁾. On the certified 98.8 % of the fibre there are none. This covers rational maps,
  algebraic involutions, shape–scale exchanges, maps involving integrals, and state-dependent
  time-reversal maps alike, with no need to guess a functional form.
* **Positive control.** For masses (1,2,2), the pipeline recovers the known swap-and-reflect
  cover at about 90 % of sheets (§11). The method is not blind to covers.

## 8. Phase 6 — integrals of motion

Over the 24 148-sheet fibre (`invariants.py`), the following each take **pairwise distinct
values on all sheets** (rounded at 10⁻⁶; r₁₂ shows a single rounding collision, 24 147 values):

H, J, I, İ, r₁₂², r₁₃², r₂₃², r₁₂, r₁₃, r₂₃.

Each of them is therefore (numerically) a primitive element of L/K, with degree ≈ N over the
sixth-jet field: none is a function of the sixth jet. Combined with d = 1 and Lemma 3.2, each is
a rational function of A,…,A⁽⁷⁾ (**STRONGLY SUPPORTED NUMERICALLY**; the explicit rational
expressions have degree ≳ 2.4·10⁴ and were not computed). No invariant was found with two
algebraic branches over the area differential field.

## 9. Phases 7–9 — singular output geometry, branches, and other regimes

**Self-intersection versus cover.** The seventh-output hypersurface P(A,…,A⁽⁷⁾) = 0 has degree
n ≥ 24 148 in A⁽⁷⁾. Above a generic sixth jet its 24 148+ sheets have distinct A⁽⁷⁾. Its
singular locus is where two sheets' A⁽⁷⁾ values coincide: a discriminant hypersurface in
sixth-jet space. That is **self-intersection of a finite-jet image**. The fourth stage's
transverse crossing at (1,3,5) and the third stage's 9th-order false pairs are of this kind,
and they separate at later orders. A **genuine complete-history cover** would need two sheets
with *identically* equal A⁽⁷⁾ near a generic target. The certificate excludes this for every
pair of certified sheets, because identically equal functions would be equal at c₀.

**Regimes away from the equilateral family.** The whole fibre was transported by parameter
homotopy to real physical targets (`regimes.py`):

| Regime (masses 1,2,3) | distinct endpoints | failed paths | real physical sheets | true source recovered | min rel. A⁽⁷⁾ gap (all sheets) | nearest other sheet's A⁽⁷⁾, rel. to source |
|---|---|---|---|---|---|---|
| strongly scalene | 23 965 | 182 | 69 | yes | 6.5·10⁻⁷ | 1.7·10⁻² |
| large angular momentum | 23 969 | 179 | 68 | yes | 5.2·10⁻⁸ | 1.2·10⁻² |
| zero angular momentum | 23 950 | 198 | 72 | yes | 1.2·10⁻⁴ | 7.4·10⁻⁴ |
| strongly expanding | 23 136 | 1 011 | 18 | yes | 1.6·10⁻⁷ | 4.1·10⁻⁶ |

In every regime, no two transported sheets share A⁽⁷⁾: the gaps are small, but they lie far
above the ~10⁻¹³ coincidence level that the equal-mass positive control produces (§11). At real
targets there are **18–72 real physical sixth-jet preimages**, all separated by A⁽⁷⁾.
(**STRONGLY SUPPORTED NUMERICALLY**: these are double-precision transports and were not
re-certified. 1–4 % of paths failed, and coincident endpoints were merged.)

**Not completed.** The near-collision (r₁₂ = 0.06), close-binary hierarchical, strongly
contracting and random regimes were started. The first ran for more than 40 minutes on stiff
paths and was stopped. Milder versions (r₁₂ = 0.2, and a binary of 0.25 with the third body at
≈3.2) were lost to a worker restart. Phase 9 is therefore only partly covered. The fibre
transport is the right tool for it; it is simply slow near collision.

**Audit of the third-stage 37 real roots** (`audit_37b.py`). We transported the fibre to the
third-stage target state (3,0,4,1,2,3,−1), using the *exact rational* sixth jet as the target.
Of the paths, 21 966 were tracked and 2 182 failed; real targets are harder. The transport
recovered 31 of the 37 roots directly and found 48 real physical endpoints. We then certified
the union of our 48 and their 37 at 212 bits against the exact target. The result:
**54 distinct real physical regular preimages, all 37 third-stage roots among them, with pairwise
disjoint A⁽⁷⁾ enclosures (EXACT COMPUTATIONAL CERTIFICATE).** The third-stage count of 37 is
confirmed as correct but incomplete; the real physical sixth-jet fibre there has at least 54
points. (A uniquely contained root of a real system in a conjugation-symmetric box is real.)

## 10. Phase 10 — generic versus uniform theorems

* **Uniform theorem.**
  * Existence of a finite k₀ (equality through k₀ implies equality of the complete history for
    all collision-free states, all masses): **PROVED**. This is Noetherianity, as in the third
    stage, and the argument is correct.
  * k₀ ≥ 10 for every fixed distinct-mass triple: **PROVED** in the fourth stage; the audit
    found no flaw, and our independent equilateral ranks agree.
  * The exact value of k₀: **UNRESOLVED**.
* **Generic theorem.** On a Zariski-open dense set, complete-area-history injectivity is
  equivalent to d = 1, and order 7 then suffices.
  * Equivalence of injectivity with d = 1: **PROVED**.
  * d = 1, hence generic injectivity at order 7: **STRONGLY SUPPORTED NUMERICALLY**.
  * Order 6 is generically insufficient: **PROVED**. There are many real physical sixth-jet
    preimages on open sets (third stage; also §9).

  So the generic order (7) and the worst-case order (≥ 10) differ, as the brief anticipated.
* **Remaining obstruction, stated precisely.** Prove that one generic sixth-jet fibre is
  complete. One route is to show that some explicit target is not an asymptotic critical value
  of F₆ (a Jelonek-set computation) and to count its fibre exactly. Another is any exact method
  giving N ≤ 48 295.

## 11. Phase 11 — consequences if d = 1, which numerics strongly support

1. Order seven is generically globally sufficient. This is equivalent to d = 1.
2. Generic sixth-jet degree for (1,2,3): 24 148 ≤ N, numerically about 24.45·10³.
   For a second, unrelated mass triple (10,17,29), the whole fibre was transported by a homotopy
   in (masses, target) space through complex masses (`mass_transport.py`). It gave 23 849
   distinct endpoints (299 failures, 0 duplicates); the source state was recovered; all A⁽⁷⁾ were
   distinct (min rel. gap 7.1·10⁻⁷); random-point membership was 749/763 (98.2 %). The fibre
   size is therefore the same, consistent with N being mass-generic. **STRONGLY SUPPORTED
   NUMERICALLY.**
3. Real physical sixth-jet branches vary by region of target space. At the sampled real targets
   the counts were 69, 68, 72 and 18 (§9), plus the third-stage target below. So the 37 of the third stage were a subset at one target.
4. Smallest generic sufficient order: **7**, conditional on d = 1. Order 6 fails on open sets.
5. **Practical algorithm.** Store one complete generic fibre (≈ 24·10³ points; here
   `data_orbit_123_ext4.npz`). For observed (A,…,A⁽⁷⁾), track all stored points by parameter
   homotopy to the observed sixth jet (minutes on 4 cores). Keep the real physical endpoints and
   select the unique one matching A⁽⁷⁾. The source state was recovered in every regime test. The A⁽⁷⁾
   separation from the nearest other sheet ranged from 1.7·10⁻² down to 4.1·10⁻⁶ (relative), so the
   final selection step needs A⁽⁷⁾ to about 6 significant digits in the worst case seen.

**Positive control: does the method detect a real cover?** For masses (1,2,2), swapping bodies 2
and 3 and reflecting is an exact signed-area-history symmetry, so d ≥ 2. The same fibre was
transported from (1,2,3) to (1,2,2) (`mass_transport_122.py`, `control_122.py`). It gave 23 823
endpoints, min relative A⁽⁷⁾ gap **3.8·10⁻¹⁶**, and **10 658 A⁽⁷⁾-coincident pairs**. All 2 000
pairs checked are exactly the swap-and-reflect images of each other. So when a genuine cover
exists, the pipeline finds it in roughly 90 % of sheets. The absence of any coincidence for
(1,2,3) and (10,17,29) is a real negative.

## 12. Phase 12

Not applicable: no nontrivial cover was found for distinct masses.

## REFEREE AUDIT OF THE THREE PREVIOUS PACKAGES

**Blind Claude package.**
* It quotients by Galilei × **O(2)**, which includes reflection, but signed area is not
  reflection-invariant. Its framing is underspecified for the signed-area problem. The area
  results are fine if read as SO(2) results.
* "Area 6-jet: 33 preimages" and "7-jet unique in all tests (1000 starts)" come from multistart
  searches. The word "exact preimages" is **overstated**. The true complex sixth-jet degree is
  ≥ 24 148 (certified here). Multistart found a tiny real subset and is no evidence about global
  injectivity.
* Its conjecture C2 (J⁶ many-to-one, J⁷ generically injective) is consistent with the evidence
  here for the area.

**Blind GPT second-stage package.**
* The universal generic sixth-order area rank certificate, Δ_A = s⁹P₄ + O(s⁷) with P₄ > 0, was
  checked algebraically. It is correct, and agrees with our exact nonzero determinants.
* The equal-mass label ambiguities (6-fold for all masses equal, 2-fold for one equal pair) are
  correct.
* Its final statement keeps distinct-mass generic uniqueness as a CONJECTURE, which is
  appropriate.
* Several "exact" labels refer to rational evaluations at specific witnesses. They are fine as
  witnesses, but they are local-rank results, not fibre results.

**Third/fourth-stage package.**
* Correct as audited: the differential-field closure, the singular-point location of false
  seventh-jet matches, the toy d = 2 example, the affine/point symmetry classification (its
  algebra was re-checked), and the equilateral rank-11 data (independently reproduced
  numerically at one mass triple).
* The degree bound "≥ 37 in A⁽⁷⁾" is correct but very weak; it is now ≥ 24 148.
* The statement that N is "the generic complex/algebraic sixth-jet ambiguity" is right, but the
  package gives no indication that N is of order 10⁴. The 37 real roots are a small real slice.
* The Noetherian existence of a uniform order is correctly credited to the third stage. It is
  not new.
* The 37 real roots themselves: all 37 are re-certified here, but they are not the whole real
  fibre (≥ 54 real physical roots at that target, §9).
* None of the packages computed or bounded N from above. None separated "closure under
  monodromy loops" from completeness. This investigation shows that the distinction matters
  (the false plateau at 19 893).

## Evidence ledger

| Claim | Status |
|---|---|
| ≥ 24 148 distinct regular complex sixth-jet preimages at one target, (1,2,3) | EXACT COMPUTATIONAL CERTIFICATE |
| Their A⁽⁷⁾ values pairwise distinct | EXACT COMPUTATIONAL CERTIFICATE |
| n ≥ 24 148, N ≥ 24 148, d ≤ N/24 148, for (1,2,3) and generic masses | PROVED |
| Differential embeddings / correspondence formulation of d | PROVED |
| [C(L):C(𝒜)] ≤ d; d = 1 ⇒ H, J ∈ ℚ(A..A⁽⁷⁾) | PROVED |
| d = 1 ⇐ completeness of one certified fibre plus properness at c₀ | PROVED |
| ≥ 54 real physical regular sixth-jet preimages (incl. all 37 of the third stage) at the third-stage target, A⁽⁷⁾-separated | EXACT COMPUTATIONAL CERTIFICATE |
| Equal masses (1,2,2): pipeline detects the swap+reflection cover (10 658 coincident pairs) | STRONGLY SUPPORTED NUMERICALLY (positive control) |
| Fibre is ≈ 98.8 % complete; N ≈ 24.4·10³ | STRONGLY SUPPORTED NUMERICALLY |
| d = 1 (h₇ primitive; generic complete-history uniqueness; order 7 generically sufficient) | STRONGLY SUPPORTED NUMERICALLY |
| Closure under one fixed loop set ⇒ complete fibre | FALSE (false plateau at 19 893) |
| Candidate symmetries (spin flip etc.) for unequal masses | FALSE (first fail at order ≤ 2) |
| Exact value of N, exact uniform order k₀ | UNRESOLVED |
| Generic uniqueness as a theorem | UNRESOLVED (reduced to completeness of one fibre) |
