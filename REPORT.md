# Emergent Gaussian channels for Pauli noise: new quantum-capacity thresholds

*Work log and results from a Claude Code session, 11 Oct 2026. Code in `pauli-emergence/` and
`depolarizing-threshold/`.*

## Resumo (PT)

Um canal quântico "ainda transmite informação quântica" até um limiar de ruído p_th que, para os canais
de Pauli mais estudados, continua em aberto. Em 1º/out/2026, Ahmed, Bhalerao, Lee, Leditzky, Leung,
Schaeffer e Smith (arXiv:2609.39747) mostraram que o canal despolarizante, agindo sobre estados
simétricos de muitos qubits, se torna um canal gaussiano bosônico. Este trabalho:

1. **generaliza a emergência para qualquer canal de Pauli**, com fórmula fechada (X, Y), e a valida por
   cálculo exato com M qubits (convergência 1/M);
2. **mostra o teto da técnica**: o canal emergente é antidegradável sse λ⊥ ≤ 1/√2, exatamente o
   limite da clonagem *covariante de fase*; para o despolarizante isso dá p ≤ 0,2197, e para BB84 e
   2-Pauli dá 0,1464;
3. **bate os recordes de maio/2026** (Agarwal et al., arXiv:2605.09138) para o canal X–Z independente
   (BB84) e o canal 2-Pauli, com estados explícitos de M = 60 qubits e **certificado rigoroso**
   (aritmética de bolas, sem usar nenhum teorema de limite):
   BB84 0,118371 → **0,126** (M = 240); 2-Pauli 0,118067 → **0,1255** (M = 240), com códigos
   de resíduo mod 3. Os pontos de
   M = 60 (0,1217188 e 0,1210938) foram certificados por dois métodos independentes e revisados por um
   agente separado, que tentou quebrá-los e não conseguiu (ver `CERTIFICATES.md`);
4. para o canal despolarizante, dá o melhor valor certificado com M explícito: 3p = **0,1977936** com
   M = 240 (o melhor certificado anterior era 0,194868 com n = 45; o recorde 0,2029 é assintótico, sem M
   explícito); no limite gaussiano, com aritmética de 128 bits, certificamos **3p = 0,19905**, que,
   com o teorema de convergência de Ahmed et al., é o maior limite com sinal final certificado que
   conhecemos (o 0,2029 deles vem de um valor em ponto flutuante, I_c ≈ 10⁻¹²; §7);
5. no limite gaussiano, códigos de resíduo (tipo gato de 4 componentes) dão positividade **certificada**
   do canal emergente do 2-Pauli em p = 0,124 (resultado condicional a um teorema de convergência; o
   certificado finito de M = 180, sem teorema algum, já vai além: 0,1251875).

6. **canal de ruído gaussiano aditivo** (deslocamento aleatório gaussiano, o modelo básico de ruído de luz e
   micro-ondas): o limite do BB84 é exatamente esse canal. Certificamos **Q > 0 para σ² ≤ 0,4075
   (σ ≤ 0,6384)**, além de σ = 1/√e ≈ 0,6065, onde a informação coerente gaussiana (Holevo–Werner, 2001)
   se anula e onde os códigos GKP param (≈ 0,598; PRA 111, 052445 (2025), que trata a faixa
   1/√e ≤ σ < 1/√2 como aberta). Conferido por uma implementação independente (§6.4).

7. **atenuador térmico com 1 fóton térmico** (o caso estudado por Mele, Catalano, Fanizza, Giovannetti e
   Lami em jul/2026, que certificaram Q > 0 até η = 0,7841 e estimaram numericamente 0,7818): certificamos
   **Q > 0 até η = 0,78** (§6.6), conferido de forma independente.

O que **não** funcionou está em §7: o recorde assintótico do despolarizante (0,2029) não foi superado; em
alta precisão, os códigos de envelope gaussiano dão I_c negativo nesse ponto (§7).

## 1. Problem

For a qubit Pauli channel N(ρ) = p_I ρ + p_X XρX + p_Y YρY + p_Z ZρZ, the quantum capacity is
Q = lim (1/n) max I_c(ρ, N^{⊗n}). Positivity of I_c(ρ, N^{⊗M}) for any M-qubit input certifies Q > 0.
The threshold p_th = sup{p : Q(N_p) > 0} is unknown for the depolarizing channel (1−3p, p, p, p),
the independent X–Z ("BB84") channel ((1−q)², q(1−q), q², q(1−q)) and the 2-Pauli channel
(1−2p, p, 0, p). Lower bounds come from degenerate codes; random (hashing) codes give 0.0631 (per-Pauli
depolarizing), 0.1100 (BB84) and 0.1135 (2-Pauli).

State of the art before this work (all 2026):

| channel | best lower bound on p_th | source |
|---|---|---|
| depolarizing (total error 3p) | 0.2029161 | Ahmed et al., arXiv:2609.39747 (Gaussian limit) |
| independent X–Z / BB84 | 0.118371 | Agarwal et al., arXiv:2605.09138 (n = 30) |
| 2-Pauli | 0.118067 | Agarwal et al., arXiv:2605.09138 (n = 30) |

## 2. Results at a glance

| channel | hashing | previous record | **this work, certified, M = 60** | **certified, M = 120** | Gaussian limit (this work) | ceiling of the route |
|---|---|---|---|---|---|---|
| BB84 (q) | 0.110028 | 0.118371 | 0.1217188 | 0.1249688 (mod 3); **0.126 at M = 240** | ≥ 0.1205 (Gaussian inputs); ≥ 0.122 (mod-4, numerical) | 0.146447 |
| 2-Pauli (p) | 0.113546 | 0.118067 | 0.1210938 | 0.1239688 (mod 3); 0.1251875 at M = 180; **0.1255 at M = 240** | ≥ 0.124 (mod-4 code, I_c certified) | 0.146447 |
| AGN channel (σ², Gaussian limit of BB84) | 1/e = 0.367879 (Gaussian inputs) | 1/e | — | — | **Q > 0 for σ² ≤ 0.4075, certified** (§6.4) | 1/2 |
| depolarizing (3p) | 0.18929 | 0.2029161 (asymptotic, M not explicit); explicit certified 0.194868 (n = 45) | 0.1956561 | 0.1970625; **0.1977936 at M = 240** (best explicit finite-size value) | **0.19905, I_c certified** (§7); 0.2029 (Ahmed et al., floating point) | 0.219670 |

Every "certified" entry is an explicit M-qubit state with a rigorous interval for I_c (list, hashes and
commands in [`CERTIFICATES.md`](CERTIFICATES.md)); none uses a limit theorem.

## 3. Emergent Gaussian channel of an arbitrary Pauli channel

Let λ_a = 1 − 2(p_b + p_c) be the Bloch shrinking factors and take the pole along axis c. For inputs in
the symmetric subspace with a bounded number k of excitations below the pole (Holstein–Primakoff mode b,
J_− ≈ √(2j) b), the output of N^{⊗M} splits into spin sectors J (Schur–Weyl), concentrated at
J ≈ λ_c M/2. Renormalising the mode in each sector:

* first moments: the collective spin transforms as J_a → λ_a J_a, so in sector J the quadratures scale as
  x → (λ_a/√λ_c) x, p → (λ_b/√λ_c) p;
* second moments: the vacuum (all spins at the pole) maps to a product of qubits with Bloch vector
  λ_c ẑ, whose sector components have geometric excitation number with mean (1−λ_c)/(2λ_c), i.e.
  quadrature variance 1/(2λ_c).

Hence the emergent single-mode Gaussian channel (vacuum variance ½)

    X = diag(λ_a, λ_b) / √λ_c ,        Y = diag(1 − λ_a², 1 − λ_b²) / (2 λ_c).

* For the depolarizing channel (λ_a = λ_b = λ_c = η) this is the thermal attenuator of Ahmed et al.
  (transmissivity η, environment 1/(2η)).
* Its complete-positivity condition Y + (i/2)(Ω − XΩXᵀ) ≥ 0 reduces to
  1 − λ_a² − λ_b² − λ_c² + 2λ_aλ_bλ_c ≥ 0, the Gram determinant, which holds for every Pauli channel.
* With the pole on the axis of phase covariance (λ_a = λ_b = λ⊥, λ_c = λ∥) the channel is
  phase-insensitive with τ = λ⊥²/λ∥ and Y = (1−λ⊥²)/(2λ∥). For BB84 with the pole on y this is the
  additive Gaussian noise channel σ² = (1−s²)/(2s²), s = 1−2q; for 2-Pauli it is an amplifier
  τ = (1−2p)²/(1−4p).

**Mechanism (route to a proof).** With the pole on the phase-covariance axis, a unital phase-covariant
Pauli channel is the convex mixture

    N = λ⊥ id + (λ⊥ − λ∥) F + (1 − 2λ⊥ + λ∥) R        (BB84: s id + s(1−s) F + (1−s)² R;  2-Pauli: (1−2p) id + 2p F)

whenever the coefficients are nonnegative, where R(ρ) = Tr ρ · 𝟙/2 (the replacement used by Ahmed et al.) and
F(ρ) = (XρX + YρY)/2 is "measure σ_z, then flip". On a symmetric state with few excitations, F applied to a
random fraction f of the qubits (i) removes the excitations it finds, i.e. pure loss T = 1 − f, and (ii)
appends a block of qubits polarised towards the opposite pole; recoupling to the total spin, the
Holstein–Primakoff mode of that block enters phase-conjugated, i.e. a quantum-limited amplifier
G = (1 − f)/(1 − 2f). For 2-Pauli, f = 2p gives A_G ∘ L_T with GT = (1−2p)²/(1−4p) = τ and added noise
G(1−T)/2 + (G−1)/2 = 2p(1−p)/(1−4p) = Y, exactly the formula above. This is the analogue of the
"loss + universal cloning" decomposition (Theorem 2) on which the depolarizing proof of Ahmed et al. rests,
with the cloner replaced by coupling to a flipped block.

**Exact validation.** `pauli-emergence/bandsym.py` evolves permutation-invariant operators exactly
(band-restricted Dicke-basis generator, float64; see §5). Thermal input on the Dicke ladder, BB84:

| q | N̄ | M = 8 | 16 | 32 | 64 | 128 | 256 | Gaussian prediction |
|---|---|---|---|---|---|---|---|---|
| 0.10 | 1 | 0.0922 | 0.1203 | 0.1427 | 0.1546 | | | 0.1668 |
| 0.1184 | 4 | | | | −0.0711 | −0.0314 | −0.0117 | +0.0077 |

The gap to the prediction is c/M (c ≈ 5.0 for N̄ = 4, ≈ 1.25 N̄ in general); Richardson extrapolation of
the q = 0.1184 column gives +0.0080.

## 4. Ceiling of the Gaussian route: phase-covariant cloning

Lami, Khatri, Adesso, Wilde (PRL 123, 050501, 2019) give the antidegradability of thermal attenuators
(τ ≤ (N+½)/(N+1)) and amplifiers (G ≥ (N+½)/N). Both are the single condition **Y ≥ τ/2** for a
phase-insensitive channel (τ, Y). For the emergent channel this is

    (1 − λ⊥²)/(2λ∥) ≥ λ⊥²/(2λ∥)   ⟺   λ⊥ ≤ 1/√2,

independent of λ∥. 1/√2 is the shrinking factor of the optimal phase-covariant 1→2 cloner. So codes that
live near one pole of the symmetric subspace (bounded excitations, M → ∞) can never certify Q > 0 beyond
the phase-covariant cloning bound: p ≤ (1 − 1/√2)/4 per Pauli = **0.21967 total** for the depolarizing
channel (whose true antidegradability point, from universal cloning, is 0.25), and q, p ≤ (1 − 1/√2)/2 =
0.146447 for BB84 and 2-Pauli. For BB84 this is the familiar 14.6% error rate of the optimal
phase-covariant cloning attack.

**Where the route can help.** For a phase-covariant channel (p_a = p_b) with the pole on the third axis c,
the ceiling λ⊥ = 1 − 2(p_b + p_c) ≤ 1/√2 lies above the hashing threshold only when p_c is not too large.
For Z-biased noise (p_X = p_Y = p/(2(b+1)), p_Z = pb/(b+1)) with the pole on z:

| bias b | hashing | Gaussian route (thermal inputs) | ceiling |
|---|---|---|---|
| 0 (≡ 2-Pauli, total error) | 0.22709 | 0.24149 | 0.29289 |
| 0.25 | 0.19260 | 0.20214 | 0.24408 |
| 0.5 (depolarizing) | 0.18929 | 0.18470 | 0.21967 |
| 1 | 0.19378 | 0.16856 | 0.19526 |
| 3 | 0.22221 | 0.15234 | 0.16737 |
| 10 | 0.27791 | 0.14643 | 0.15342 |

So the method is useless for strongly dephasing-biased channels (b ≳ 1, ceiling below hashing) and
useful between the 2-Pauli and depolarizing channels, i.e. two equal error rates larger than the third
(§6.2).

## 5. Exact finite-M computation

* The y-frame BB84 and 2-Pauli channels are unital and phase covariant about the pole, so
  N = exp(L) with L = g D[σ−] + g D[σ+] + (γ_D/4) D[σ_z], g = −ln(λ∥)/2, γ_D = −2 ln λ⊥ + ln λ∥
  (γ_D < 0 for 2-Pauli; the composite map is still the CPTP channel N).
* On permutation-invariant operators the generator is the PIQS generator (Shammah et al., PRA 98,
  063815). `qutip.piqs` evaluates its coefficients in **float32** (relative error ~10⁻⁸); `dicke_lindblad.py`
  and `bandsym.py` re-implement them in float64.
* L conserves δ = m − m′, so each band is an invariant subspace of size ~M²/4 (instead of ~M³/6).
* Validation against brute-force 2^M simulation, random rank-2 inputs, M = 3…6, BB84 / depolarizing /
  2-Pauli: |ΔI_c| ≤ 1.5×10⁻¹⁵.
* `finite_opt.py`: the channel is linear, so outputs of all |k⟩⟨l| are computed once per (M, p), and
  rank-r inputs are optimised with exact gradients (L-BFGS). Reproduces the n = 30 records of Agarwal et
  al. (I_c > 0 at 0.1185 for BB84 and 0.118 for 2-Pauli) from a parity-code start in ~10 s.

## 6. New thresholds

### 6.1 Exact, M = 60 (rank-2 parity codes around the y pole)

Bisection with warm starts (`threshold_M.py`, I_c > 10⁻⁶ required):

| channel | last positive p | I_c (bits, M qubits) | first failure |
|---|---|---|---|
| 2-Pauli | 0.1210938 | 4.41×10⁻⁶ | 0.121125 |
| BB84 | 0.1217188 | 6.34×10⁻⁶ | 0.12175 |
| depolarizing (per Pauli; mod-3 codes, z pole) | 0.0652187 (3p = 0.1956561) | 1.55×10⁻⁶ | 0.06525 |
| depolarizing, M = 90, K = 80 | 0.0654687 | 3.08×10⁻⁶ | 0.0655 |

Independent cross-check of the first two rows with `qutip.piqs` (different implementation of the
generator, float32 coefficients): I_c = 4.30×10⁻⁶ and 6.24×10⁻⁶ (difference 1×10⁻⁷, explained by the
float32 coefficients).

M = 120, K = 90 (bisection in float64, then certified with `certify_big.py`, see below):

| channel | last positive p | I_c | first failure | `qutip.piqs` cross-check |
|---|---|---|---|---|
| 2-Pauli | 0.1232188 | 2.10×10⁻⁶ | 0.12325 | 2.02×10⁻⁶ |
| BB84 | 0.1239062 | 1.45×10⁻⁶ | 0.1239375 | 1.36×10⁻⁶ |
| depolarizing (per Pauli) | 0.0656875 (3p = 0.1970625) | 1.27×10⁻⁶ | 0.0657187 | — |

For the depolarizing channel the explicit-M thresholds saturate (0.06522 → 0.06547 → 0.06569 per Pauli
for M = 60, 90, 120), well below the asymptotic 0.06764 (= 0.2029/3): the record codes of Ahmed et al.
need ~56 excitations, out of reach of explicit M here.

The sequence M = 30 → 60 → 120 (BB84 0.1185 → 0.1217 → 0.1239; 2-Pauli 0.1180 → 0.1211 → 0.1232)
still rises with roughly geometric increments (ratio ≈ 0.7), which would extrapolate to ≈ 0.128–0.129,
above the Gaussian-limit codes found so far. The best limit codes are therefore not yet known.

States: `pauli-emergence/runs/cert_2p_M60_p0.1210938.npy`, `runs/cert_bb_M60_thr.npy`,
`runs/cert_dep_M60_thr.npy` (rows = the two unnormalised codewords in the excitation basis of the Dicke
ladder of the pole: y for BB84/2-Pauli, z for depolarizing).

**Certificates.** Two independent rigorous certifiers (details and all intervals, rounded outward, in
[`CERTIFICATES.md`](CERTIFICATES.md)):

* `certify_finite.py`: generator entries, rates and square roots in Arb balls, `arb_mat.exp` per band,
  Weyl/Ostrowski eigenvalue enclosures. Cost ~n³ per band (≈ 40 min at M = 60).
* `certify_big.py`: Arb generator, float64 scaling-and-Taylor integration with a rigorous 1-norm error
  budget (the band generators of BB84 and depolarizing have logarithmic 1-norm μ₁ ≤ 0, i.e. the evolution
  is a contraction in the coordinates; for 2-Pauli μ₁ ≈ 6.8 at M = 120 and the budget includes e^{μ₁}),
  Duhamel bound for the float rounding of the generator, Arb entropy of the computed blocks plus the
  continuity term ε′ ln(D/ε′). Validated on the three M = 60 states against `certify_finite.py`
  (its intervals contain the Arb-exp ones); ≈ 10 min at M = 120.
  For 2-Pauli the dephasing rate of the generator is negative and μ₁ grows linearly with M (≈ 3.2, 6.8,
  10, 14 at M = 60, 120, 180, 240), so float64 rounding errors can be amplified by up to e^{μ₁}. With
  `CERT_LD=1` the band integration runs in x87 extended precision (u = 2⁻⁶⁴): at M = 60 its interval
  [4.4142155, 4.4142560]×10⁻⁶ is 1000× tighter than in float64 and contains the Arb-exp interval; it is
  ~7× slower. A self-audit found that the first version of this mode formed the Taylor step scalar h/n
  as a float64 quotient (error 2⁻⁵⁴, not covered by the 2⁻⁶⁴ budget); it is now formed in extended
  precision, the M = 60 check was repeated with the same interval, and the two M = 180 runs started
  with the faulty version were discarded (`runs/*_unsound.log`). Float64 mode, used for every other
  certificate, was not affected (its 4u term covers that rounding).

| channel | p (exact) | M | I_c (bits, outward) |
|---|---|---|---|
| 2-Pauli | 605469/5000000 = 0.1210938 | 60 | [4.41423406, 4.41423749]×10⁻⁶ |
| BB84 | 304297/2500000 = 0.1217188 | 60 | [6.34444364, 6.34445247]×10⁻⁶ |
| depolarizing | 652187/10⁷ = 0.0652187 (3p = 0.1956561) | 60 | [1.54889076, 1.54890269]×10⁻⁶ |
| BB84 | 619531/5000000 = 0.1239062 | 120 | [1.44239924, 1.44929662]×10⁻⁶ |
| 2-Pauli | 308047/2500000 = 0.1232188 | 120 | [1.73815583, 2.46289307]×10⁻⁶ |
| depolarizing | 1051/16000 = 0.0656875 (3p = 0.1970625) | 120 | [1.26617376, 1.27106873]×10⁻⁶ |

Mod-3 residue codes (codewords on k ≡ 0 and 1 mod 3 around the pole) beat parity codes for BB84 and
2-Pauli as well:

| channel | p (exact) | M | code | I_c (bits, outward) |
|---|---|---|---|---|
| BB84 | 156211/1250000 = 0.1249688 | 120 | mod 3 | [3.32197192, 3.32992534]×10⁻⁶ |
| 2-Pauli | 154961/1250000 = 0.1239688 | 120 | mod 3 | [4.05046902, 4.96028384]×10⁻⁶ |
| 2-Pauli | 2003/16000 = 0.1251875 | 180 | mod 3 (extended precision) | [1.17881571, 1.19653069]×10⁻⁶ |
| 2-Pauli | 251/2000 = 0.1255 | 240 | mod 3 (extended precision) | [1.88681426, 2.70734623]×10⁻⁶ |
| BB84 | 63/500 = 0.126 | 240 | mod 3 | [5.54571350, 5.56189889]×10⁻⁶ |

Hence Q > 0 at these points: p_th(BB84) ≥ 0.126, p_th(2-Pauli) ≥ 0.1255 and, for the depolarizing
channel, 3p_th ≥ 0.1977936 with an explicit M = 240 state (previous best certified explicit value:
0.194868 at n = 45, Krohn-Grimberghe, arXiv:2608.15870). BB84 and depolarizing are monotone in p
(composition), so Q > 0 on the whole interval below; for 2-Pauli the bound is in the sup sense.

**Independent review.** A separate agent tried to break the M = 60 certificates (frame change, generator
vs qutip.piqs line by line, an exact Schur–Weyl check of the band evolution at M = 60, its own brute force
for M ≤ 10 against both certifiers, the residue pattern of the states, a fuzz test of the entropy
enclosure). Verdict: all valid; the minor issues it found (outward rounding of quoted intervals, a
fail-safe NaN case in the enclosure, provenance of logs) are fixed. A second agent then reviewed
`certify_big.py`, on which every M ≥ 120 entry rests (error budget line by line, a new brute force for
M ≤ 10, true integration errors against Arb references, at most 1% of the bound). Verdict: sound; its
robustness suggestions (platform check for extended precision, explicit errors, input validation) were
applied, and all certificates are being re-run with the final code (`CERTIFICATES.md`).

### 6.2 A family between 2-Pauli and depolarizing

(p_X, p_Y, p_Z) = p (1, 1, r)/(2 + r), total error p, pole on z. r = 0 is the 2-Pauli channel (p = 2 p_2P),
r = 1 the depolarizing channel (p = 3 p_dep). Certified points (M = 60, `certify_big.py`; best of parity and
mod-3 codes found so far):

| r | hashing | certified threshold (M = 60) | code | I_c (bits) |
|---|---|---|---|---|
| 0 | 0.227092 | 0.2421876 parity; 0.2438124 mod 3 (M = 120: 0.2479376) | | see §6.1 |
| 0.25 | 0.199895 | 0.2228375 = 17827/80000 | mod 3 | [1.32667, 1.33010]×10⁻⁶ |
| 0.5 | 0.192599 | 0.2099375 = 3359/16000 | mod 3 | [3.73849, 3.74093]×10⁻⁶ |
| 0.75 | 0.189926 | 0.2015875 = 16127/80000 | mod 3 | [1.63392, 1.63633]×10⁻⁶ |
| 1 | 0.189290 | 0.1956561 (M = 120: 0.1970625) | mod 3 | see §6.1 |

The gain over hashing is 0.006–0.023 across the family and the certified thresholds decrease
monotonically from 2-Pauli to depolarizing. (A first mod-3 run at r = 0.5 lost the state at one step and
gave only 0.1987; parity gave 0.201725 = 8069/40000, also certified; the retry gave 0.2099375.) No point-by-point records for the
interior of this family were found in the literature (Bausch–Leditzky 2021 cover the Pauli simplex with
graph-state codes; we did not compare point by point).

**Mod-3 codes beat parity also at the endpoints.** Re-running M = 60 with mod-3 codes: 2-Pauli
0.1219062 (parity 0.1210938), BB84 0.1227813 (parity 0.1217188). Larger-M mod-3 runs are in
`runs/t*3_M120.log`.

### 6.3 Gaussian limit

* Gaussian (thermal) inputs, infinite energy: BB84 0.120488, 2-Pauli 0.120743, depolarizing 0.18470 (below
  hashing, as Ahmed et al. found; their mod-3 codes give 0.2029).
* Residue codes (codewords on n ≡ a_i mod k): 2-Pauli, mod-4 {0,2}, p = 31/250 = 0.124:
  **I_c ∈ [1.96141999×10⁻⁴, 1.96142031×10⁻⁴] bits, certified** (`certify_limit.py`; output truncated at
  D = 280 with a flag state, flag mass 2×10⁻¹⁷; T = 94/125, G = 94/63 exact).
  Turning this into p_th(2-Pauli) ≥ 0.124 would need the convergence theorem of Ahmed et al. extended from
  the depolarizing channel to this channel (§8); the explicit M = 180 certificate (0.1251875) makes that
  step unnecessary for the record.

### 6.4 The additive Gaussian noise channel: Q > 0 beyond σ = 1/√e

For the BB84 (independent X–Z) channel the emergent channel of §3 has X = 𝟙 (τ = λ⊥²/λ∥ = 1): it is the
additive Gaussian noise channel N_σ (classical noise / Gaussian random-displacement channel), which adds
variance σ² = Y = (1 − (1−2q)²)/(2(1−2q)²) to each quadrature (q = (a + a†)/√2, vacuum variance ½). Its
quantum capacity is unknown. Gaussian inputs give Q ≥ max(0, log₂(1/(eσ²))) (Holevo–Werner), positive
only for σ < 1/√e ≈ 0.6065; GKP codes with maximum-likelihood decoding reach σ ≈ 0.598, and the window
1/√e ≤ σ < 1/√2 (σ² ≥ ½ is antidegradable, which is also the ceiling of §4) is stated as open (Lin and
Noh, arXiv:2411.04277, PRA 111, 052445 (2025)). Since N_σ₂ ∘ N_σ₁ = N_√(σ₁²+σ₂²), Q is monotone in σ.

Result: the mod-3 residue code obtained at M = 240 for BB84 and re-optimised in the Gaussian limit
(codewords on n ≡ 0 and n ≡ 1 mod 3, mean photon number ≈ 15, support n ≤ 140) gives, at σ² = 2/5:

    I_c ∈ [3.35715244, 3.35738269] × 10⁻⁶ bits      (certify_limit.py agn 2/5, D = 320; also D = 380)

hence Q(N_σ) > 0 for every σ² ≤ 2/5. Re-optimising the same family with a larger cutoff (N_in = 300,
D = 600; `agn_scan3.sh`) and certifying at σ² = 81/200:

    I_c ∈ [6.33597116, 6.40767279] × 10⁻⁸ bits      (certify_limit.py agn 81/200, N_in = 300, D = 600)

hence **Q(N_σ) > 0 for every σ² ≤ 81/200 = 0.405, i.e. σ ≤ 0.6364**. Rigour: the channel is A_G ∘ L_T with
G = 7/5, T = 5/7 exact; Fock matrix elements in Arb; output levels above D collapsed into a flag (a channel
on Bob's side, so the certified number is a lower bound by data processing; flag mass ≈ 4×10⁻⁴⁷);
entropies by Weyl/Ostrowski enclosures. Other residue families with the same cutoff (float thresholds, I_c > 10⁻⁶ required; `runs/agnL_*.log`):
mod 4 {0,2} (≈ 72 photons) 0.40175, parity 0.39975, mod 5 {0,1} 0.39975, mod 5 {0,2} 0.39575,
mod 4 {0,1} 0.39375, mod 6 {0,3} 0.39175; mod 3 {0,1} 0.4035 is the best residue family. Rank-3 mod-3
inputs do not help (the optimum stays rank 2).

Breaking the symmetry helps. Starting from the mod-3 state plus noise and optimising over all rank-2
inputs (N_in = 200, D = 400, `agn_general.py`, `agn_general_thr.py`), the optimum stays ≈ 99.9% on the
residues but its small off-residue part (≈ 0.1% per codeword, overlap ⟨c₀|c₁⟩ ≈ −7×10⁻⁴) multiplies the
coherent information: 4.37×10⁻⁷ instead of 6.4×10⁻⁸ at σ² = 0.405. Certified with a single R⊗B block
(k = 1; for such inputs the R⊗flag block is not diagonal and its entropy is enclosed by
[h(t), h(t) + t ln 2], t = flag mass):

    σ² = 81/200:   I_c ∈ [4.36362113, 4.37080935] × 10⁻⁷ bits
    σ² = 407/1000: I_c ∈ [5.43799903, 5.50010769] × 10⁻⁸ bits
    σ² = 163/400:  I_c ∈ [5.97733793, 6.88032591] × 10⁻⁹ bits     (N_in = 300, D = 600; negative at 0.408)

hence **Q(N_σ) > 0 for every σ² ≤ 0.4075, i.e. σ ≤ 0.6384**. A second continuation from the high-energy
mod-4 {0,2} state (≈ 77 photons) with its symmetry broken stops at 0.406 (`runs/agnB.log`), so for
single-mode rank-2 inputs the frontier found is σ² ≈ 0.4075–0.408. The same independent agent recomputed both
symmetry-free states with full matrices (no symmetry assumed) and two different purifications:
4.3682405×10⁻⁷ and 5.4779650×10⁻⁸ bits, inside the intervals (and 6.58088×10⁻⁹ for the σ² = 163/400
state, by three routes); it audited the general flag bound
(and stress-tested it on toy inputs with flag mass up to 2×10⁻² and strongly non-diagonal flag blocks).
Its one remark, that `certify()` did not check residue support when k > 1, is fixed: such inputs are
now rejected.

What the codes look like. They are not three-component cat codes: superposing |α⟩, |αω⟩, |αω²⟩ into
the n ≡ 0, 1 (mod 3) codewords and optimising α and the weights gives I_c < 10⁻¹² for every σ² in
[0.36, 0.39], i.e. nothing. The optimal codewords have a positive main lobe (n ≈ 3–18) followed by a tail
of alternating sign, with photon-number variance ≈ 2.7 times the Poisson value (`runs/agn_k3_sigma2_0.4.npy`).

A simple physical family explains most of the gain: **three-component squeezed cats**, i.e. the n ≡ 0 and
n ≡ 1 (mod 3) parts of a displaced *phase-squeezed* state D(α)S(r)|0⟩ (three parameters α, r, w;
`agn_sqcat.py`). Plain cats (r = 0) give nothing; with squeezing r ≈ −1 and α ≈ 5.7 (⟨n⟩ = α² + sinh²r ≈ 34) the
input stays positive up to σ² = 0.392, certified at σ² = 49/125:

    I_c ∈ [1.73865751, 1.73880842] × 10⁻⁶ bits      (certify_limit.py agn 49/125, N_in = 140, D = 320)

so about 60% of the gain over 1/e already comes from a state of a kind prepared in optics labs; the fully
optimised codes (super-Poissonian, with an alternating tail, then a slight Z₃ breaking) add the rest.

At finite M the same symmetry breaking gains little: for BB84 at M = 60 the threshold moves from
0.1227813 (mod 3) to 0.1228125 (`runs/tbbgen_M60.log`, float). The large effect is specific to the limit.

Consequence. At σ² = 2/5 the coherent information optimised over Gaussian inputs is 0 while
Q ≥ 3.357×10⁻⁶, so the Gaussian one-shot coherent information is *not* the quantum capacity of this
channel. Harrington and Preskill (PRA 64, 062301 (2001), Conclusions) wrote that "it is tempting to
speculate that this optimal one-shot coherent information actually is the quantum capacity of the
channel"; the certificate settles that in the negative. Prior claims beyond 1/√e: Zhang, Wu and Guo (PRA
107, 062408 (2023)) report a code-capacity threshold σ ≈ 0.67 for rectangular GKP + XZZX surface codes
from matching-decoder simulations up to distance 9; Lin and Noh (arXiv:2411.04277, App. F) find the
crossings drift downward with distance and do not accept it as a threshold above 1/√e. That is a
simulation, not a proof; to our knowledge the certificate here is the first rigorous Q > 0 beyond 1/√e.
Literature check (12 Oct 2026, all quant-ph listings May–Oct 2026, cs.IT Sep–Oct, citations of
arXiv:2609.39747): no prior or concurrent result for any of the claims in this report.

Corollary (all phase-insensitive channels). A phase-insensitive Gaussian channel with transmissivity or
gain τ and added noise Y becomes N_σ with σ² = (Y + |1−τ|/2)/τ after a quantum-limited amplifier of gain
1/τ (τ < 1) or a pure loss of transmissivity 1/τ (τ > 1), i.e. N_σ is a degraded version of it; hence
Q > 0 whenever (Y + |1−τ|/2)/τ ≤ 0.4075. For the thermal attenuator (Y = (1−η)(N+½)) this reads
(1−η)(N+1)/η ≤ 0.4075, and for the thermal amplifier (Y = (τ−1)(N+½)) it reads (τ−1)(N+1)/τ ≤ 0.4075. Both
go beyond the Holevo–Werner (Gaussian-input) thresholds when the environment has N > 4.23 thermal
photons: e.g. at N = 10 the new windows are η ∈ [0.96428, 0.96614) and τ ∈ (1.03632, 1.03847]. (Mele et
al. treat N = 1, where this corollary gives nothing new; conversely their ν = 1, η = 0.7841 point implies
AGN positivity only at σ² = (1−η)ν/η ≈ 0.275 < 1/e, via E = L_η ∘ N_σ.)

Independent verification (separate agent, no repository code imported): the channel implemented as the
diffusion Lindbladian exp(σ²(D[a] + D[a†])), checked to add exactly σ² to both quadrature variances and to
map the vacuum to a thermal state with n̄ = σ²; matrix elements cross-checked four ways (60-digit quadrature
of the displacement average, the Lindblad propagator, A_G∘L_T and L_T′∘A_G′) to 10⁻¹⁵; I_c of the stored
state 3.3572969×10⁻⁶ bits at σ² = 2/5 (stable for cutoffs 320/400/520), inside the certified interval;
same agreement at σ² = 0.38, 0.403, 0.404 and correct failure at 0.405 for that state. For the
σ² = 81/200 state it obtained I_c = 6.38363×10⁻⁸ bits (cutoffs 600/700/800 agree to 3×10⁻¹³), inside the
certified interval, and reproduced the certificate exactly; the sign change of that state is at
σ² ≈ 0.40509. Its audit of the certifier found no
soundness issue; its one robustness remark (a flag mass below the working precision made the certifier
return NaN, failing safe) is fixed. Closest related work: Mele, Catalano, Fanizza, Giovannetti, Lami
(arXiv:2607.27449) do the analogous thing for the thermal attenuator; their certified points do not imply
this one by degradation ((1−η)ν < ησ² there).

### 6.5 A thermal amplifier beyond the Gaussian threshold

The emergent channel of the 2-Pauli family is A_G ∘ L_T with T = 1−2p, G = (1−2p)/(1−4p): a thermal
amplifier. The mod-3 code transferred from M = 240 (§6.4, `runs/lim_2p3_fromM240.log`) and continued in
the limit up to p = 0.12725 gives, at p = 51/400 (T = 149/200, G = 149/98):

    I_c ∈ [9.10836870, 9.11017611] × 10⁻⁷ bits      (certify_limit.py 2pauli 51/400, N_in = 140, D = 320)

The channel is the thermal amplifier with gain τ = 22201/19600 ≈ 1.13270 and N = 149/51 ≈ 2.9216 thermal
photons. At this gain Gaussian inputs give Q > 0 only for N < 2.6533 (log₂(τ/(τ−1)) − g(N) = −0.1187 at
N = 149/51), and the channel is antidegradable only for N ≥ 3.768 (Y ≥ τ/2), so the certificate covers
about a quarter of the open window at this gain. The AGN corollary of §6.4 does not reach this point
(its AGN-equivalent noise is σ² = 0.459 > 0.4075). No non-Gaussian positivity result for thermal amplifiers
was found in the literature check (Mele et al. treat the attenuator).

Independent check (same separate agent as §6.4, no repository code): three routes (the Lindblad
semigroup exp(κ↑D[a†] + κ↓D[a]) with κ↑ = 0.488658, κ↓ = 0.364050; L_T′ ∘ A_G′ with T′ = 12001/19600,
G′ = 22201/12001; its own A_G ∘ L_T), cross-checked against 50-digit matrix elements, all give
I_c = 9.109425×10⁻⁷ bits (±2×10⁻¹⁴ over routes and cutoffs 320–600), inside the certified interval; the
certificate reproduces exactly (and at D = 400). It also swept all single-mode Gaussian inputs (thermal
and squeezed thermal, 60-digit arithmetic): coherent information ≤ 0 everywhere, as expected.

Breaking the residue symmetry (as in §6.4; `amp_scan.sh`) extends this a little: a symmetry-free rank-2
input stays positive up to p = 0.128, certified at p = 16/125 (T = 93/125, G = 93/61; gain
τ = 8649/7625 ≈ 1.1343, N = 93/32 ≈ 2.906 thermal photons; k = 1, N_in = 200, D = 400):

    I_c ∈ [1.43395294, 1.43876228] × 10⁻⁷ bits

At this gain thermal inputs need N < 2.6206 and the channel is antidegradable for N ≥ 3.7231, so the
certificate covers about 26% of the open window there (Gaussian value −0.1273). Independently recomputed by three routes
(1.4369750×10⁻⁷ bits, inside the interval); thermal and squeezed-thermal inputs all negative there.

### 6.6 Thermal attenuator with one thermal photon: beyond Mele et al.

Mele, Catalano, Fanizza, Giovannetti and Lami (arXiv:2607.27449, Jul 2026) proved Q > 0 for the thermal
attenuator E_{η,ν} with ν = 1 down to η = 0.7841, where Gaussian inputs stop at η = 0.8 and the channel is
antidegradable for η ≤ 0.75 (their numerics, not certified, suggest positivity to η ≈ 0.7818). The same
tools apply directly (channel `th1`: A_G ∘ L_T with G = (1−η)(ν+½) + (1+η)/2, T = η/G, so X = √η and
Y = (1−η)(ν+½)). A residue scan (`th1_scan.sh`) finds mod 3 {0,1} best (≈ 17 photons; parity stops at
η = 0.78525, mod 4 {0,1} at 0.78375, mod 4 {0,2} at 0.78725), and at η = 3133/4000 = 0.78325:

    I_c ∈ [8.27408111, 8.27425105] × 10⁻⁶ bits      (certify_limit.py th1 867/4000, N_in = 80, D = 240)

Since E_{η,ν} = E_{η/η′,ν} ∘ E_{η′,ν}, Q is non-decreasing in η, so **Q(E_{η,1}) > 0 for every η ≥ 0.78325**,
improving the rigorous bound 0.7841 (at their point they certify I_c ≈ 1.4×10⁻⁷; ours at 0.78325 is
8.3×10⁻⁶). Independent check (same separate agent): the channel matches their definition (ν = environment
mean photon number; their η_AD = 0.75 and η_G = 0.8 reproduced), matrix elements checked at 50 digits
against the beam-splitter definition, four routes (thermal master equation, AGN ∘ loss, A_G ∘ L_T,
L_T′ ∘ A_G′) give I_c = 8.2742014×10⁻⁶ bits (spread 1.6×10⁻¹³), inside the interval; all thermal and
squeezed-thermal inputs give I_c < 0 there.

Further, with a larger cutoff the mod-3 family stays positive to η = 0.7815 (float, `runs/th1_k3_big.log`),
and breaking the residue symmetry (as in §6.4; `th1_scan2.sh`) keeps a rank-2 input positive to η = 0.78
(negative at 0.7795). Certified at η = 39/50 (T = 39/61, G = 61/50; k = 1, N_in = 120, D = 300):

    I_c ∈ [2.12631890, 2.16095917] × 10⁻⁸ bits

so **Q(E_{η,1}) > 0 for every η ≥ 0.78**, below both their certified 0.7841 and their uncertified numerical
0.7818. Independent check: four routes give 2.150005×10⁻⁸ bits (spread 1.7×10⁻¹³), inside the interval;
Gaussian inputs all negative (−0.174 at infinite energy). The window 0.75 < η < 0.78 remains open.

**Map over ν** (`nu_scan.sh`, `nu_scan2.sh`, `plot_attenuator_map.py` → `attenuator_map.png`). Mele et al. treat
ν = 1 only. The same codes (mod 3 {0,1}, then symmetry breaking) certify Q > 0 well inside the open window
(η_AD, η_G) for every ν tried:

| ν | η_AD (antideg.) | η_G (thermal inputs) | certified η | I_c (bits) | window closed |
|---|---|---|---|---|---|
| 1/4 | 0.6 | 0.65147 | 313/500 = 0.626 | [2.48378509, 2.48746745]×10⁻⁷ | 50% |
| 1/2 | 0.66667 | 0.72207 | 1393/2000 = 0.6965 | [6.31451396, 6.36306427]×10⁻⁸ | 46% |
| 1 | 0.75 | 0.8 | 39/50 = 0.78 | [2.12631890, 2.16095917]×10⁻⁸ | 40% |
| 2 | 0.83333 | 0.87097 | 8587/10000 = 0.8587 | [1.14709176, 1.14892934]×10⁻⁷ | 33% |
| 3 | 0.875 | 0.90459 | 1791/2000 = 0.8955 | [8.31897855, 8.33853440]×10⁻⁸ | 31% |
| 5 | 0.91667 | 0.93723 | 7459/8000 = 0.932375 | [1.17390097, 1.17400578]×10⁻⁶ | 24% |

(k = 1 symmetry-free inputs with N_in = 160, D = 360, except ν = 5: mod 3, k = 3, N_in = 260, D = 520.) Since Q is monotone in η at
fixed ν, each row gives Q > 0 for all larger η. For ν ≳ 4.2 the AGN corollary also applies (η ≥ (ν+1)/(ν+1.4075)),
but the direct certificates go further at every ν tried. Independent check (same separate agent, three
routes per point: A_G ∘ L_T, loss then AGN, thermal master equation): all five new values inside their
intervals (spread ≤ 2.2×10⁻¹³), all Gaussian inputs negative at each point; it also regenerated the
squeezed-cat state from (α, r, w) at 80 digits (max deviation 3.6×10⁻¹⁶) and found I_c = 1.7387559×10⁻⁶,
inside its interval. The ν = 1/4 and ν = 5 states pressed against N_in = 160; with N_in = 260 (`nu_scan3.sh`)
ν = 5 moved from 0.932675 to 0.932375 (certified, 47 s; the independent re-implementation gives 1.1739661×10⁻⁶, inside, routes agreeing to 1.6×10⁻¹⁴), while for ν = 1/4 the symmetry-free continuation did not
reach a positive value at η = 0.6255 (I_c = −5.8×10⁻⁸), so 0.626 stands there.

## 7. Depolarizing channel in the Gaussian limit

**Certified value.** `depolarizing-threshold/hp_cont.py` continues the mod-3 code in η with the coherent
information and its gradient computed in 128-bit ball arithmetic (`hp_ic.py`: Arb channel coefficients,
eigenbases refined to the working precision, cluster-aware Jacobi rotations; `hp_par.py` runs the six
charge blocks in parallel). Double precision has a noise floor of ~3×10⁻¹⁴ bits for this problem; 128 bits
gives values to ~10⁻²⁶ (a rigorous enclosure, `hp_ic.certify_hp`, agreed with the approximate value to
that level). From the M = 240 limit state (N_in = 240, D = 360):

| η | 3p | I_c (bits) |
|---|---|---|
| 0.73525 | 0.1985625 | +1.52×10⁻⁶ |
| 2939/4000 = 0.73475 | 0.1989375 | **[1.07683723, 1.07994494]×10⁻⁷ certified** |
| 2351/3200 = 0.7346875 | 0.198984375 | **[1.82651528, 1.85766717]×10⁻⁸ certified** |
| 0.73465625 | 0.1990078 | −1.8×10⁻⁸ (optimized) |
| 0.734625 | 0.19903125 | −1.3×10⁻⁷ (optimized) |
| 0.7345 | 0.1991250 | −4.1×10⁻⁷ (optimized) |
| 0.73425 | 0.1993125 | −9.1×10⁻⁷ (optimized) |

Certificates by `certify_limit.py depol` (the reviewed certifier L; flag masses ~10⁻⁵⁶). An independent
re-implementation written from scratch (pure loss then amplifier with long-double Kraus elements, the
thermal master equation by uniformization and by dense expm, and a beam splitter with a Fock-mixed
thermal environment; output cutoff 420, tail < 3×10⁻⁷⁵) gives 1.0787368×10⁻⁷ and 1.8455798×10⁻⁸, routes
agreeing to ≤ 7.5×10⁻¹⁴, both inside the intervals; Gaussian inputs are negative at both points
(log₂(η/(1−η)) − g(ν) = −0.1665, −0.1671). With Theorem 6 of Ahmed et al. (proved for the depolarizing
channel) this gives p_th ≥ 0.198984375 in total error for sufficiently large M. The branch has a
genuine zero crossing between η = 0.73465625 and 0.7346875 (values ~10⁻⁸, far above rounding), with the mean photon number
staying near 20: it is not the record branch, and high precision was not what limited it.

**A branch with more photons.** Stretching that state by 2 in the residue index (n = 3m + r, m → 2m;
`stretch_init.py`, mean 38.5 photons) and optimizing at η = 3673/5000 = 0.7346, where the 20-photon branch
is negative, gives a positive code: **3p = 0.19905** (p = 1327/20000 per Pauli; T = 13490929/21682500,
G = 8673/7346), certified I_c ∈ **[1.52256332, 1.54665840]×10⁻⁸** bits (`stretch2b_eta0.7346.npy`,
N_in = 240, D = 360). The independent re-implementation gives 1.5376335×10⁻⁸ (routes agree to 1.5×10⁻¹³;
output tail beyond 360 photons 1.4×10⁻⁵⁷), and Gaussian inputs are negative there (−0.1678). With
Theorem 6 of Ahmed et al.: **p_th ≥ 0.19905**. More photons, lower threshold, smaller I_c: the same
trade-off that leads to the 56-photon record codes. The continuation of this branch is running
(`runs/hpcont38.log`).

**The record point.** In 128-bit arithmetic (`pauli-emergence/dep_hp_envelope.py`, `dep_hp_scan.py`), mod-3
codes with Gaussian envelopes are negative at η = 0.7294452 for every mean (40, 56.5, 70), width (2–6)
and weight (0.30–0.62) scanned, with |I_c| from 10⁻⁵ down to 10⁻²⁸ (the optimized plateau state of
earlier work is enclosed rigorously as −4.66625545720×10⁻¹² bits); they stay negative at η = 0.735,
0.74 and 0.75. Squeezed-cat codes D(α)S(r)|0⟩ (α = 5–8, |r| ≤ 0.4) are negative too, down to −10⁻²⁸ at
64 photons. High-precision L-BFGS from the Gaussian plateau state went from −4.67×10⁻¹² to −3.69×10⁻¹²
in 20 evaluations (stopped). The record code (mean ≈ 56 photons, I_c ≈ 10⁻¹²) therefore has structure
beyond a Gaussian envelope that we did not find; without the authors' state we can neither confirm nor
refute their floating-point value.

**Earlier attempts (double precision).**

* Gaussian-limit continuation of the M = 240 mod-3 code (`dep_limit_scan.sh`, N_in = 300, D = 600): the
  limit threshold of that family is p = 0.0661875 per Pauli (3p = 0.19856), below the record 0.0676387
  (superseded by the certified 0.19905 above).
  Breaking the residue symmetry, which multiplied I_c by 7 for the AGN channel (§6.4), did not help here:
  from a perturbed start the optimiser did not even recover a positive value at p = 0.066
  (`runs/limdep_gen.log`). The record regime (≈ 56 photons, I_c ≈ 10⁻¹²) stays out of reach.
* At the record point itself (η = 0.7294452) the best mod-3 state found here has I_c = −4.7×10⁻¹² (float);
  symmetry-free rank-2 optimisation from it (two seeds, `dep_record_scan.sh`) ends at −4×10⁻¹¹ and
  −9×10⁻¹¹, i.e. slightly worse. So neither a rigorous certificate of the 0.2029 point nor an improvement
  on it was obtained; that claim of Ahmed et al. rests on their floating-point value.

* Full L-BFGS over residue codes (k = 2…6, ranks 2–3) in float64 reaches I_c > 0 only down to η ≈ 0.737
  (p ≈ 0.197); the record point η = 0.7294452 needs the regime the authors describe (mean photon ≈ 56).
* That regime is a near-classical plateau: the three phase blobs are almost orthogonal for both Bob and
  the environment, I_c ≈ 0 with curvature ≈ −0.16 and gradient ~10⁻⁶; only ~10⁻¹¹ is available there,
  which matches the 10⁻¹² margin reported (at the edge of double precision).
* Large-amplitude coherent-state cat codes: the quantum-Chernoff overlap exponents of Bob and of the
  environment cross at η = 0.7757, so this limit is worse than the record.

## 8. Status of rigour

* §3 formula: derived from exact first/second moments plus the Gaussianity of the limit; the latter is
  proven for the depolarizing channel by Ahmed et al. and here supported by exact numerics (1/M
  convergence) for BB84 and 2-Pauli. A proof for general Pauli channels is open.
* §4 ceiling: rigorous for the emergent channel (Lami et al. + algebra).
* §6.1: exact finite-M statements, all certified (M = 60 by two independent methods, M = 120–240 by
  `certify_big.py`), using no limit theorem. Remaining trust assumptions: correctness of the certifier
  code (two independent reviews) and of python-flint/Arb; IEEE round-to-nearest and scipy's CSR
  product semantics; x87 extended precision for the two 2-Pauli entries computed with `CERT_LD=1`
  (checked at run time); the PIQS generator formulas (Shammah et al.),
  checked here against brute force (M ≤ 10) and, by the reviewer, against an exact Schur–Weyl formula at
  M = 60.
* §6.3: I_c of the emergent channel is certified; the step to the Pauli channel uses the open theorem.

## 9. Open problems

* Prove the convergence theorem (Ahmed et al., Thm. 6) for phase-covariant Pauli channels (no longer
  needed for the 2-Pauli record, which the explicit M = 180 certificate already exceeds).
* Find the optimal Gaussian-limit codes: the finite-M trend points to ≈ 0.128 for BB84 and 2-Pauli, and the
  ceiling is 0.146447.
* Larger M: thresholds keep rising (M = 240 runs in `runs/t*_M240.log`); `certify_big.py` scales as the
  float computation does.
* Depolarizing: an explicit-M certificate approaching 0.2029 needs states with ~56 excitations, i.e.
  M in the thousands, or a quantitative version of the convergence theorem.

## 10. Reproduce

    pip install numpy scipy python-flint qutip
    cd pauli-emergence
    python3 piqs_ic.py                       # brute force vs permutation-invariant evolution
    python3 thermal_test.py --q 0.10 --Nbar 1 --K 20 --Ms 8,16,32,64
    python3 threshold_M.py --channel bb84 --M 60 --K 60 --p0 0.117 --tag tbb_M60
    python3 certify_finite.py bb84 304297/2500000 runs/cert_bb_M60_thr.npy 60 2 0,1
    python3 certify_big.py bb84 619531/5000000 runs/state_bb_M120_thr.npy 120 2 0,1   # all in CERTIFICATES.md
    python3 certify_limit.py 2pauli 31/250 runs/2p_k4_p0.124.npy 0,2 4 160 280
