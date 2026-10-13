# Certificate manifest

Each line certifies I_c(ρ, N^{⊗M}) > 0 for an explicit permutation-symmetric rank-2 state ρ on M qubits,
hence Q(N) ≥ I_c/M > 0 and p_th ≥ p (p_th := sup{p : Q(N_p) > 0}). Intervals are in bits and rounded
outward. Paths are relative to `physics/pauli-emergence/`.

Channels (Pauli probabilities): BB84 / independent X–Z ((1−q)², q(1−q), q², q(1−q)); 2-Pauli
(1−2p, p, 0, p); depolarizing (1−3p, p, p, p); family r: (1−p, p/(2+r), p/(2+r), pr/(2+r)) with the
pole on z (`pr<r>` in the scripts). Codes: rows of the `.npy` file are the two unnormalised
codewords in the excitation basis of the Dicke ladder around the pole (y axis for BB84 / 2-Pauli, z for
depolarizing); residues mod 2 (parity) or mod 3 as listed, with exact zeros elsewhere.

## Two independent certifiers

* `certify_finite.py` — every band of the generator exponentiated with `arb_mat.exp` (Arb ball
  arithmetic); cost ~n³ per band, practical up to M = 60.
* `certify_big.py` — generator entries in Arb; float64 scaling-and-Taylor integration with a rigorous
  1-norm error budget (truncation, rounding of every product and sum, logarithmic norm μ₁ for propagation,
  Duhamel bound for L vs its float rounding); entropy of the computed blocks enclosed in Arb, plus the
  continuity correction ε′ ln(D/ε′). Validated against `certify_finite.py` on all three M = 60 states
  (its intervals contain the Arb-exp ones). With `CERT_LD=1` the integration runs in x87 extended
  precision (u = 2⁻⁶⁴), needed for 2-Pauli at large M where μ₁ > 0 (≈ 10.6 at M = 180). The first version
  of that mode formed the step scalar h/n in float64 (error 2⁻⁵⁴, outside the 2⁻⁶⁴ budget); found in a
  self-audit and fixed (commit b23e357) before any extended-precision certificate entered this table;
  the M = 60 cross-check against Arb exp was repeated after the fix
  (`runs/certB_2p_M60_LDfix.log`: [4.41421550e-6, 4.41425605e-6] ⊇ the Arb-exp interval). Float64 mode is
  unaffected (its 4u term covers that rounding).

## Certificates

| channel | p (exact) | p | M | K | residues | state file | sha256 (first 16) | certifier | I_c interval (bits) |
|---|---|---|---|---|---|---|---|---|---|
| 2-Pauli | 605469/5000000 | 0.1210938 | 60 | 60 | mod 2 {0,1} | `runs/cert_2p_M60_p0.1210938.npy` | 1f5e1e9af7f55fea | finite (Arb exp) | [4.41423406e-6, 4.41423749e-6] |
| 2-Pauli | 605469/5000000 | 0.1210938 | 60 | 60 | mod 2 {0,1} | same | 1f5e1e9af7f55fea | big | [4.40445984e-6, 4.42401170e-6] |
| BB84 | 304297/2500000 | 0.1217188 | 60 | 60 | mod 2 {0,1} | `runs/cert_bb_M60_thr.npy` | 3eec80bc6be5f2eb | finite (Arb exp) | [6.34444364e-6, 6.34445247e-6] |
| BB84 | 304297/2500000 | 0.1217188 | 60 | 60 | mod 2 {0,1} | same | 3eec80bc6be5f2eb | big | [6.34280587e-6, 6.34609023e-6] |
| depolarizing | 652187/10000000 | 0.0652187 (3p = 0.1956561) | 60 | 60 | mod 3 {0,1} | `runs/cert_dep_M60_thr.npy` | 95b7b71c59885ac9 | finite (Arb exp) | [1.54889076e-6, 1.54890269e-6] |
| depolarizing | 652187/10000000 | 0.0652187 | 60 | 60 | mod 3 {0,1} | same | 95b7b71c59885ac9 | big | [1.54770191e-6, 1.55009153e-6] |
| BB84 | 619531/5000000 | 0.1239062 | 120 | 90 | mod 2 {0,1} | `runs/state_bb_M120_thr.npy` | 6248359d72a58e80 | big | [1.44239924e-6, 1.44929662e-6] |
| 2-Pauli | 308047/2500000 | 0.1232188 | 120 | 90 | mod 2 {0,1} | `runs/state_2p_M120_thr.npy` | c54db19b8afc4d82 | big | [1.73815583e-6, 2.46289307e-6] |
| depolarizing | 1051/16000 | 0.0656875 (3p = 0.1970625) | 120 | 90 | mod 3 {0,1} | `runs/state_dep_M120_thr.npy` | 2067c9f0cd871ed3 | big | [1.26617376e-6, 1.27106873e-6] |
| BB84 | 156211/1250000 | 0.1249688 | 120 | 90 | mod 3 {0,1} | `runs/state_bb3_M120_thr.npy` | e3e7a1aff721ad41 | big | [3.32197192e-6, 3.32992534e-6] |
| 2-Pauli | 154961/1250000 | 0.1239688 | 120 | 90 | mod 3 {0,1} | `runs/state_2p3_M120_thr.npy` | 65593ca6babd71c1 | big | [4.05046902e-6, 4.96028384e-6] |
| family r = 1/4 | 17827/80000 | 0.2228375 (total) | 60 | 60 | mod 3 {0,1} | `runs/state_pr0.25_k3_M60_thr.npy` | 1546bd790432b60b | big | [1.32667340e-6, 1.33009405e-6] |
| family r = 1/2 | 3359/16000 | 0.2099375 (total) | 60 | 60 | mod 3 {0,1} | `runs/state_pr0.5_k3b_M60_thr.npy` | 3a05874737beccad | big | [3.73849269e-6, 3.74092802e-6] |
| family r = 3/4 | 16127/80000 | 0.2015875 (total) | 60 | 60 | mod 3 {0,1} | `runs/state_pr0.75_k3_M60_thr.npy` | 0ea6c8003a9478ba | big | [1.63392704e-6, 1.63632317e-6] |
| depolarizing | 41207/625000 | 0.0659312 (3p = 0.1977936) | 240 | 110 | mod 3 {0,1} | `runs/state_dep_M240_thr.npy` | 9dffeb3cd424e413 | big | [1.32046252e-6, 1.33093482e-6] |
| BB84 | 63/500 | 0.126 | 240 | 110 | mod 3 {0,1} | `runs/state_bb3_M240_q0.126.npy` | b947cad42e44ebca | big | [5.54571320e-6, 5.56189919e-6] |
| 2-Pauli | 2003/16000 | 0.1251875 | 180 | 90 | mod 3 {0,1} | `runs/state_2p3_M180_thr.npy` | 29f739e0098b53a8 | big, `CERT_LD=1` | [1.17881571e-6, 1.19653069e-6] |
| 2-Pauli | 251/2000 | 0.1255 | 240 | 110 | mod 3 {0,1} | `runs/state_2p3_M240_p0.1255.npy` | 682e0a3d42e431d1 | big, `CERT_LD=1` | [1.88681426e-6, 2.70734623e-6] |

Logs: `runs/certF_*_M60_thr.log` (Arb exp; these predate the provenance lines, the state files are
unchanged, see hashes) and `runs/certB_*_M60.log`, `runs/certB_*_M120.log` (with argv, input sha256 and
python-flint version).

Commands:

    python3 certify_finite.py 2pauli 605469/5000000 runs/cert_2p_M60_p0.1210938.npy 60 2 0,1
    python3 certify_finite.py bb84 304297/2500000 runs/cert_bb_M60_thr.npy 60 2 0,1
    python3 certify_finite.py depol 652187/10000000 runs/cert_dep_M60_thr.npy 60 3 0,1
    python3 certify_big.py bb84 619531/5000000 runs/state_bb_M120_thr.npy 120 2 0,1
    python3 certify_big.py 2pauli 308047/2500000 runs/state_2p_M120_thr.npy 120 2 0,1
    python3 certify_big.py depol 1051/16000 runs/state_dep_M120_thr.npy 120 3 0,1
    python3 certify_big.py bb84 156211/1250000 runs/state_bb3_M120_thr.npy 120 3 0,1
    python3 certify_big.py 2pauli 154961/1250000 runs/state_2p3_M120_thr.npy 120 3 0,1
    python3 certify_big.py pr0.25 17827/80000 runs/state_pr0.25_k3_M60_thr.npy 60 3 0,1
    python3 certify_big.py pr0.5 3359/16000 runs/state_pr0.5_k3b_M60_thr.npy 60 3 0,1
    python3 certify_big.py pr0.75 16127/80000 runs/state_pr0.75_k3_M60_thr.npy 60 3 0,1
    python3 certify_big.py depol 41207/625000 runs/state_dep_M240_thr.npy 240 3 0,1
    python3 certify_big.py bb84 63/500 runs/state_bb3_M240_q0.126.npy 240 3 0,1
    CERT_LD=1 CERT_PROCS=4 python3 certify_big.py 2pauli 2003/16000 runs/state_2p3_M180_thr.npy 180 3 0,1
    CERT_LD=1 CERT_PROCS=4 python3 certify_big.py 2pauli 251/2000 runs/state_2p3_M240_p0.1255.npy 240 3 0,1

## Gaussian-limit (bosonic) certificates

`certify_limit.py` (Arb, Fock basis, output truncated at D with a flag state; see README §6.3–6.4):

| channel | parameter (exact) | code | state file | sha256 (first 16) | Nin, D | I_c interval (bits) |
|---|---|---|---|---|---|---|
| additive Gaussian noise | σ² = 2/5 (T = 5/7, G = 7/5) | mod 3 {0,1} | `runs/agn_k3_sigma2_0.4.npy` | 7917b7de2ccf2add | 140, 320 | [3.35715244e-6, 3.35738269e-6] |
| additive Gaussian noise | σ² = 81/200 (T = 200/281, G = 281/200) | mod 3 {0,1} | `runs/agn_k3_sigma2_0.405.npy` | 7e39e6ffc2144ccc | 300, 600 | [6.33597116e-8, 6.40767279e-8] |
| additive Gaussian noise | σ² = 81/200 | general rank 2 (k = 1) | `runs/agn_gen_sigma2_0.405.npy` | d9a9371dde9ac4ee | 200, 400 | [4.36362113e-7, 4.37080935e-7] |
| additive Gaussian noise | σ² = 407/1000 (T = 1000/1407, G = 1407/1000) | general rank 2 (k = 1) | `runs/agn_gen_sigma2_0.407.npy` | cb6c059e736810ac | 200, 400 | [5.43799903e-8, 5.50010769e-8] |
| additive Gaussian noise | σ² = 163/400 (T = 400/563, G = 563/400) | general rank 2 (k = 1) | `runs/agn_gen_sigma2_0.4075.npy` | 756c07a44641776b | 300, 600 | [5.97733793e-9, 6.88032591e-9] |
| additive Gaussian noise, stretched mod 3 (22 photons) | σ² = 51/125 (T = 125/176, G = 176/125) | mod 3 {0,1} | `runs/agn_k3s15_sigma2_0.408.npy` | e404a53c99f7edae | 300, 600 | [9.57696812e-8, 9.62308372e-8] |
| additive Gaussian noise, stretched mod 3 (22 photons) | σ² = 41/100 (T = 100/141, G = 141/100) | mod 3 {0,1} | `runs/agn_k3s15_sigma2_0.41.npy` | bf70a21fa8c3637e | 300, 600 | [1.85979392e-8, 1.92549153e-8] |
| additive Gaussian noise, stretched mod 3 (29 photons) | σ² = 411/1000 (T = 1000/1411, G = 1411/1000) | mod 3 {0,1} | `runs/agn_k3s2_sigma2_0.411.npy` | 58c42750d00e3e2d | 300, 600 | [1.17506701e-9, 1.58865651e-9] |
| additive Gaussian noise, stretched mod 3 (29 photons) | σ² = 3289/8000 (T = 8000/11289, G = 11289/8000) | mod 3 {0,1} | `runs/agn_k3s2_sigma2_0.411125.npy` | 4cec514bc0b35e4f | 300, 600 | [1.00509101e-9, 1.40751494e-9] |
| additive Gaussian noise, stretched mod 3 (29 photons) | σ² = 6579/16000 (T = 16000/22579, G = 22579/16000) | mod 3 {0,1} | `runs/agn_k3s2_sigma2_0.4111875.npy` | c7b1356a3fdff992 | 300, 600 | [9.61014601e-10, 1.32820991e-9] |
| additive Gaussian noise, stretched mod 3 (29 photons) | σ² = 823/2000 (T = 2000/2823, G = 2823/2000) | mod 3 {0,1} | `runs/agn_k3s2_sigma2_0.4115.npy` | 91359b8bf852b9b9 | 300, 600 | [2.47918646e-10, 6.12580391e-10] |
| additive Gaussian noise, stretched mod 3 (35 photons) | σ² = 4121/10000 (T = 10000/14121, G = 14121/10000) | mod 3 {0,1} | `runs/agn_k3s24_sigma2_0.4121.npy` | 6df867127041c4cf | 300, 600 | [2.28227285e-11, 3.70100931e-10] |
| thermal attenuator, ν = 1/4 | η = 313/500 (p = 187/500) | general rank 2 (k = 1) | `runs/th0.25_best.npy` | 8c1af8b44b2422ca | 160, 360 | [2.48378509e-7, 2.48746745e-7] |
| thermal attenuator, ν = 1/2 | η = 1393/2000 (p = 607/2000) | general rank 2 (k = 1) | `runs/th0.5_best.npy` | fb26e9b7951cefd0 | 160, 360 | [6.31451396e-8, 6.36306427e-8] |
| thermal attenuator, ν = 2 | η = 8587/10000 (p = 1413/10000) | general rank 2 (k = 1) | `runs/th2_best.npy` | 413fcd7ec3489aac | 160, 360 | [1.14709176e-7, 1.14892934e-7] |
| thermal attenuator, ν = 3 | η = 1791/2000 (p = 209/2000) | general rank 2 (k = 1) | `runs/th3_best.npy` | 6cd8d3af49f65888 | 160, 360 | [8.31897855e-8, 8.33853440e-8] |
| thermal attenuator, ν = 1/4, larger cutoff | η = 1251/2000 (p = 749/2000; T = 5004/8749, G = 8749/8000) | general rank 2 (k = 1) | `runs/thgen025_N400_s0.3745.npy` | 47f7a1772b80acd8 | 400, 800 | [4.73283428e-8, 4.89118794e-8] |
| thermal attenuator, ν = 1/2, larger cutoff | η = 6961/10000 (p = 3039/10000; T = 13922/23039, G = 23039/20000) | general rank 2 (k = 1) | `runs/thgen05_N240_s0.3039.npy` | c475267b5d317b4c | 240, 480 | [9.79254604e-9, 1.05297073e-8] |
| thermal attenuator, ν = 2, larger cutoff | η = 3433/4000 (p = 567/4000; T = 3433/5134, G = 2567/2000) | general rank 2 (k = 1) | `runs/thgen2_N280_s0.14175.npy` | e00251e7e14a8da9 | 280, 560 | [3.83709010e-8, 3.86723241e-8] |
| thermal attenuator, ν = 3, larger cutoff | η = 17907/20000 (p = 2093/20000; T = 17907/26279, G = 26279/20000) | general rank 2 (k = 1) | `runs/thgen3_N320_s0.10465.npy` | b5374ff708ab43e4 | 320, 640 | [1.42721100e-7, 1.43166597e-7] |
| thermal attenuator, ν = 5 | η = 7459/8000 (p = 541/8000; T = 7459/10705, G = 2141/1600) | mod 3 {0,1} | `runs/thnu_5_big_best.npy` | 2b85f0323ce13231 | 260, 520 | [1.17390097e-6, 1.17400578e-6] |
| thermal attenuator, ν = 5, larger cutoff | η = 14913/16000 (p = 1087/16000; T = 4971/7145, G = 4287/3200) | mod 3 {0,1} | `runs/thnu5_N360_eta0.9320625.npy` | 2fabfd89705e91e9 | 360, 720 | [4.44226553e-8, 4.46329925e-8] |
| emergent depolarizing (thermal attenuator, ν = 1/(2η)) | η = 3756117/5120000 (p = 1363883/20480000 per Pauli, 3p = 0.199787548828125; T = 14108414917689/22722859520000, G = 8876117/7512234) | mod 3 {0,1}, mean 57.5 photons | `../depolarizing-threshold/runs/hpcont57b_eta0.7336166.npy` | 8d735c0f9da97002 | 240, 360 | [4.99212570e-10, 8.85924622e-10] |
| emergent depolarizing | η = 939063/1280000 (p = 340937/5120000, 3p = 0.1997677734375) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont57b_eta0.7336430.npy` | bdb0a58f083b0bdb | 240, 360 | [4.03487761e-10, 7.45272255e-10] |
| emergent depolarizing | η = 234777/320000 (p = 85223/1280000, 3p = 0.19974140625) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont57b_eta0.7336781.npy` | 68f76503fd738bc6 | 240, 360 | [5.24050621e-10, 9.41855682e-10] |
| emergent depolarizing | η = 29349/40000 (p = 10651/160000, 3p = 0.19970625) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont57b_eta0.7337250.npy` | 8fdea697037bd795 | 240, 360 | [1.43547608e-9, 1.77347936e-9] |
| emergent depolarizing | η = 587/800 (p = 213/3200, 3p = 0.1996875) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont57_eta0.7337500.npy` | c6fadf5fc86bca2b | 240, 360 | [1.04727317e-9, 1.45793917e-9] |
| emergent depolarizing | η = 23481/32000 (p = 8519/128000 per Pauli, 3p = 0.1996640625; T = 551357361/887696000, G = 55481/46962) | mod 3 {0,1}, mean 57.5 photons | `../depolarizing-threshold/runs/hpcont57_eta0.7337813.npy` | 78fe258b720f161f | 240, 360 | [8.36623805e-10, 1.20960365e-9] |
| emergent depolarizing | η = 11741/16000 (p = 4259/64000, 3p = 0.199640625; T = 137851081/221928000, G = 27741/23482) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont57_eta0.7338125.npy` | 0d305290c225c092 | 240, 360 | [3.32999404e-10, 7.46110460e-10] |
| emergent depolarizing | η = 5871/8000 (p = 2129/32000, 3p = 0.19959375; T = 34468641/55484000, G = 13871/11742) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont57_eta0.7338750.npy` | 24c40426445e91a7 | 240, 360 | [1.81914406e-9, 2.22028303e-9] |
| emergent depolarizing | η = 11743/16000 (p = 4257/64000 per Pauli, 3p = 0.199546875; T = 137898049/221944000, G = 27743/23486) | mod 3 {0,1}, mean 57.5 photons | `../depolarizing-threshold/runs/hpcont57_eta0.7339375.npy` | 609cfc3e72d4bd11 | 240, 360 | [2.57112611e-9, 2.97906949e-9] |
| emergent depolarizing | η = 367/500 (p = 133/2000, 3p = 0.1995; T = 134689/216750, G = 867/734) | mod 3 {0,1}, mean 57.5 photons | `../depolarizing-threshold/runs/s57_eta0.734_best.npy` | c730e18c56981346 | 240, 360 | [6.64143550e-10, 1.00461950e-9] |
| emergent depolarizing | η = 7341/10000 (p = 2659/40000 per Pauli, 3p = 0.199425; T = 53890281/86705000, G = 17341/14682) | mod 3 {0,1}, mean 38.5 photons | `../depolarizing-threshold/runs/hpcont38_eta0.7341000.npy` | f389ca632e62cd5f | 240, 360 | [4.69714037e-8, 4.72320460e-8] |
| emergent depolarizing | η = 29369/40000 (p = 10631/160000, 3p = 0.19933125; T = 862538161/1387380000, G = 69369/58738) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont38_eta0.7342250.npy` | fda1524d3ff1acf3 | 240, 360 | [4.50609099e-8, 4.53214898e-8] |
| emergent depolarizing | η = 14687/20000 (p = 5313/80000, 3p = 0.1992375; T = 215707969/346870000, G = 34687/29374) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont38_eta0.7343500.npy` | 98f3035d24d6a29f | 240, 360 | [2.36712801e-8, 2.40039012e-8] |
| emergent depolarizing | η = 3673/5000 (p = 1327/20000 per Pauli, 3p = 0.19905; T = 13490929/21682500, G = 8673/7346) | mod 3 {0,1}, mean 38.5 photons | `../depolarizing-threshold/runs/stretch2b_eta0.7346.npy` | 1b44cc551b74c8c1 | 240, 360 | [1.52256332e-8, 1.54665840e-8] |
| emergent depolarizing | η = 2351/3200 (p = 849/12800 per Pauli, 3p = 0.198984375; T = 5527201/8881600, G = 5551/4702) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont_eta0.7346875.npy` | c4f3ae62d4e8ce2b | 240, 360 | [1.82651528e-8, 1.85766717e-8] |
| emergent depolarizing | η = 2939/4000 (p = 1061/16000, 3p = 0.1989375; T = 8637721/13878000, G = 6939/5878) | mod 3 {0,1} | `../depolarizing-threshold/runs/hpcont_eta0.7347500.npy` | d01d54a7afcd92b9 | 240, 360 | [1.07683723e-7, 1.07994494e-7] |
| thermal attenuator, ν = 5 (superseded) | η = 37307/40000 (p = 2693/40000) | mod 3 {0,1} | `runs/th5_best.npy` | 3a0ad4e6310ec0cf | 160, 360 | [2.50657502e-6, 2.50662460e-6] |
| additive Gaussian noise, squeezed cat | σ² = 49/125 | mod 3 {0,1} (residues of D(α)S(r) on vacuum) | `runs/agn_sqcat_sigma2_0.392.npy` | eb13a99ef8255ad0 | 140, 320 | [1.73865751e-6, 1.73880842e-6] |
| thermal attenuator, ν = 1 | η = 39/50 (p = 11/50; T = 39/61, G = 61/50) | general rank 2 (k = 1) | `runs/th1_gen_eta0.78.npy` | 6c91dfe91ec43e9f | 120, 300 | [2.12631890e-8, 2.16095917e-8] |
| thermal attenuator, ν = 1 | η = 3133/4000 (p = 867/4000; T = 3133/4867, G = 4867/4000) | mod 3 {0,1} | `runs/th1_k3_p0.21675.npy` | eff871b4adcf91ce | 80, 240 | [8.27408111e-6, 8.27425105e-6] |
| thermal attenuator, ν = 1, stretched mod 3 (30 photons) | η = 12461/16000 (p = 3539/16000; T = 12461/19539, G = 19539/16000) | mod 3 {0,1} | `runs/th1_s15_eta0.7788125.npy` | d974f5456576354b | 240, 480 | [2.44659596e-9, 2.90742624e-9] |
| thermal attenuator, ν = 1, stretched mod 3 (30 photons) | η = 623/800 (p = 177/800; T = 623/977, G = 977/800) | mod 3 {0,1} | `runs/th1_s15N32_eta0.77875.npy` | bcd0016475563bc1 | 320, 640 | [2.58920431e-10, 7.64931369e-10] |
| thermal amplifier (emergent 2-Pauli) | p = 16/125 (T = 93/125, G = 93/61; τ ≈ 1.1343, N = 93/32) | general rank 2 (k = 1) | `runs/amp_gen_p0.128.npy` | 9fb96e4e40e7e33a | 200, 400 | [1.43395294e-7, 1.43876228e-7] |
| thermal amplifier (emergent 2-Pauli) | p = 51/400 (T = 149/200, G = 149/98; τ ≈ 1.1327, N = 149/51) | mod 3 {0,1} | `runs/amp_2p3_p0.1275.npy` | 57ba63a1717f0653 | 140, 320 | [9.10836870e-7, 9.11017611e-7] |
| emergent 2-Pauli channel | p = 31/250 (T = 94/125, G = 94/63) | mod 4 {0,2} | `runs/2p_k4_p0.124.npy` | b7ffd6aa8ad2f999 | 160, 280 | [1.96141999e-4, 1.96142031e-4] |

    python3 certify_limit.py agn 2/5 runs/agn_k3_sigma2_0.4.npy 0,1 3 140 320
    python3 certify_limit.py agn 81/200 runs/agn_k3_sigma2_0.405.npy 0,1 3 300 600
    python3 certify_limit.py agn 81/200 runs/agn_gen_sigma2_0.405.npy 0,0 1 200 400
    python3 certify_limit.py agn 407/1000 runs/agn_gen_sigma2_0.407.npy 0,0 1 200 400
    python3 certify_limit.py agn 163/400 runs/agn_gen_sigma2_0.4075.npy 0,0 1 300 600
    python3 certify_limit.py agn 51/125 runs/agn_k3s15_sigma2_0.408.npy 0,1 3 300 600
    python3 certify_limit.py agn 41/100 runs/agn_k3s15_sigma2_0.41.npy 0,1 3 300 600
    python3 certify_limit.py agn 411/1000 runs/agn_k3s2_sigma2_0.411.npy 0,1 3 300 600
    python3 certify_limit.py agn 3289/8000 runs/agn_k3s2_sigma2_0.411125.npy 0,1 3 300 600
    python3 certify_limit.py agn 6579/16000 runs/agn_k3s2_sigma2_0.4111875.npy 0,1 3 300 600
    python3 certify_limit.py agn 823/2000 runs/agn_k3s2_sigma2_0.4115.npy 0,1 3 300 600
    python3 certify_limit.py agn 4121/10000 runs/agn_k3s24_sigma2_0.4121.npy 0,1 3 300 600
    python3 certify_limit.py 2pauli 31/250 runs/2p_k4_p0.124.npy 0,2 4 160 280
    python3 certify_limit.py 2pauli 51/400 runs/amp_2p3_p0.1275.npy 0,1 3 140 320
    python3 certify_limit.py 2pauli 16/125 runs/amp_gen_p0.128.npy 0,0 1 200 400
    python3 certify_limit.py th1 867/4000 runs/th1_k3_p0.21675.npy 0,1 3 80 240
    python3 certify_limit.py th1 11/50 runs/th1_gen_eta0.78.npy 0,0 1 120 300
    python3 certify_limit.py th1 3539/16000 runs/th1_s15_eta0.7788125.npy 0,1 3 240 480
    python3 certify_limit.py th1 177/800 runs/th1_s15N32_eta0.77875.npy 0,1 3 320 640
    python3 certify_limit.py th0.25 187/500 runs/th0.25_best.npy 0,0 1 160 360
    python3 certify_limit.py th0.5 607/2000 runs/th0.5_best.npy 0,0 1 160 360
    python3 certify_limit.py th2 1413/10000 runs/th2_best.npy 0,0 1 160 360
    python3 certify_limit.py th3 209/2000 runs/th3_best.npy 0,0 1 160 360
    python3 certify_limit.py th0.25 749/2000 runs/thgen025_N400_s0.3745.npy 0,0 1 400 800
    python3 certify_limit.py th0.5 3039/10000 runs/thgen05_N240_s0.3039.npy 0,0 1 240 480
    python3 certify_limit.py th2 567/4000 runs/thgen2_N280_s0.14175.npy 0,0 1 280 560
    python3 certify_limit.py th3 2093/20000 runs/thgen3_N320_s0.10465.npy 0,0 1 320 640
    python3 certify_limit.py th5 2693/40000 runs/th5_best.npy 0,1 3 160 360
    python3 certify_limit.py th5 541/8000 runs/thnu_5_big_best.npy 0,1 3 260 520
    python3 certify_limit.py th5 1087/16000 runs/thnu5_N360_eta0.9320625.npy 0,1 3 360 720
    python3 certify_limit.py depol 1363883/20480000 ../depolarizing-threshold/runs/hpcont57b_eta0.7336166.npy 0,1 3 240 360
    python3 certify_limit.py depol 340937/5120000 ../depolarizing-threshold/runs/hpcont57b_eta0.7336430.npy 0,1 3 240 360
    python3 certify_limit.py depol 85223/1280000 ../depolarizing-threshold/runs/hpcont57b_eta0.7336781.npy 0,1 3 240 360
    python3 certify_limit.py depol 10651/160000 ../depolarizing-threshold/runs/hpcont57b_eta0.7337250.npy 0,1 3 240 360
    python3 certify_limit.py depol 213/3200 ../depolarizing-threshold/runs/hpcont57_eta0.7337500.npy 0,1 3 240 360
    python3 certify_limit.py depol 8519/128000 ../depolarizing-threshold/runs/hpcont57_eta0.7337813.npy 0,1 3 240 360
    python3 certify_limit.py depol 4259/64000 ../depolarizing-threshold/runs/hpcont57_eta0.7338125.npy 0,1 3 240 360
    python3 certify_limit.py depol 2129/32000 ../depolarizing-threshold/runs/hpcont57_eta0.7338750.npy 0,1 3 240 360
    python3 certify_limit.py depol 4257/64000 ../depolarizing-threshold/runs/hpcont57_eta0.7339375.npy 0,1 3 240 360
    python3 certify_limit.py depol 133/2000 ../depolarizing-threshold/runs/s57_eta0.734_best.npy 0,1 3 240 360
    python3 certify_limit.py depol 2659/40000 ../depolarizing-threshold/runs/hpcont38_eta0.7341000.npy 0,1 3 240 360
    python3 certify_limit.py depol 10631/160000 ../depolarizing-threshold/runs/hpcont38_eta0.7342250.npy 0,1 3 240 360
    python3 certify_limit.py depol 5313/80000 ../depolarizing-threshold/runs/hpcont38_eta0.7343500.npy 0,1 3 240 360
    python3 certify_limit.py depol 1327/20000 ../depolarizing-threshold/runs/stretch2b_eta0.7346.npy 0,1 3 240 360
    python3 certify_limit.py depol 849/12800 ../depolarizing-threshold/runs/hpcont_eta0.7346875.npy 0,1 3 240 360
    python3 certify_limit.py depol 1061/16000 ../depolarizing-threshold/runs/hpcont_eta0.7347500.npy 0,1 3 240 360
    python3 certify_limit.py agn 49/125 runs/agn_sqcat_sigma2_0.392.npy 0,1 3 140 320

The AGN entry was checked by an independent agent with its own implementation (diffusion Lindbladian, no
repository code): I_c = 3.3572969×10⁻⁶ bits, inside the interval (README §6.4). The seven stretched AGN
entries were checked by another from-scratch implementation with five routes (loss then amplifier, amplifier
then loss, the diffusion semigroup by uniformization and by expm, unitary dilations; cutoffs 600–800):
9.6074153×10⁻⁸, 1.9038668×10⁻⁸, 1.4606876×10⁻⁹, 1.2823360×10⁻⁹, 1.2180839×10⁻⁹, 4.9355320×10⁻¹⁰ and 2.7278801×10⁻¹⁰ bits, all inside. The larger-cutoff ν = 1/4 entry gives 4.8512689×10⁻⁸ bits, the ν = 1/2 entry gives 1.0321304×10⁻⁸ bits, the ν = 5 entry gives 4.4564070×10⁻⁸ bits and the ν = 2 and ν = 3 entries give 3.8572287×10⁻⁸ and 1.4302399×10⁻⁷ bits with the same code, inside.
The two stretched ν = 1 entries
were checked with the same code (four routes, including the thermal beam splitter and the thermal master
equation): 2.7722611×10⁻⁹ and 6.3155836×10⁻¹⁰ bits, inside. The emergent
depolarizing entries and the ν = 5 entry (η = 7459/8000) were checked by another from-scratch
implementation with three routes (pure loss then amplifier, thermal master equation, beam splitter with a
Fock-mixed thermal environment; output cutoffs 420/560, tails < 3×10⁻⁷⁵): 7.472121×10⁻¹⁰ (η = 3756117/5120000), 6.313772×10⁻¹⁰ (939063/1280000),
7.921672×10⁻¹⁰ (234777/320000), 1.6618809×10⁻⁹ (0.733725), 1.3148727×10⁻⁹ (0.73375),
1.0813608×10⁻⁹ (0.73378125), 6.061502×10⁻¹⁰ (0.7338125), 2.0771682×10⁻⁹ (0.733875),
2.8369733×10⁻⁹ (0.7339375), 8.794423×10⁻¹⁰ (0.734), 4.7133075×10⁻⁸ (0.7341), 4.5225625×10⁻⁸ (0.734225),
2.3876871×10⁻⁸ (0.73435), 1.5376335×10⁻⁸ (0.7346), 1.8455798×10⁻⁸, 1.0787368×10⁻⁷ and 1.1739661×10⁻⁶ bits, routes agreeing to ≤ 7.5×10⁻¹⁴, all inside (README §7, §6.6). Intervals from the
Weyl/Ostrowski enclosure depend on the float eigenbasis, so another LAPACK build can print slightly
different (equally valid) endpoints.

## Independent review

A separate agent (fresh context) tried to break the three M = 60 certificates. It re-derived the frame
change and Bloch factors, checked exp(L) = N including the negative-dephasing case, compared the generator
line by line with qutip.piqs, verified the band evolution at M = 60 against an exact Schur–Weyl formula for
product inputs (agreement ~1e-14 in every band and sector), compared its own 2^M brute force with both
certifiers for M ≤ 10 (brute-force value inside the Arb interval in all cases), checked the exact-zero
residue pattern of the stored states, and fuzzed `entropy_interval` (219 cases, none unsound). Verdict:
all three VALID; remaining trust assumption: correctness of python-flint's `arb_mat.exp`.
Minor findings, all fixed afterwards: intervals had been quoted rounded to nearest (now outward);
`entropy_interval` could return NaN for balls containing 0 (fail-safe; now uses |x|⁺²); logs lacked
provenance (now printed).

A second independent agent (fresh context, 12 Oct 2026) reviewed `certify_big.py`, the certifier behind
every M ≥ 120 entry, line by line against its error budget: Taylor truncation, rounding of the sparse
products (γ_q, no duplicate CSR entries), summation, propagation with μ₁, the Duhamel model error, the
assembly rounding and weighting, the trace-norm/Mirsky/Fannes chain with D summed over sectors, the Z lower
bound and the direct-sum/multiplicity bookkeeping. Experiments: its own from-scratch 2^M brute force for
13 random residue codes (M = 4…10, all channel types), always inside the intervals of both certifiers in
both precisions; certify_big vs certify_finite at M = 24, 30; the true integration error against an Arb
reference for 2-Pauli bands with μ₁ up to 76 and ‖L‖₁ up to 1381, at most 1% of the claimed bound in
float64 (29% in extended precision); certify_finite on the family r = 1/4 point (previously certified only
by `certify_big`): [1.3283821, 1.3283854]×10⁻⁶, inside the published interval. Verdict: float64 mode sound
as implemented, including the older code that produced the earlier logs; extended-precision mode sound on
this host after the step-scalar fix. Its suggestions, all applied before the final runs: refuse
`CERT_LD=1` unless `np.longdouble` has a 64-bit mantissa (on platforms where it is plain double the
budget would be 4–6× too small, demonstrated); an explicit error instead of an `assert` for the
continuity-bound range; reject codes with weight off their residues or K > M (the lower bound survives
such inputs by a twirling argument, the upper end of the interval does not); re-run every certificate with
the final code so the manifest is reproducible bit for bit (`recertify_all.sh`, logs in `runs/certB_*.log`;
status: all re-runs are done, M = 60 to 240 including both `CERT_LD=1` entries, all still positive,
intervals at most ~4×10⁻¹³ wider, and the table above carries the re-run intervals).

Caveat (also from the review): BB84 and depolarizing channels are monotone in p by composition, so
Q > 0 holds for every smaller p as well. For the 2-Pauli family no such composition exists, so
"p_th ≥ 0.1255" holds in the sup sense; it does not by itself show Q > 0 at every smaller p.
