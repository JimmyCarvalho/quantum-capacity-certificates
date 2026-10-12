# Quantum capacity certificates

Computer-assisted proofs that several standard noise models have **positive quantum capacity**
(Q > 0) at noise levels where this was not known before. Each claim is a statement
I_c(ρ, N) > 0 for an explicit input state ρ, enclosed by a rigorous interval computed in ball arithmetic
(Arb, via python-flint). Since Q ≥ I_c, a positive interval proves Q > 0.

Paper: [`paper/main.pdf`](paper/main.pdf) (J. K. E. Carvalho). Full work log: [`REPORT.md`](REPORT.md).
Every certificate, with exact parameters, state file, SHA-256 and command: [`CERTIFICATES.md`](CERTIFICATES.md).

## Results

| channel | previously known | certified here |
|---|---|---|
| additive Gaussian noise (σ², vacuum variance ½) | Q > 0 for σ² < 1/e = 0.3679 (Holevo–Werner 2001); GKP codes reach σ ≈ 0.598 | **Q > 0 for σ² ≤ 163/400 = 0.4075** (σ ≤ 0.6384) |
| thermal attenuator, ν = 1 thermal photon | Q > 0 for η ≥ 0.7841 (Mele et al., Jul 2026) | **Q > 0 for η ≥ 0.78** |
| thermal attenuator, ν = 1/4 … 5 | Gaussian-input thresholds only (ν ≠ 1) | **Q > 0 inside the open window at every ν tried, closing 22–50% of it** (`attenuator_map.png`) |
| thermal amplifier, gain 1.1343, 2.906 thermal photons | Gaussian inputs give 0 | **Q > 0** |
| independent X–Z (BB84) Pauli channel | q_th ≥ 0.118371 (May 2026) | **q_th ≥ 0.126** (explicit 240-qubit state) |
| 2-Pauli channel | p_th ≥ 0.118067 (May 2026) | **p_th ≥ 0.1255** (explicit 240-qubit state) |
| depolarizing channel | 3p_th ≥ 0.2029 (asymptotic, floating point); best explicit certified 0.194868 | best explicit certified **0.1977936** (240 qubits); in the Gaussian limit **0.19905 with I_c certified** (with the convergence theorem of Ahmed et al., the largest bound we know whose final sign is certified); the 0.2029 record is *not* beaten |

The additive-Gaussian-noise result also implies Q > 0 for every phase-insensitive Gaussian channel with
(Y + |1−τ|/2)/τ ≤ 0.4075, which goes beyond the Gaussian-input thresholds of thermal attenuators and
amplifiers with more than 4.23 thermal photons.

Physical picture: about 60% of the gain over 1/e on the additive-noise channel already comes from a
three-parameter **squeezed cat** (the n ≡ 0, 1 mod 3 parts of a displaced phase-squeezed state), certified up
to σ² = 0.392.

The key ingredients are the observation that symmetric states of many qubits under a Pauli channel behave
like a single bosonic mode under a Gaussian channel (generalising Ahmed et al., arXiv:2609.39747), rank-two
residue ("cat-like") codes, and a small deliberate breaking of their Z₃ symmetry, which multiplies the
coherent information several-fold near threshold.

## Reproduce

Requirements (tested): Python 3.13, numpy 2.5, scipy 1.18, python-flint 0.9.0; qutip 5.3 and matplotlib
only for cross-checks and the figure.

    pip install -r requirements.txt
    ./verify_certificates.sh quick     # Gaussian-limit certificates and M <= 60 (minutes)
    ./verify_certificates.sh full      # everything, including M = 120-240 (many CPU-hours)

The script checks the SHA-256 of each state file against `CERTIFICATES.md` and prints the certified
interval of each claim. Single certificates, e.g.

    cd pauli-emergence
    python3 certify_limit.py agn 163/400 runs/agn_gen_sigma2_0.4075.npy 0,0 1 300 600
    python3 certify_limit.py th1 11/50 runs/th1_gen_eta0.78.npy 0,0 1 120 300

Layout: `pauli-emergence/` (finite-M Dicke-basis evolution, optimisers, certifiers `certify_finite.py`,
`certify_big.py`, `certify_limit.py`, input states and logs in `runs/`), `depolarizing-threshold/`
(Gaussian-channel tools and the Arb Fock-space certifier `certify.py`), `paper/`.

## Status and disclosure

The computations, code and first draft were produced with the assistance of an AI system (Claude,
Anthropic); the author takes responsibility for the content. The certifiers were reviewed, and the main
certificates re-computed with independent code, by separate instances of the same AI system (summaries in
`CERTIFICATES.md` and `REPORT.md`). **Independent human verification is pending**; issues and checks are
welcome.

## License

Code: MIT (`LICENSE`). Text, figures and data (paper, reports, state files, logs): CC BY 4.0
(`LICENSE-CC-BY-4.0.md`).
