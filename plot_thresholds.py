"""Certified threshold lower bounds versus the number of qubits M (values from CERTIFICATES.md)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (M, p) of each certified entry; best code per M
BB84 = [(60, 0.1217188), (120, 0.1249688), (240, 0.126)]
TWO_PAULI = [(60, 0.1210938), (120, 0.1239688), (180, 0.1251875), (240, 0.1255)]
DEPOL_3P = [(60, 0.1956561), (120, 0.1970625), (240, 0.1977936)]

fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4))
for pts, name, c in [(BB84, "BB84 (independent X–Z)", "C0"), (TWO_PAULI, "2-Pauli", "C1")]:
    a.plot(*zip(*pts), "o-", color=c, label=f"{name}, certified (this work)")
a.axhline(0.118371, color="C0", ls="--", lw=1, label="previous record BB84 (n = 30)")
a.axhline(0.118067, color="C1", ls="--", lw=1, label="previous record 2-Pauli (n = 30)")
a.axhline(0.110028, color="C0", ls=":", lw=1, label="hashing BB84")
a.axhline(0.113546, color="C1", ls=":", lw=1, label="hashing 2-Pauli")
a.set_xscale("log")
a.minorticks_off()
a.set_xticks([60, 120, 180, 240], ["60", "120", "180", "240"])
a.set_xlabel("qubits M of the explicit input state")
a.set_ylabel("certified lower bound on p_th")
a.legend(fontsize=7, loc="lower right")
a.set_title("Q > 0 certified up to")

b.plot(*zip(*DEPOL_3P), "o-", color="C2", label="certified, explicit M (this work)")
b.axhline(0.2029161, color="k", ls="--", lw=1, label="asymptotic record (Gaussian limit, floating point)")
b.axhline(0.199425, color="C3", ls="-.", lw=1, label="Gaussian limit, I_c certified (this work)")
b.axhline(0.194868, color="C2", ls="--", lw=1, label="previous explicit certified (n = 45)")
b.axhline(0.18929, color="C2", ls=":", lw=1, label="hashing")
b.set_xscale("log")
b.minorticks_off()
b.set_xticks([60, 120, 240], ["60", "120", "240"])
b.set_xlabel("qubits M of the explicit input state")
b.set_ylabel("certified lower bound on 3p_th")
b.legend(fontsize=7, loc="lower right")
b.set_title("Depolarizing channel")

fig.tight_layout()
fig.savefig("thresholds.png", dpi=150)
