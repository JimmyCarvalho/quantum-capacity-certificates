"""Check qutip.piqs conventions: local emission+pumping (+dephasing) Liouvillian vs. brute-force Pauli channel."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, itertools
from qutip import piqs, Qobj
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm

def pauli_channel_bruteforce(rho, probs, N):
    I = np.eye(2); X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1])
    P = [I, X, Y, Z]
    out = rho.astype(complex)
    for site in range(N):
        new = np.zeros_like(out)
        for a in range(4):
            op = np.array([[1.0]])
            for s in range(N):
                op = np.kron(op, P[a] if s == site else I)
            new += probs[a] * op @ out @ op.conj().T
        out = new
    return out

def single_qubit_lambdas(gE, gP, gD, t):
    d = piqs.Dicke(1, emission=gE, pumping=gP, dephasing=gD)
    L = d.liouvillian().full()
    U = expm(L * t)
    # basis: dicke basis for N=1 is |j=1/2,m=1/2>, |j=1/2,m=-1/2> ; |m=+1/2> = excited?
    lam = {}
    for name, op in [("x", np.array([[0, 1], [1, 0]])), ("y", np.array([[0, -1j], [1j, 0]])), ("z", np.diag([1, -1]))]:
        v = U @ (op / 2).reshape(-1, order="F")  # qutip superoperators act on column-stacked vectors
        lam[name] = np.real(np.trace(op @ v.reshape(2, 2, order="F")))
    return lam

if __name__ == "__main__":
  print("emission=pumping=1, t=0.3 :", single_qubit_lambdas(1, 1, 0, 0.3), " expect lam_perp=e^{-t*?}")
  print("dephasing=1, t=0.3       :", single_qubit_lambdas(0, 0, 1, 0.3))
