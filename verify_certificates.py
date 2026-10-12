"""Re-run the certificates listed in CERTIFICATES.md and check the input files.

usage: python3 verify_certificates.py quick|full

quick: hash checks, all Gaussian-limit certificates (certify_limit.py) and certify_big.py with M <= 60.
full:  everything, including certify_finite.py (Arb matrix exponential, ~40 min each at M = 60) and
       certify_big.py at M = 120-240 (CPU-hours; extended-precision runs need x87 long double).
"""
import hashlib
import os
import re
import shlex
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(ROOT, "pauli-emergence")


def hash_checks(text):
    bad = 0
    seen = set()
    for line in text.splitlines():
        if not line.startswith("|"):
            continue
        path = re.search(r"`(runs/[^`]+\.npy)`", line)
        digest = re.search(r"\|\s*([0-9a-f]{16})\s*\|", line)
        if not (path and digest) or path.group(1) in seen:
            continue
        seen.add(path.group(1))
        full = os.path.join(CODE, path.group(1))
        got = hashlib.sha256(open(full, "rb").read()).hexdigest()[:16] if os.path.exists(full) else "missing"
        ok = got == digest.group(1)
        bad += not ok
        print(f"{'ok ' if ok else 'BAD'} {path.group(1)} {got}")
    return bad


def commands(text, mode):
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not (s.startswith("python3 certify") or s.startswith("CERT_LD=1")):
            continue
        args = shlex.split(s)
        env = {}
        while "=" in args[0]:
            k, v = args.pop(0).split("=", 1)
            env[k] = v
        script = args[1]
        if mode == "quick":
            if script == "certify_finite.py":
                continue
            if script == "certify_big.py" and int(args[5]) > 60:
                continue
        out.append((env, args))
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "quick"
    text = open(os.path.join(ROOT, "CERTIFICATES.md")).read()
    print("== input files")
    bad = hash_checks(text)
    print("== certificates")
    failed = 0
    for env, args in commands(text, mode):
        e = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", **env)
        r = subprocess.run(args, cwd=CODE, env=e, capture_output=True, text=True)
        res = [l for l in r.stdout.splitlines() if "Ic in" in l]
        ok = "CERTIFIED POSITIVE" in r.stdout
        failed += not ok
        print(("ok  " if ok else "FAIL") + " " + " ".join(args[1:]) + "\n     " + (res[-1] if res else r.stderr[-300:]))
    print(f"== {bad} hash mismatches, {failed} certificates not reproduced")
    sys.exit(1 if bad or failed else 0)


if __name__ == "__main__":
    main()
