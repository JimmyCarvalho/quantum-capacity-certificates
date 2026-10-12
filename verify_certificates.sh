#!/bin/bash
# Re-run the certificates of CERTIFICATES.md (quick or full); see verify_certificates.py.
cd "$(dirname "$0")"
exec python3 verify_certificates.py "${1:-quick}"
