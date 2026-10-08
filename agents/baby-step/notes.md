Height 36 bits:
- Solved practice instance (p=44382896699, n=4438300897) in 0.12s, k=1664885196.
- Negation BSGS stride M = 2*m + 1, m = isqrt(n)//2 + 1 (~33k entries).
- Baby dict maps x -> (j << 1) | (y & 1) saving object allocation and lookup overhead.
- Inversion using pow(d, -1, p) (extended Euclidean algorithm in CPython) is fast and reliable.

Plan for 40 bits:
- At 40 bits (n ~ 1.1e12), m = isqrt(n)//2 + 1 is ~520,000 entries. Dict fits in ~40MB RAM.
- Expected runtime for 40 bits with this BSGS implementation is ~1.5 - 2.0s.
- If necessary to speed up, Montgomery batch inversion in baby steps or giant steps can shave 20-30% off, though pure pow(..., -1, p) in Python may already be fast enough.
- For 44+ bits (m > 2-3 million entries, potential MemoryError in sandbox), prepare van Oorschot-Wiener parallel Pollard rho or kangaroo with distinguished points.
