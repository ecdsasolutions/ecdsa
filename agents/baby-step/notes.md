Height 40 bits:
- Solved practice instance (p=792318143303, n=396158503739) in 0.71s wall clock, k=80157930310. Verified k*G == P.
- Negation BSGS stride M = 2m + 1, m = isqrt(n)//2 + 1 (~315k baby steps).
- Dict stores x -> (j << 1) | (y & 1).
- Direct pow(dx, -1, p) (CPython C-level egcd) achieves ~600k-750k point additions per second, outperforming Python-level Montgomery batch inversion at 40 bits due to zero interpreter overhead.
- Added zero-overhead `try: pow(...) except ValueError:` for exact multiples of M.
- Verified against 15 boundary cases (1, 2, m±1, M±1, 2M±1, n-1, n//2).

Plan for 44 bits:
- At 44 bits (n ~ 1.7e13), m = isqrt(n)//2 + 1 is ~1.0 - 1.5 million entries.
- Python dict for 1.5M entries takes ~100MB RAM, which should fit within sandbox memory limit.
- Expected runtime for 44 bits with negation BSGS is ~3-5 seconds.
- At 48 bits, memory reaches >6M entries (causing MemoryError with dict); will need open-addressing array('I')/array('H') packed tables or parallel Pollard rho with distinguished points.
