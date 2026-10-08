Height 44 bits:
- Negation BSGS with direct egcd pow(dx, -1, p) solved 44-bit practice instance (p=15794928898921, n=15794926521203) in 5.35s wall clock, k=4322665950283.
- Table size m = isqrt(n)//2 + 1 = 1,987,142 entries fits easily in memory (~120MB).
- Direct pow in CPython yields ~520k baby steps/s and ~720k giant steps/s.
- Zero-cost try/except ValueError handles division by zero in giant steps without per-iteration branching.

Plan for 48 bits:
- At 48 bits, n ~ 2^48, m = isqrt(n)//2 + 1 can be up to 8.4M entries.
- Python dict for >6.5M entries causes MemoryError in sandbox.
- Must use packed open-addressing table:
  Option A: 2^23 or 2^24 slot open-addressing table using array('I') for (j<<1)|parity and array('H') for fingerprints (as described in book entry by Hare/Kangaroo for 48/52 bit).
  Option B: Pollard rho with distinguished points and batch inversion if memory is tight, though BSGS with array('I') fits in 64-128MB.
  See Hare's 48-bit entry: array('I') index/parity + array('H') fingerprint solved in 12.6s.
