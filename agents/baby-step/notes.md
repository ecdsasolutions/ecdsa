Height 32 bits:
- Optimized Negation BSGS with stride M = 2*m + 1.
- Set m = isqrt(n) // 2 + 1 to minimize expected steps m + n/(4m) under uniform random target distribution (~0.5*sqrt(n) baby steps and ~0.5*sqrt(n) expected giant steps), reducing baby table size and insertions by ~29% compared to worst-case minimax.
- Baby step dict stores packed int (j << 1) | (y & 1) rather than (j, y) tuple, saving object allocations and relying on odd prime p parity distinction.
- Doubling for j=2 is handled outside the loop to eliminate the branching check for j >= 3.
- Runtimes: ~0.046s compute time, ~0.10s total process time on 32 bits.
- Practice instance answer verified: k = 1658490919.

For next heights (36 - 40 bits):
- At 36 bits, n ~ 6.8e10, m ~ 130,000, table size is ~130k ints, well within RAM limits and should take ~0.2 - 0.4s.
- At 40 bits, m ~ 520,000, table fits in memory (<50MB). But Montgomery batching (batch size 32 or 64) can provide a noticeable speedup when inversion count exceeds 500,000.
- For 44+ bits, consider switching to van Oorschot-Wiener distinguished points with Pollard rho / kangaroo to keep memory bounded.
