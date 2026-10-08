40-bit solved: negation BSGS. m=isqrt(n//2)+1, stride 2m+1 = 569127.
k = i*M +/- j, j in [1,m], j=0 via infinity.

Two working impls:
- Plain affine (baby[x]=(j<<1)|(y&1), pow(x,-1,p) per step): ~1.0s, simpler.
- Jacobian + Montgomery batch_inv over all Z: ~1.3s, deterministic, handles infinity (Z==0 => k=i*M mod n). Not faster at 40 bits.

Both verified k*G==P. At 40 bits pow(x,-1,p) (egcd) is cheap; Jacobian not worth it yet.

Next (42/44-bit): memory is main concern. m doubles per 2 bits: 42-bit m~570k, 44-bit m~1.1M. Python dict still OK but use packed arrays (array('Q') open-addressing) per book for 44+ bits. Book shows batch64 inversion helps at 40+; sequential BSGS steps can't be batched (sequential dependency), batching only helps Pollard rho multi-walker.

Reference rho at 40 bits ≈ sqrt(n)≈402k ops in pure python (10s+), so BSGS is ~10x+ faster.
