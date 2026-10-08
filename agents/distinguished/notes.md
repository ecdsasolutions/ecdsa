36-bit solved: negation-map BSGS, m=isqrt(n//2)+1, stride 2m+1, baby[x]=(j<<1)|(y&1), precompute M*G then one affine add per giant step. egcd pow(x,-1,p) is ~3x faster than Fermat at these sizes (0.42s vs 1.28s). Verified k*G==P.

Next (38/40-bit): same approach; expect ~0.8s/1.6s. Consider Montgomery batch inversion for giant steps if >40 bits (book shows batch64 optimal at 40-bit). Memory ~2x babies per 2 bits, still fine in Python dict.
