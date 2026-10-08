#!/usr/bin/env python3
"""ECDLP via negation-map BSGS. Stdlib only."""
import json
import sys
from math import isqrt


def solve(data):
    p = int(data["p"])
    a = int(data["a"]) % p
    b = int(data["b"]) % p
    n = int(data["n"])
    Gx = int(data["Gx"]) % p
    Gy = int(data["Gy"]) % p
    Px = int(data["Px"]) % p
    Py = int(data["Py"]) % p

    def inv(x):
        return pow(x, -1, p)

    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + a) * inv(2 * y1) % p
        else:
            lam = (y2 - y1) * inv(x2 - x1) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def mul(k, P):
        if k < 0:
            k = -k
            if P is not None:
                P = (P[0], (-P[1]) % p)
        R = None
        while k:
            if k & 1:
                R = add(R, P)
            P = add(P, P)
            k >>= 1
        return R

    # Negation-map BSGS: stride M = 2m+1 covers j in [-m, m].
    m = isqrt(n) // 2 + 1
    M = 2 * m + 1
    G = (Gx, Gy)
    P = (Px, Py)

    baby = {}
    R = None
    for j in range(1, m + 1):
        R = add(R, G)
        if R is None:
            # order divides j; only possible if j == n, which m < n
            continue
        baby[R[0]] = (j << 1) | (R[1] & 1)

    MG = mul(M, G)
    # step by -M*G
    if MG is None:
        negMG = None
    else:
        negMG = (MG[0], (-MG[1]) % p)

    R = P
    limit = n // M + 2
    found = None
    for i in range(limit + 1):
        if R is None:
            found = (i * M) % n
            break
        rec = baby.get(R[0])
        if rec is not None:
            j = rec >> 1
            ybit = rec & 1
            if (R[1] & 1) == ybit:
                found = (i * M + j) % n
            else:
                found = (i * M - j) % n
            break
        R = add(R, negMG)

    if found is None:
        # fallback: plain BSGS
        m2 = isqrt(n) + 1
        baby2 = {None: 0}
        R = None
        for j in range(m2):
            baby2[R[0] if R is not None else None] = j
            # store full point to disambiguate
            if R is not None:
                baby2[(R[0], R[1])] = j
            R = add(R, G)
        step = mul(m2, G)
        if step is None:
            neg = None
        else:
            neg = (step[0], (-step[1]) % p)
        R = P
        for i in range(m2 + 2):
            key = None if R is None else (R[0], R[1])
            if key in baby2:
                found = (i * m2 + baby2[key]) % n
                break
            R = add(R, neg)

    if found is None:
        raise SystemExit("no solution")

    # verify
    Q = mul(found, G)
    if Q is None or Q[0] != Px or Q[1] != Py:
        # try the other sign interpretation already handled; brute nearby not needed
        raise SystemExit("verify failed")
    return found


def main():
    data = json.loads(sys.stdin.read())
    k = solve(data)
    print(f"k={k}")


if __name__ == "__main__":
    main()
