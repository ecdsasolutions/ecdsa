import sys, json, math

def inv(x, p):
    return pow(x, -1, p)

def add(P, Q, a, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a) * inv(2 * y1, p) % p
    else:
        lam = (y2 - y1) * inv(x2 - x1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)

def mul(k, P, a, p):
    R = None
    while k:
        if k & 1:
            R = add(R, P, a, p)
        P = add(P, P, a, p)
        k >>= 1
    return R

def solve():
    d = json.load(sys.stdin)
    p = d["p"]; a = d["a"]; n = d["n"]
    G = (d["Gx"], d["Gy"])
    P = (d["Px"], d["Py"])

    m = math.isqrt(n // 2) + 1
    M = 2 * m + 1

    # baby steps j*G for j=1..m, key by x with parity
    baby = {}
    cur = G
    for j in range(1, m + 1):
        x, y = cur
        baby[x] = (j << 1) | (y & 1)
        if j < m:
            cur = add(cur, G, a, p)

    # M*G precomputed
    MG = mul(M, G, a, p)
    negMG = (MG[0], (-MG[1]) % p)

    # giant steps R = P - i*M*G
    R = P
    i = 0
    while True:
        x, y = R
        e = baby.get(x)
        if e is not None:
            j = e >> 1
            ypar = e & 1
            if (y & 1) == ypar:
                k = (i * M + j) % n
            else:
                k = (i * M - j) % n
            if k != 0:
                print(f"k={k}", flush=True)
                return
        R = add(R, negMG, a, p)
        i += 1

if __name__ == "__main__":
    solve()