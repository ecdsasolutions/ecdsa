import sys, json, math

def solve(inp):
    p = inp["p"]
    a = inp["a"]
    n = inp["n"]
    Gx, Gy = inp["Gx"], inp["Gy"]
    Px, Py = inp["Px"], inp["Py"]

    G = (Gx, Gy)
    P = (Px, Py)

    def inv(v):
        return pow(v, -1, p)

    def neg(pt):
        x, y = pt
        return (x, (-y) % p)

    def add(Q, R):
        if Q is None:
            return R
        if R is None:
            return Q
        x1, y1 = Q
        x2, y2 = R
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + a) * inv((2 * y1) % p) % p
        else:
            lam = ((y2 - y1) % p) * inv((x2 - x1) % p) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def scalarmul(k, pt):
        acc = None
        while k:
            if k & 1:
                acc = add(acc, pt)
            pt = add(pt, pt)
            k >>= 1
        return acc

    m = math.isqrt(n) // 2 + 1
    M = 2 * m + 1

    # baby steps: store x -> (j<<1)|(y&1) for j in [1, m]
    baby = {}
    pt = None
    for j in range(1, m + 1):
        pt = add(pt, G)
        x, y = pt
        baby[x] = (j << 1) | (y & 1)

    # giant steps: R = P - i*M*G
    MG = scalarmul(M, G)
    negMG = neg(MG)
    R = P
    limit = n // M + 2
    k = 0
    for i in range(limit):
        if R is None:
            k = (i * M) % n
            break
        x, y = R
        v = baby.get(x)
        if v is not None:
            j = v >> 1
            if (v & 1) == (y & 1):
                k = (i * M + j) % n
            else:
                k = (i * M - j) % n
            break
        R = add(R, negMG)
    return k

if __name__ == "__main__":
    inp = json.loads(sys.stdin.read())
    k = solve(inp)
    print("k=%d" % k)