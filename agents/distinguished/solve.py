import sys, json
from math import isqrt

def solve(p, a, b, Gx, Gy, n, Px, Py):
    # Negation-map baby-step giant-step (Bernstein-Lange / standard).
    # m = isqrt(n//2)+1, M = 2m+1.  k = i*M +/- j with j in [0, m].
    m = isqrt(n // 2) + 1
    M = 2 * m + 1

    def inv(x):
        return pow(x % p, -1, p)   # extended-gcd inverse, ~6x faster than Fermat

    def add(x1, y1, x2, y2):
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + a) * inv(2 * y1) % p
        else:
            lam = (y2 - y1) * inv(x2 - x1) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    # baby steps: j*G for j = 1..m, store x -> (j<<1)|(y&1)
    baby = {}
    x, y = Gx, Gy
    for j in range(1, m + 1):
        baby[x] = (j << 1) | (y & 1)
        x, y = add(x, y, Gx, Gy)

    # M*G via double-and-add, then negate for giant step -= M*G
    def mul(k):
        rx = ry = None
        tx, ty = Gx, Gy
        while k:
            if k & 1:
                if rx is None:
                    rx, ry = tx, ty
                else:
                    rr = add(rx, ry, tx, ty)
                    if rr is None:
                        rx = ry = None
                    else:
                        rx, ry = rr
            k >>= 1
            r = add(tx, ty, tx, ty)
            tx, ty = r if r is not None else (None, None)
        return rx, ry

    mx, my = mul(M)
    mx, my = mx, (-my) % p  # -M*G

    x, y = Px, Py
    for i in range(m + 1):
        if x is None:
            return (i * M) % n          # k = i*M exactly (j = 0)
        if x in baby:
            e = baby[x]
            j = e >> 1
            if (y & 1) == (e & 1):
                return (i * M + j) % n
            else:
                return (i * M - j) % n
        x, y = add(x, y, mx, my)
    return None

def main():
    d = json.load(sys.stdin)
    k = solve(int(d["p"]), int(d["a"]), int(d["b"]),
              int(d["Gx"]), int(d["Gy"]), int(d["n"]), int(d["Px"]), int(d["Py"]))
    print("k=%d" % k)

if __name__ == "__main__":
    main()