import sys, json
from math import isqrt

def batch_inv(vals, p):
    n = len(vals)
    res = [0]*n
    c = [0]*n
    acc = 1
    for i in range(n):
        c[i] = acc
        acc = acc * vals[i] % p
    invacc = pow(acc, -1, p)
    for i in range(n-1, -1, -1):
        res[i] = invacc * c[i] % p
        invacc = invacc * vals[i] % p
    return res

def solve(data):
    p = data['p']; a = data['a']; b = data['b']
    Gx = data['Gx']; Gy = data['Gy']
    n = data['n']; Px = data['Px']; Py = data['Py']

    m = isqrt(n // 2) + 1
    M = 2*m + 1

    # add affine point (x2,y2) [Z=1] to Jacobian (X1,Y1,Z1); returns Jacobian
    def add_affine(X1, Y1, Z1, x2, y2):
        if Z1 == 0:
            return x2, y2, 1
        Z1Z1 = Z1 * Z1 % p
        U2 = x2 * Z1Z1 % p
        S2 = y2 * Z1 % p * Z1Z1 % p
        U1 = X1; S1 = Y1
        H = (U2 - U1) % p
        R = (S2 - S1) % p
        if H == 0:
            if R == 0:   # P == Q: double
                A = X1 * X1 % p
                B = Y1 * Y1 % p
                C = B * B % p
                D = (2*((X1+B)*(X1+B) - A - C)) % p
                ZZ = Z1*Z1 % p
                E = (3*A + a*ZZ*ZZ) % p
                F = E*E % p
                X3 = (F - 2*D) % p
                Y3 = (E*(D - X3) - 8*C) % p
                Z3 = (2*Y1*Z1) % p
                return X3, Y3, Z3
            else:         # P == -Q: infinity
                return 0, 0, 0
        H2 = H*H % p
        H3 = H*H2 % p
        U1H2 = U1*H2 % p
        X3 = (R*R - H3 - 2*U1H2) % p
        Y3 = (R*(U1H2 - X3) - S1*H3) % p
        Z3 = (Z1*H) % p
        return X3, Y3, Z3

    # ---- baby steps: j*G, j = 1..m (Jacobian, no inversions) ----
    Xs = [0]*m; Ys = [0]*m; Zs = [0]*m
    X, Y, Z = Gx, Gy, 1    # j = 1
    for j in range(m):
        Xs[j] = X; Ys[j] = Y; Zs[j] = Z
        X, Y, Z = add_affine(X, Y, Z, Gx, Gy)

    # ---- S = M*G (affine, double-and-add) ----
    def affine_mul(k, X, Y):
        rx, ry = 0, 0
        while k:
            if k & 1:
                if rx == 0:
                    rx, ry = X, Y
                else:
                    lam = (Y - ry) * pow(X - rx, -1, p) % p
                    xr = (lam*lam - rx - X) % p
                    yr = (lam*(rx - xr) - ry) % p
                    rx, ry = xr, yr
            k >>= 1
            lam = (3*X*X + a) * pow(2*Y, -1, p) % p
            xr = (lam*lam - 2*X) % p
            yr = (lam*(X - xr) - Y) % p
            X, Y = xr, yr
        return rx, ry

    Sx, Sy = affine_mul(M, Gx, Gy)
    nSy = (p - Sy) % p   # -S

    # ---- giant steps: Q_i = P - i*M*G, i = 0..m (Jacobian) ----
    Xg = [0]*(m+1); Yg = [0]*(m+1); Zg = [0]*(m+1)
    X, Y, Z = Px, Py, 1
    for i in range(m+1):
        if Z == 0:
            return (i * M) % n
        Xg[i] = X; Yg[i] = Y; Zg[i] = Z
        X, Y, Z = add_affine(X, Y, Z, Sx, nSy)

    # ---- batch invert all Z's ----
    allZ = Zs + Zg
    invZ = batch_inv(allZ, p)
    invZ_s = invZ[:m]
    invZ_g = invZ[m:]

    # ---- build baby dict: x -> (j<<1)|parity ----
    baby = {}
    for j in range(m):
        zinv = invZ_s[j]
        zinv2 = zinv*zinv % p
        x = Xs[j] * zinv2 % p
        y = Ys[j] * zinv2 % p * zinv % p
        if x not in baby:
            baby[x] = ((j+1) << 1) | (y & 1)

    # ---- scan giant steps ----
    for i in range(m+1):
        zinv = invZ_g[i]
        zinv2 = zinv*zinv % p
        x = Xg[i] * zinv2 % p
        y = Yg[i] * zinv2 % p * zinv % p
        e = baby.get(x)
        if e is not None:
            j = e >> 1
            par = e & 1
            if (y & 1) == par:
                return (i * M + j) % n
            else:
                return (i * M - j) % n
    return None

def main():
    data = json.load(sys.stdin)
    k = solve(data)
    print(f"k={k}")

if __name__ == "__main__":
    main()