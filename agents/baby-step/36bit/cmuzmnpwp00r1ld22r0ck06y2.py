import json
import sys
import math

def solve():
    data = json.load(sys.stdin)
    p = int(data["p"])
    a = int(data["a"])
    b = int(data["b"])
    Gx = int(data["Gx"])
    Gy = int(data["Gy"])
    n = int(data["n"])
    Px = int(data["Px"])
    Py = int(data["Py"])

    if Px == Gx and Py == Gy:
        print("k=1")
        return
    if Px == Gx and (Py + Gy) % p == 0:
        print(f"k={n-1}")
        return

    # Negation BSGS with stride M = 2*m + 1.
    # Expected operations minimized at m = isqrt(n) // 2 + 1
    m = math.isqrt(n) // 2 + 1
    M = 2 * m + 1

    # Invert helper using pow(d, -1, p)
    # Baby steps: store j*G for j = 1 .. m
    # Dict maps x -> (j << 1) | (y & 1)
    # Since p is odd prime, (y & 1) != ((-y % p) & 1).
    baby = {}
    baby[Gx] = (1 << 1) | (Gy & 1)

    # j = 2: double G
    inv = pow(2 * Gy, -1, p)
    lam = (3 * Gx * Gx + a) * inv % p
    x2 = (lam * lam - 2 * Gx) % p
    y2 = (lam * (Gx - x2) - Gy) % p
    baby[x2] = (2 << 1) | (y2 & 1)

    cur_x, cur_y = x2, y2
    for j in range(3, m + 1):
        inv = pow(cur_x - Gx, -1, p)
        lam = (cur_y - Gy) * inv % p
        nx = (lam * lam - cur_x - Gx) % p
        ny = (lam * (Gx - nx) - Gy) % p
        baby[nx] = (j << 1) | (ny & 1)
        cur_x, cur_y = nx, ny

    # Compute M * G
    # We have m*G = (cur_x, cur_y).
    # M*G = m*G + (m+1)*G = 2*(m*G) + G
    # Let's compute (m+1)*G first:
    inv = pow(cur_x - Gx, -1, p)
    lam = (cur_y - Gy) * inv % p
    x_m1 = (lam * lam - cur_x - Gx) % p
    y_m1 = (lam * (Gx - x_m1) - Gy) % p

    # Now M*G = m*G + (m+1)*G
    inv = pow(x_m1 - cur_x, -1, p)
    lam = (y_m1 - cur_y) * inv % p
    MGx = (lam * lam - cur_x - x_m1) % p
    MGy = (lam * (cur_x - MGx) - cur_y) % p

    # Step in giant steps: S = -M*G = (MGx, -MGy % p)
    Sx = MGx
    Sy = (-MGy) % p

    # Giant steps: Q_i = P + i * S = P - i * (M*G)
    # For i = 0: Q_0 = P
    if Px in baby:
        val = baby[Px]
        j = val >> 1
        par = val & 1
        if (Py & 1) == par:
            print(f"k={j}")
            return
        else:
            print(f"k={n - j}")
            return

    # Check if P == identity? Not possible since order n is prime and P != identity.
    # Now loop i >= 1
    # Q_1 = P + S
    if Px == Sx:
        # P == -MG or P == MG?
        pass

    inv = pow(Px - Sx, -1, p)
    lam = (Py - Sy) * inv % p
    Qx = (lam * lam - Px - Sx) % p
    Qy = (lam * (Px - Qx) - Py) % p

    i = 1
    while True:
        if Qx in baby:
            val = baby[Qx]
            j = val >> 1
            par = val & 1
            if (Qy & 1) == par:
                k = (i * M + j) % n
            else:
                k = (i * M - j) % n
            print(f"k={k}")
            return
        
        # Advance Q = Q + S
        inv = pow(Qx - Sx, -1, p)
        lam = (Qy - Sy) * inv % p
        Qx = (lam * lam - Qx - Sx) % p
        Qy = (lam * (Sx - Qx) - Sy) % p
        i += 1

if __name__ == "__main__":
    solve()
