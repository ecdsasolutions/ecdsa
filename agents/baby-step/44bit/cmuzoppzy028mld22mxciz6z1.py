import sys, json, math, time

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

    if Px is None:
        print("k=0")
        return

    if Px == Gx:
        if Py == Gy:
            print("k=1")
            return
        elif (Py + Gy) % p == 0:
            print(f"k={n - 1}")
            return

    def double(x, y):
        inv = pow(2 * y, -1, p)
        l = (3 * x * x + a) * inv % p
        nx = (l * l - 2 * x) % p
        ny = (l * (x - nx) - y) % p
        return nx, ny

    m = math.isqrt(n) // 2 + 1
    M = 2 * m + 1

    baby = {}
    baby[Gx] = (1 << 1) | (Gy & 1)

    cur_x, cur_y = double(Gx, Gy)
    baby[cur_x] = (2 << 1) | (cur_y & 1)

    for j in range(3, m + 1):
        inv = pow(cur_x - Gx, -1, p)
        lam = ((cur_y - Gy) * inv) % p
        cur_x = (lam * lam - cur_x - Gx) % p
        cur_y = (lam * (Gx - cur_x) - Gy) % p
        baby[cur_x] = (j << 1) | (cur_y & 1)

    def point_add(x1, y1, x2, y2):
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            return double(x1, y1)
        inv = pow(x2 - x1, -1, p)
        lam = ((y2 - y1) * inv) % p
        nx = (lam * lam - x1 - x2) % p
        ny = (lam * (x1 - nx) - y1) % p
        return nx, ny

    def scalar_mul(k_val, bx, by):
        rx, ry = None, None
        qx, qy = bx, by
        while k_val > 0:
            if k_val & 1:
                if rx is None:
                    rx, ry = qx, qy
                else:
                    rx, ry = point_add(rx, ry, qx, qy)
            qx, qy = double(qx, qy)
            k_val >>= 1
        return rx, ry

    MGx, MGy = scalar_mul(M, Gx, Gy)
    neg_MGx, neg_MGy = MGx, (-MGy) % p

    gamma_x, gamma_y = Px, Py
    max_i = n // M + 2
    k_found = None

    for i in range(max_i + 1):
        if gamma_x is None:
            k_found = (i * M) % n
            break

        if gamma_x in baby:
            val = baby[gamma_x]
            j = val >> 1
            if (gamma_y & 1) == (val & 1):
                k_found = (i * M + j) % n
            else:
                k_found = (i * M - j) % n
            break

        try:
            inv = pow(gamma_x - neg_MGx, -1, p)
            lam = ((gamma_y - neg_MGy) * inv) % p
            gamma_x = (lam * lam - gamma_x - neg_MGx) % p
            gamma_y = (lam * (neg_MGx - gamma_x) - neg_MGy) % p
        except ValueError:
            if (gamma_y + neg_MGy) % p == 0:
                gamma_x, gamma_y = None, None
            else:
                gamma_x, gamma_y = double(gamma_x, gamma_y)

    check_x, check_y = scalar_mul(k_found, Gx, Gy)
    assert (check_x, check_y) == (Px, Py), "Verification failed!"
    print(f"k={k_found}")

if __name__ == "__main__":
    solve()
