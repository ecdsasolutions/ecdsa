# Hare: 48-bit session

- Height: 48 bits
- Model: x-ai/grok-4.7
- Outcome: limit_time (time limit reached)
- Cost: $0.2645
- Turns: 11
- Duration: 304.2s
- Started: 2026-10-08 15:27:05 UTC
- Session: cmuzowkuq02cnld22gqoyc6gu

## Notes

(no notes saved this session)

## Actions

1. book_search '48-bit packed negation BSGS': 1 entries
2. book_search '48 bit BSGS batch inversion speedup': 0 entries
3. book_search 'packed uint64 open-addressed negation BSGS Montgomery': 0 entries
4. run solve.py on practice: exit 1 in 6.10s
5. write_file dump.py (1 lines)
6. run dump.py on practice: exit 0 in 0.03s
7. write_file dump.py (3 lines)
8. run dump.py on practice: exit 0 in 0.02s
9. write_file dump.py (4 lines)
10. run dump.py on practice: exit 0 in 0.03s
11. browser_search 'Pollard rho elliptic curve distinguished points Bernstein Lange negation map BSGS': 4 results
12. write_file solve.py (285 lines)
13. browser_open https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm_for_logarithms: Pollard's rho algorithm for logarithms - Wikipedia
14. run solve.py on practice: exit 0 in 18.21s
15. write_file verify.py (38 lines)

## Final solve.py

```python
#!/usr/bin/env python3
"""48-bit ECDLP: negation-map BSGS with packed uint64 open addressing.

k = i*M ± j, M = 2m+1, m = isqrt(n)//2+1.
Baby x-coordinates live in an open-addressed array('Q'):
  value = (fp << 25) | ((j << 1) | y_parity), 0 = empty.
fp = low64(x * golden) >> 25. Fingerprint hits are checked with a scalar mul.
"""
import json
import sys
from array import array

GOLD = 0x9E3779B97F4A7C15
MASK64 = (1 << 64) - 1
LOW = (1 << 25) - 1
BATCH = 128


def isqrt(n):
    x = int(n ** 0.5)
    while x * x > n:
        x -= 1
    while (x + 1) * (x + 1) <= n:
        x += 1
    return x


def solve(p, a, Gx, Gy, n, Px, Py):
    if Px == 0 and Py == 0:
        return 0

    Gx %= p
    Gy %= p
    Px %= p
    Py %= p

    m = isqrt(n) // 2 + 1
    M = 2 * m + 1

    slots = 1
    need = m + (m >> 1) + 8
    while slots < need:
        slots <<= 1
    mask = slots - 1
    table = array("Q", [0]) * slots

    def inv_batch(zs, L):
        pref = [0] * L
        acc = zs[0]
        pref[0] = acc
        for i in range(1, L):
            acc = (acc * zs[i]) % p
            pref[i] = acc
        iv = pow(acc, -1, p)
        out = [0] * L
        for i in range(L - 1, 0, -1):
            out[i] = (iv * pref[i - 1]) % p
            iv = (iv * zs[i]) % p
        out[0] = iv
        return out

    def mul(k):
        if k % n == 0:
            return 0, 0
        k %= n
        rx = ry = None
        ax, ay = Gx, Gy
        while k:
            if k & 1:
                if rx is None:
                    rx, ry = ax, ay
                else:
                    dx = (ax - rx) % p
                    if dx == 0:
                        if (ay + ry) % p == 0:
                            rx = ry = None
                        else:
                            lam = ((3 * ax * ax + a) * pow((2 * ay) % p, -1, p)) % p
                            nx = (lam * lam - 2 * ax) % p
                            ny = (lam * (ax - nx) - ay) % p
                            rx, ry = nx, ny
                    else:
                        lam = ((ay - ry) * pow(dx, -1, p)) % p
                        nx = (lam * lam - rx - ax) % p
                        ny = (lam * (rx - nx) - ry) % p
                        rx, ry = nx, ny
            lam = ((3 * ax * ax + a) * pow((2 * ay) % p, -1, p)) % p
            nx = (lam * lam - 2 * ax) % p
            ny = (lam * (ax - nx) - ay) % p
            ax, ay = nx, ny
            k >>= 1
        if rx is None:
            return 0, 0
        return rx, ry

    def consider(x, y, i_base):
        mix = (x * GOLD) & MASK64
        h = mix & mask
        fp = mix >> 25
        tab = table
        while True:
            v = tab[h]
            if v == 0:
                return None
            if (v >> 25) == fp:
                low = v & LOW
                jv = low >> 1
                if (y & 1) == (low & 1):
                    cand = (i_base * M + jv) % n
                else:
                    cand = (i_base * M - jv) % n
                rx, ry = mul(cand)
                if rx == Px and ry == Py:
                    return cand
            h = (h + 1) & mask

    # S[i] = i*G for i = 1..BATCH
    Sx = [0] * (BATCH + 1)
    Sy = [0] * (BATCH + 1)
    Sx[1], Sy[1] = Gx, Gy
    lam = ((3 * Gx * Gx + a) * pow((2 * Gy) % p, -1, p)) % p
    nx = (lam * lam - 2 * Gx) % p
    ny = (lam * (Gx - nx) - Gy) % p
    Sx[2], Sy[2] = nx, ny
    for i in range(3, BATCH + 1):
        dx = (Gx - Sx[i - 1]) % p
        lam = ((Gy - Sy[i - 1]) * pow(dx, -1, p)) % p
        nx = (lam * lam - Sx[i - 1] - Gx) % p
        ny = (lam * (Sx[i - 1] - nx) - Sy[i - 1]) % p
        Sx[i], Sy[i] = nx, ny

    tab = table
    g = GOLD
    m64 = MASK64
    msk = mask

    def insert(x, j, yp):
        mix = (x * g) & m64
        h = mix & msk
        v = ((mix >> 25) << 25) | ((j << 1) | yp)
        while tab[h]:
            h = (h + 1) & msk
        tab[h] = v

    limit = m if m < BATCH else BATCH
    for i in range(1, limit + 1):
        insert(Sx[i], i, Sy[i] & 1)

    rx, ry = Sx[limit], Sy[limit]
    j = limit
    # reuse buffers
    dxs = [1] * BATCH
    while j < m:
        cnt = BATCH if j + BATCH <= m else m - j
        dbl = 0
        rj = rx
        for t in range(cnt):
            dx = (Sx[t + 1] - rj) % p
            if dx == 0:
                dbl |= 1 << t
                dxs[t] = 1
            else:
                dxs[t] = dx
        ivs = inv_batch(dxs, cnt)
        base = j
        for t in range(cnt):
            i = t + 1
            if (dbl >> t) & 1:
                lam = ((3 * rx * rx + a) * pow((2 * ry) % p, -1, p)) % p
                nx = (lam * lam - 2 * rx) % p
                ny = (lam * (rx - nx) - ry) % p
            else:
                lam = ((Sy[i] - ry) * ivs[t]) % p
                nx = (lam * lam - rx - Sx[i]) % p
                ny = (lam * (rx - nx) - ry) % p
            insert(nx, base + i, ny & 1)
            if i == cnt:
                rx, ry = nx, ny
        j += cnt

    # D = M*G, U[i] = -i*D
    rx_, ry_ = None, None
    ax, ay = Gx, Gy
    kk = M
    while kk:
        if kk & 1:
            if rx_ is None:
                rx_, ry_ = ax, ay
            else:
                dx = (ax - rx_) % p
                lam = ((ay - ry_) * pow(dx, -1, p)) % p
                nx = (lam * lam - rx_ - ax) % p
                ny = (lam * (rx_ - nx) - ry_) % p
                rx_, ry_ = nx, ny
        lam = ((3 * ax * ax + a) * pow((2 * ay) % p, -1, p)) % p
        nx = (lam * lam - 2 * ax) % p
        ny = (lam * (ax - nx) - ay) % p
        ax, ay = nx, ny
        kk >>= 1

    Ux = [0] * (BATCH + 1)
    Uy = [0] * (BATCH + 1)
    Ux[1], Uy[1] = rx_, (-ry_) % p
    lam = ((3 * Ux[1] * Ux[1] + a) * pow((2 * Uy[1]) % p, -1, p)) % p
    nx = (lam * lam - 2 * Ux[1]) % p
    ny = (lam * (Ux[1] - nx) - Uy[1]) % p
    Ux[2], Uy[2] = nx, ny
    ux1, uy1 = Ux[1], Uy[1]
    for i in range(3, BATCH + 1):
        dx = (ux1 - Ux[i - 1]) % p
        lam = ((uy1 - Uy[i - 1]) * pow(dx, -1, p)) % p
        nx = (lam * lam - Ux[i - 1] - ux1) % p
        ny = (lam * (Ux[i - 1] - nx) - Uy[i - 1]) % p
        Ux[i], Uy[i] = nx, ny

    qx, qy = Px, Py
    max_i = (n - 1 + m) // M
    i_base = 0
    while i_base <= max_i:
        hit = consider(qx, qy, i_base)
        if hit is not None:
            return hit
        remain = max_i - i_base
        if remain <= 0:
            break
        cnt = BATCH if remain > BATCH else remain
        infm = 0
        dbl = 0
        for t in range(cnt):
            dx = (Ux[t + 1] - qx) % p
            if dx == 0:
                if (qy + Uy[t + 1]) % p == 0:
                    infm |= 1 << t
                    dxs[t] = 1
                else:
                    dbl |= 1 << t
                    dxs[t] = 1
            else:
                dxs[t] = dx
        ivs = inv_batch(dxs, cnt)
        for t in range(cnt):
            i = t + 1
            idx = i_base + i
            if (infm >> t) & 1:
                cand = (idx * M) % n
                rx, ry = mul(cand)
                if rx == Px and ry == Py:
                    return cand
                continue
            ux = Ux[i]
            uy = Uy[i]
            if (dbl >> t) & 1:
                lam = ((3 * qx * qx + a) * pow((2 * qy) % p, -1, p)) % p
                nx = (lam * lam - 2 * qx) % p
                ny = (lam * (qx - nx) - qy) % p
            else:
                lam = ((uy - qy) * ivs[t]) % p
                nx = (lam * lam - qx - ux) % p
                ny = (lam * (qx - nx) - qy) % p
            hit = consider(nx, ny, idx)
            if hit is not None:
                return hit
            if i == cnt:
                qx, qy = nx, ny
        i_base += cnt
    return 0


def main():
    data = json.loads(sys.stdin.read())
    k = solve(
        int(data["p"]),
        int(data["a"]),
        int(data["Gx"]),
        int(data["Gy"]),
        int(data["n"]),
        int(data["Px"]),
        int(data["Py"]),
    )
    print(f"k={k}")


if __name__ == "__main__":
    main()
```
