import math
from typing import List

import numpy as np

from crypto import (
    crt,
    divisors,
    euler_phi,
    extended_euclid,
    extended_gcd,
    miller_rabin,
    prime,
    prime_divisors,
    shanks,
)


def main():
    print("3a", shanks(11, 21, 71))
    print("3b", shanks(156, 116, 593))
    print("3c", shanks(650, 2213, 3571))

    print("4b", crt(137, 423, 87, 191))
    print("4d", crt(*crt(5, 9, 6, 10), 7, 11))

    print(
        "5c",
        divisors(1159),
        prime(19),
        prime(61),
        73 * 577 % 1080,
        pow(614, 577, 1159),
        pow(158, 73, 1159),
    )
    print(
        "5d",
        divisors(8023),
        prime(71),
        prime(113),
        7151 * 751 % 7840,
        pow(677, 7151, 8023),
        pow(1355, 751, 8023),
    )

    q7b("7b.i", 577, 60, 1463)
    q7b("7b.ii", 959, 1583, 1625)
    q7b("7b.iii", 133957, 224689, 2134440)

    q8("8b", [10988423, 25910155], [16784693, 11514115], 38749709)
    q8(
        "8c",
        [70583995, 173111957, 180311381],
        [4911157, 7346999, 29597249],
        225022969,
    )
    q8(
        "8d",
        [1103927639, 1022313977, 387632407],
        [76923209, 106791263, 7764043],
        1291233941,
    )

    q9b("9b.i", 1729)
    q9b("9b.i", 10585)
    q9b("9b.i", 75361)
    q9b("9b.i", 1024651)

    q9e("9e.i", 29341)
    q9e("9e.ii", 172947529)

    q10("10a", 1105)
    q10("10b", 294409)
    q10("10c", 294439)
    q10("10d", 118901509)
    q10("10e", 118901521)
    q10("10f", 118901527)
    q10("10g", 118915387)

    q13("13a", 1739)
    q13("13b", 220459)
    q13("13c", 48356747)

    print("14a")
    for i in range(2, 1 + 10):
        print("2^%d - 1 =" % i, p := 2**i - 1, prime_divisors(p))


def q7b(n: str, e: int, c: int, m: int):
    print("%s\nphi(N)" % n, phi := euler_phi(m))
    print("primes", prime_divisors(m))
    g, _, y = extended_gcd(phi, e)
    assert g == 1
    extended_euclid(phi, e)
    print(y, x := pow(c, y, m))
    assert pow(x, e, m) == c


def q8(text: str, e: List[int], d: List[int], N: int):
    print(text, *(pq := divisors(N, True)), end=", ")
    p, q = pq
    print("g =", g := math.gcd(p - 1, q - 1), end=", ")
    t = math.gcd(*[e0 * d0 - 1 for e0, d0 in zip(e, d)])
    print("T =", t, end=", ")
    g0 = (q - 1) * (p - 1) // g
    assert t % g0 == 0
    print("k =", k := t // g0, end=", ")
    print("n =", n := k * g, end=", ")
    print("roots =", pn := np.roots([1, n * t - N - 1, N]))
    assert set(map(round, pn)) == pq


def q9b(n: str, q: int):
    print(
        n,
        lcm := math.lcm(*[p - 1 for p in prime_divisors(q) - {1, q}]),
        (q - 1) // lcm,
    )
    assert (q - 1) % lcm == 0


def q9e(n: str, q: int):
    print(n, pd := prime_divisors(q) - {1, q}, end=" ")
    for p in pd:
        assert (q - 1) % (p - 1) == 0
        print(p - 1, "*", (q - 1) // (p - 1), "= p", end=" ")
    print("")


def q10(n: str, q: int):
    wi = [*range(10, 1 + 20)]
    ws = [miller_rabin(q, a) for a in wi]
    print(
        n,
        q,
        "probably prime" if (w := all(ws)) else "composite",
        wi if w else wi[ws.index(False)],
    )
    assert w == prime(q)


def q13(n: str, q: int):
    """Uses the Pollard p-1 method to factor n"""
    print(n, q)
    j = a = 2
    while True:
        a = pow(a, j, q)
        d = math.gcd(a - 1, q)
        print("  ", j, "a^j mod n =", a, "gcd(a-1, q) =", d)
        if 1 < d < q:
            print(q, "=", q // d, "*", d)
            print(prime_divisors(q))
            assert q % d == 0
            break
        j += 1


if __name__ == "__main__":
    main()
