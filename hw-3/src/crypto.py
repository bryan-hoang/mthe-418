import functools
import math


def euler_phi(x: int) -> int:
    return sum(math.gcd(x, a) == 1 for a in range(1, x))


def not_coprime(x: int):
    return [1] + [a for a in range(1, 1 + x) if math.gcd(x, a) != 1]


def extended_gcd(a: int, b: int):
    """
    Extended Euclidean Algorithm. For any a, b, returns values gcd(a, b), x,
    and y such that a*x + b*y = gcd(a, b)
    """
    if a == 0:
        return b, 0, 1
    g, x, y = extended_gcd(b % a, a)
    return g, y - (b // a) * x, x


def extended_euclid(a: int, b: int):
    r0, a0, b0 = a, 1, 0
    r1, a1, b1 = b, 0, 1
    print("   q   |     r |     a |     b")
    print("       | %5d | %5d | %5d" % (r0, a0, b0))
    print("       | %5d | %5d | %5d" % (r1, a1, b1))
    while r1 != 1:
        q = r0 // r1
        r0, a0, b0, r1, a1, b1 = r1, a1, b1, r0 % r1, a0 - q * a1, b0 - q * b1
        print(" %5d | %5d | %5d | %5d" % (q, r1, a1, b1))


def crt(a1: int, m1: int, a2: int, m2: int):
    """Generalized Chinese Remainder Theorem (CRT).
    Find x, where x = a1 mod m1, x = a2 mod m2, if such a solution exists. It
    always exists when m1 and m2 are coprime, or if g = gcd(m1, m2), and a1 ==
    a2 mod g Returns an x, and the modulus (m1 * m2 for the coprime case)
    """
    g, x, y = extended_gcd(m1, m2)
    if a1 % g == a2 % g:
        ar, mr = (a2 * x * m1 + a1 * y * m2) // g, m1 * m2 // g
    else:
        raise ValueError(
            "The system x = %d mod %d, x = %d mod %d has no solution"
            % (a1, m1, a2, m2)
        )
    return ar % mr, mr


def dlp(g: int, h: int, p: int):
    """Return x s.t. g ^ x = h mod p"""
    return next(x for x in range(1, 1 + p) if pow(g, x, p) == h)


def shanks(a: int, b: int, n: int) -> int:
    # Return x s.t. a^x = b mod n
    m = math.ceil(math.sqrt(n))
    table = {pow(a, j, n): j for j in range(0, m)}
    ai = pow(a, -m, n)
    y = b
    for i in range(0, m):
        if y in table:
            return i * m + table[y]
        y = (y * ai) % n
    raise ValueError


def miller_rabin(n: int, a: int) -> bool:
    """Test n for primality, return False if composite, True if inconclusive"""
    if n % 2 == 0 or 1 < math.gcd(a, n) < n:
        return False
    q, k = n - 1, 0
    while q % 2 == 0:
        q, k = q // 2, k + 1
    a = pow(a, q, n)
    if a == 1:
        return True
    for i in range(k):
        if a == n - 1:
            return True
        a = pow(a, 2, n)
    return False


def pollard_p_minus_1(n: int, bound: int) -> int:
    """Uses the Pollard p-1 method to factor n"""
    a = 2
    for j in range(2, bound):
        a = pow(a, j, n)
        if 1 < (d := math.gcd(a - 1, n)) < n:
            return d
    return 0


@functools.lru_cache(None)
def divisors(n: int, proper: bool = False):
    if proper:
        return divisors(n, False) - {1, n}
    return set(
        functools.reduce(
            list.__add__,
            ([i, n // i] for i in range(1, math.isqrt(n) + 1) if n % i == 0),
        )
    )


@functools.lru_cache(None)
def prime_divisors(n: int):
    return {x for x in divisors(n) if prime(x)}


@functools.lru_cache(None)
def prime(n: int):
    return divisors(n) == {1, n}
