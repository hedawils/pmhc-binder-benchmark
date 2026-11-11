from math import comb

def fisher(a, b, c, d):
    """Two-sided Fisher exact on [[a,b],[c,d]]."""
    n = a + b + c + d
    r1, r2, c1 = a + b, c + d, a + c
    def p(x):
        return comb(r1, x) * comb(r2, c1 - x) / comb(n, c1)
    obs = p(a)
    lo = max(0, c1 - r2)
    return sum(p(x) for x in range(lo, min(r1, c1) + 1) if p(x) <= obs * (1 + 1e-9))

def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    ph = k / n
    d = 1 + z*z/n
    c = (ph + z*z/(2*n)) / d
    h = z * ((ph*(1-ph)/n + z*z/(4*n*n)) ** 0.5) / d
    return (max(0, c - h) * 100, min(1, c + h) * 100)

tests = [
    ("interface pAE < 7 (marginal, per sequence)",      2, 60,  63, 449),
    ("three-filter survivors (sequential, per seq)",    2, 60,  37, 449),
    ("backbones with any hit",                          2, 10,  33, 141),
    ("binder pLDDT > 92 (marginal, per sequence)",     37, 60, 187, 449),
    ("scRMSD < 2.0 (marginal, per sequence)",          50, 60, 361, 449),
]

print("arm 1 = fold-conditioned on KH46's own fold, NY-ESO-1, n=60 seqs / 10 backbones")
print("arm 3 = Liu et al. privileged scaffolds, NY-ESO-1, n=449 seqs / 141 backbones (53%% scored)")
print()
print("%-46s %18s %18s %10s %s" % ("", "arm 1", "arm 3", "p", "verdict"))
for name, k1, n1, k3, n3 in tests:
    r1, r3 = 100*k1/n1, 100*k3/n3
    l1, u1 = wilson(k1, n1)
    l3, u3 = wilson(k3, n3)
    p = fisher(k1, n1 - k1, k3, n3 - k3)
    overlap = not (u1 < l3 or u3 < l1)
    v = "n.s." if p >= 0.05 else "DIFFERENT"
    print("%-46s %5.1f%% [%4.1f-%4.1f] %5.1f%% [%4.1f-%4.1f] %10.4f %s%s"
          % (name, r1, l1, u1, r3, l3, u3, p, v, "  (CIs overlap)" if overlap else ""))
print()
print("fold =", "%.2fx" % (( 63/449)/(2/60)), "on pAE marginal;",
            "%.2fx" % ((37/449)/(2/60)), "on three-filter survivors")
