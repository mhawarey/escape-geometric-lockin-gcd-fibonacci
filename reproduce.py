#!/usr/bin/env python3
"""reproduce.py — reproduction code for
"Escape from geometric lock-in in GCD-augmented Fibonacci recurrences" (M. Hawarey, AIR Journal of Mathematics and
Computational Sciences). Version 1.0.0.

C(n) = C(n-1) + C(n-2) + gcd(C(n-1), C(n-2)). States (alpha, beta), map R, step multiplier g = gcd(beta, alpha+1);
shifted state (u, v) = (alpha+1, beta+1); norm N = v^2 - u v - u^2. Labels follow the paper.

  python reproduce.py --quick      Theorem 3 remark, Proposition 4, Theorem 6, Theorem 7, Table 1, Figure 1 (minutes)
  python reproduce.py --census     Table 2 and every Section 5 number (box 200 census: several hours on 2 cores;
                                   checkpointed to census_box200.jsonl, rerun to resume). Robustness pass and the
                                   m*sqrt(5) family need cypari2 (pip install cypari2).
Requires Python >= 3.10, sympy, matplotlib.
"""
import argparse, json, math, os, statistics as st, sys, time
from collections import Counter
from math import gcd
from multiprocessing import Pool
from sympy import factorint, sqrt_mod, n_order, discrete_log, primerange, multiplicity

RITER, RLOG = 3_000_000, 10**12        # orbit test: iteration up to RITER; discrete log for split primes up to RLOG
class Undecided(Exception): pass

def R(s):
    a, b = s; g = gcd(b, a + 1); return (b // g, (a + b + 1) // g), g
def norm(u, v): return v * v - u * v - u * u

# ---- Theorem 3: J_r = { j >= 0 : Q^j(u, v) = (0, 1) mod r } -------------------------------------------------------
def hit_iter(u, v, r):
    a, b = u % r, v % r; start = (a, b); j = 0; one = 1 % r
    while True:
        if a == 0 and b == one: return j
        a, b = b, (a + b) % r; j += 1
        if (a, b) == start: return None
def hit_split(u, v, r):                  # r = +-1 mod 5: z1 = rho^-j with (-1)^j = N mod r
    s5 = sqrt_mod(5, r); i2 = pow(2, -1, r); rho, rho2 = (1 + s5) * i2 % r, (1 - s5) * i2 % r
    z1, z2 = (u * rho + v - u) % r, (u * rho2 + v - u) % r
    if z1 == 0 or z2 == 0: return None
    try: j1 = discrete_log(r, pow(z1, -1, r), rho)
    except ValueError: return None
    o = n_order(rho, r); want = 0 if z1 * z2 % r == 1 else 1
    for j in (j1, j1 + o):
        if j % 2 == want: return j
    return None
def hit(u, v, r):
    if r <= RITER: return hit_iter(u, v, r)
    if r % 5 in (1, 4) and r <= RLOG: return hit_split(u, v, r)
    raise Undecided(r)
def certificate(u, v, digit_cap=30):
    """Least amplifying step from the shifted state (u, v), or None (never amplified again): Theorem 3 / Lemma 2(c)."""
    N0 = norm(u, v)
    if N0 * N0 == 1:
        a, b = u, v
        for j in range(3):
            if gcd(b - 1, a) > 1: return j, None
            a, b = b, a + b
    if len(str(abs(N0) + 1)) > digit_cap: raise Undecided("digits")
    P = sorted(set(factorint(abs(N0 - 1))) | set(factorint(abs(N0 + 1))))
    hs = [h for h in (hit(u, v, r) for r in P) if h is not None]
    return (min(hs) if hs else None), P
def escape_data(a, b, max_steps=1200, digit_cap=30):
    s, d, n, amps, check, seen = (a, b), 1, 2, [], True, set()
    for _ in range(max_steps):
        if s in seen: return "cycle"
        seen.add(s)
        if check:
            try:
                j, P = certificate(s[0] + 1, s[1] + 1, digit_cap)
                if j is None: return dict(n0=n, state=s, d=d, N=norm(s[0] + 1, s[1] + 1), primes=P, amps=amps)
            except Undecided: pass
            check = False
        s2, g = R(s)
        if g > 1: amps.append((n - 2, g)); check = True
        s, d, n = s2, d * g, n + 1
    return None
def pisano(r):
    a, b, k = 0, 1, 0
    while True:
        a, b = b, (a + b) % r; k += 1
        if (a, b) == (0, 1 % r): return k

def quick():
    F = [0, 1]; L = [2, 1]
    while len(L) < 3000: F.append(F[-1] + F[-2]); L.append(L[-1] + L[-2])
    # Theorem 3 remark: late first amplifications among states with entries <= 300
    late = []
    for a in range(1, 301):
        for b in range(1, 301):
            if gcd(a, b) == 1:
                j, P = certificate(a + 1, b + 1)
                if j is not None and j >= 3000: late.append((j, (a, b), [r for r in P if hit_iter(a + 1, b + 1, r) == j]))
    late.sort(); print(f"Theorem 3 remark: {len(late)} states (entries <= 300) first amplify after step 3,000; latest: {late[-1]}")
    full = [r for r in primerange(2, 10000) if r in (2, 3, 5) or (r % 5 in (2, 3) and pisano(r) == 2 * (r + 1))]
    print(f"Proposition 4: full primes below 100: {[r for r in full if r < 100]}; below 10^4: {len(full)}")
    # Theorem 6
    c = [None, 11, 5]
    while len(c) < 2502: c.append(c[-1] + c[-2] + gcd(c[-1], c[-2]))
    print("Theorem 6(a): C(n) = 6 L_(n-1) - 1 for n <= 2501:", all(c[n] == 6 * L[n - 1] - 1 for n in range(1, 2502)))
    for r in (179, 181):
        s5 = int(sqrt_mod(5, r)); s5 = min(s5, r - s5); rho = (1 + s5) * pow(2, -1, r) % r; z1 = 6 * s5 % r
        print(f"  r = {r}: sqrt5 = {s5}, rho = {rho} (order {n_order(rho, r)}), z1 = {z1} (order {n_order(z1, r)}), hit = {hit_iter(12, 6, r)}")
    c = [None, 1, 135]
    while len(c) < 2502: c.append(c[-1] + c[-2] + gcd(c[-1], c[-2]))
    print("Theorem 6(b): C(1..14) =", c[1:15], "; C(n) = 27(6 L_(n-3) - 1) for 11 <= n <= 2501:",
          all(c[n] == 27 * (6 * L[n - 3] - 1) for n in range(11, 2502)))
    print("Theorem 7: seeds (48894t - 1, 78961t), t <= 1000, map to (281, 455) with multiplier 281t:",
          all(gcd(48894 * t - 1, 78961 * t) == 1 and R((48894 * t - 1, 78961 * t)) == ((281, 455), 281 * t) for t in range(1, 1001)))
    print("Table 1 (p, v3(p), amplifications, n0, state, N0, primes):")
    for p in range(3, 400, 6):
        E = escape_data(1, p, 1200)
        if isinstance(E, dict): print(f"  {p}, {multiplicity(3, p)}, {len(E['amps'])}, {E['n0']}, {E['state']}, {E['N']}, {E['primes']}")
    figure1()

def figure1():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    def lb(a, b, T=2400):
        s = (a, b); out = []
        for _ in range(T + 1):
            x = str(s[1]); out.append(math.log10(int(x[:17])) + max(0, len(x) - 17)); s, _ = R(s)
        return out
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 10})
    fig, ax = plt.subplots(figsize=(6.3, 3.4)); n = list(range(2401))
    for (a, b), lab, col, lw in (((1, 45), "seed (1, 45): undecided", "#7a7a7a", 1.0), ((1, 135), "seed (1, 135): escape", "#1f4e9c", 1.4),
                                 ((1, 9), "seed (1, 9): 4-cycle (29, 47)", "#b5442c", 1.4)):
        ax.plot(n, lb(a, b), label=lab, color=col, lw=lw)
    ax.set_xlabel(r"step $n$"); ax.set_ylabel(r"$\log_{10}\beta_n$"); ax.set_xlim(0, 2400); ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left"); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig("figure1.png", dpi=400); print("Figure 1 written: figure1.png")

KNOWN = {(1, 2): "fixed point (1,2)", (3, 5): "2-cycle (3,5)", (29, 47): "4-cycle (29,47)", (97, 157): "period-22 (97,157)",
         (4003, 6477): "period-65 (4003,6477)"}
def run_seed(ab, T=2400):
    a, b = ab; s = (a, b); seen = {s: 0}; amps = 0; check = True; digits = {}; amps_at = {}
    for k in range(1, T + 1):
        if check:
            try:
                j, P = certificate(s[0] + 1, s[1] + 1, 30)
                if j is None: return dict(a=a, b=b, fate="escape", cert=dict(step=k - 1, state=list(s), N=norm(s[0] + 1, s[1] + 1)))
            except Undecided: pass
            check = False
        s, g = R(s)
        if g > 1: amps += 1; check = True
        if s in seen:
            cyc = [s]; t, _ = R(s)
            while t != s: cyc.append(t); t, _ = R(t)
            m = min(cyc); return dict(a=a, b=b, fate="cycle", repeat_step=k, cycle=KNOWN.get(m, f"period-{len(cyc)} {m}"))
        seen[s] = k
        if k in (1200, 2400): digits[k] = len(str(s[1])); amps_at[k] = amps
    return dict(a=a, b=b, fate="undecided", digits=digits, amps_at=amps_at, last_state=list(s))
def census(B=200, out="census_box200.jsonl"):
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try: r = json.loads(line); done.add((r["a"], r["b"]))
            except Exception: pass
    seeds = [(a, b) for a in range(1, B + 1) for b in range(1, B + 1) if gcd(a, b) == 1 and (a, b) not in done]
    print(f"census box {B}: {len(done)} done, {len(seeds)} to go"); t0 = time.time()
    with Pool() as pool, open(out, "a") as f:
        for i, r in enumerate(pool.imap_unordered(run_seed, seeds, chunksize=20)):
            f.write(json.dumps(r) + "\n"); f.flush()
            if (i + 1) % 2000 == 0: print(f"  {i + 1}/{len(seeds)} ({time.time() - t0:.0f} s)", flush=True)
    rows = [json.loads(l) for l in open(out)]
    order = ["fixed point (1,2)", "2-cycle (3,5)", "4-cycle (29,47)", "period-22 (97,157)", "period-65 (4003,6477)", "certified escape", "undecided"]
    for Bx in sorted({b for b in (50, 100, 120, 150, 200) if b <= B} | {B}):
        t = Counter()
        for r in rows:
            if r["a"] > Bx or r["b"] > Bx: continue
            if r["fate"] == "cycle" and r["repeat_step"] <= 1200: t[r["cycle"]] += 1
            elif r["fate"] == "escape" and r["cert"]["step"] <= 1200: t["certified escape"] += 1
            else: t["undecided"] += 1
        n = sum(t.values())
        print(f"Table 2 / Section 5, box {Bx} (n = {n}): " + "; ".join(f"{k}: {t[k]} ({100 * t[k] / n:.1f}%)" for k in order))
    und = [r for r in rows if r["fate"] == "undecided"]
    for T in ("2400",):
        print(f"undecided (box {B}), after {T} steps: digits min/median/max {min(r['digits'][T] for r in und)}/"
              f"{int(st.median(r['digits'][T] for r in und))}/{max(r['digits'][T] for r in und)}; amplifications median "
              f"{int(st.median(r['amps_at'][T] for r in und))}")
    print("every certified escape has 30 | N (Proposition 4):", all(r["cert"]["N"] % 30 == 0 for r in rows if r["fate"] == "escape"))
    esc = [r for r in rows if r["fate"] == "escape"]
    def canon(u, v):
        while v > u > 0: u, v = v - u, u
        return (u, v)
    cc = Counter(canon(r["cert"]["state"][0] + 1, r["cert"]["state"][1] + 1) for r in esc)
    print(f"escaping orbits (box {B}, 30-digit certificates): {len(esc)} escapes on {len(cc)} orbits; orbit of Theorem 6(a): {cc[(12, 6)]}")
    s1 = [run_seed((1, p)) for p in range(1, 2001)]; e1 = [r["b"] for r in s1 if r["fate"] == "escape"]
    print(f"seeds (1, p), p <= 2000: {len(e1)} escape; counts at 250/500/1000/1500/2000: {[sum(1 for p in e1 if p <= P) for P in (250, 500, 1000, 1500, 2000)]}")
    try:
        import cypari2
    except ImportError:
        print("cypari2 not installed: robustness pass and m*sqrt(5) family skipped"); return
    pari = cypari2.Pari(); pari.allocatemem(1 << 30, silent=True)
    def hit_any(u, v, r, N0):          # membership of z in <phi> mod r via element orders
        if r <= RITER: return hit_iter(u, v, r) is not None
        if r % 5 in (1, 4):
            s5 = pari(f"sqrt(Mod(5,{r}))"); rho = (1 + s5) / 2; z1 = u * rho + (v - u)
            if z1 == 0: return False
            e = int(pari.znorder(rho))
            if z1 ** e != 1: return False
            if e % 2: return True
            return (0 if z1 ** (e // 2) == 1 else 1) == (0 if N0 % r == 1 else 1)
        t = pari(f"ffgen(Mod(1,{r})*(x^2-x-1),'t)"); return (u * t + (v - u)) ** int(pari.fforder(t)) == 1
    def primes_of(n, secs=20):
        try: return [int(p) for p in pari(f"alarm({secs}, factor({abs(n)})[,1]~)")]
        except BaseException: return None
    res = Counter()
    for r in und:                      # robustness pass: state after the last amplification within 2,400 steps
        s = (r["a"], r["b"]); last = s
        for _ in range(2400):
            s, g = R(s)
            if g > 1: last = s
        u, v = last[0] + 1, last[1] + 1; N0 = norm(u, v)
        if N0 * N0 == 1 or len(str(abs(N0))) > 100: res["unresolved"] += 1; continue
        P1, P2 = primes_of(N0 - 1), primes_of(N0 + 1)
        if P1 is None or P2 is None: res["unresolved"] += 1; continue
        res["certified escape" if not any(hit_any(u, v, p, N0) for p in set(P1) | set(P2)) else "amplifies again"] += 1
    print(f"robustness pass (box {B} undecided seeds):", dict(res))
    ms = []
    for m in range(2, 3001):
        u, v = 2 * m, m; N0 = -5 * m * m
        P = set(int(p) for p in pari.factor(abs(N0 - 1))[0]) | set(int(p) for p in pari.factor(abs(N0 + 1))[0])
        if not any(hit_any(u, v, r, N0) for r in P): ms.append(m)
    print(f"m*sqrt(5) family: {len(ms)} values of m <= 3000 never amplified; first {ms[:5]}; counts at 500/1000/2000/3000: "
          f"{[sum(1 for m in ms if m <= M) for M in (500, 1000, 2000, 3000)]}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--census", action="store_true")
    ap.add_argument("--box", type=int, default=200)
    a = ap.parse_args()
    if not (a.quick or a.census): ap.print_help(); sys.exit()
    if a.quick: quick()
    if a.census: census(a.box, f"census_box{a.box}.jsonl")
