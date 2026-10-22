# Escape from geometric lock-in in GCD-augmented Fibonacci recurrences — reproduction code

This repository reproduces every numerical result, both tables and Figure 1 of

> Hawarey, M. (2026). Escape from geometric lock-in in GCD-augmented Fibonacci recurrences. *AIR Journal of Mathematics and Computational Sciences*.

The recurrence is C(n) = C(n−1) + C(n−2) + gcd(C(n−1), C(n−2)). The script works with the states (α, β), the map R and the step multiplier g = gcd(β, α + 1). The function names and printed labels follow the paper (Theorem 3, Proposition 4, Theorem 6, Theorem 7, Table 1, Table 2, Figure 1).

## Requirements

- Python 3.10 or later, with `sympy` and `matplotlib`:
  `pip install sympy matplotlib`
- Optional, for the robustness pass and the m√5 family in Section 5: `pip install cypari2`

## Usage

```
python reproduce.py --quick    # Theorem 3 remark, Proposition 4, Theorems 6-7, Table 1, Figure 1 (a few minutes)
python reproduce.py --census   # Table 2 and all Section 5 numbers (box 200; several hours on 2 cores)
```

- **Checkpointing:** the census writes `census_box200.jsonl` one seed at a time. If it is interrupted, run the same command again and it resumes where it stopped.
- **Smaller box:** `--box B` runs the census for a smaller box B. Table 2 of the paper uses `--box 200`, which also contains box 120.

## Certificate

Theorem 3 decides whether a state is ever amplified again. Only primes dividing N² − 1 can amplify, where N = v² − uv − u² is the norm of the shifted state (u, v) = (α + 1, β + 1). For each such prime r, the orbit of (u, v) under the Fibonacci matrix modulo r is checked for the pair (0, 1).
- For r ≤ 3·10⁶ the orbit is iterated.
- For larger r ≡ ±1 (mod 5) up to 10¹², a discrete logarithm is used.
- Other cases are reported as undecided.

## Citation

See `CITATION.cff`. Please cite the paper and the archived version of this code.

## License

MIT
