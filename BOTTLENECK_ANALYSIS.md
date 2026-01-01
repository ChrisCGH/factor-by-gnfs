# Performance Analysis Results

## Key Finding: check_for_remaining_relations is the Bottleneck

Analysis of the LatticeSiever timing data reveals that **check_for_remaining_relations() takes 49% of total sieving time** - this is by far the dominant bottleneck.

### What This Function Does

The function iterates through potentially smooth points and calls `is_actually_smooth()` which:

1. **Primality testing** - Uses probabilistic primality tests on remaining quotients
2. **Full factorization** - Calls `factorise_no_trial()` to completely factor composite numbers
3. **Smooth verification** - Checks if all factors are within the smoothness bounds

### Why It's Expensive

From the code analysis in `LatticeSiever.cpp`:

```cpp
bool LatticeSiever::is_actually_smooth(FastVeryLong& remaining_quotient, ...) {
    // ...
    bool is_prime = remaining_quotient.is_probable_prime();  // EXPENSIVE
    // ...
    if (!remaining_quotient.factorise_no_trial(&fac))        // VERY EXPENSIVE
        return false;
    // ...
}
```

**Primality testing** and **factorization** are computationally intensive operations that can't be optimized by the compiler. They require:
- Modular exponentiation for primality tests
- Trial division or more sophisticated factorization algorithms
- Multiple iterations for probabilistic primality

### Time Breakdown from Profiling

```
Phase                          Time    Percentage
----------------------------------------------------
check relations               2.37s    49.0%  ← BOTTLENECK
sieve by vectors 1            0.90s    18.7%
sieve by vectors 2            0.51s    10.6%
eliminate rational            0.34s     7.1%
check interval 1              0.33s     6.9%
check interval 2              0.20s     4.2%
eliminate algebraic           0.09s     1.8%
remove factors (all)          0.08s     1.6%
```

## Optimization Opportunities

### 1. Early Exit Strategies (Easy, High Impact)

Before calling expensive `is_actually_smooth()`, add quick rejection tests:

```cpp
// Check if remaining_quotient is obviously too large
if (remaining_quotient > L * L) continue;

// Check if remaining_quotient has small factors we already know about
// (could use a Bloom filter or hash table)
```

### 2. Batch Primality Testing (Medium Difficulty)

Instead of testing one number at a time, accumulate candidates and test them in batches:
- GMP has optimized batch operations
- Reduces function call overhead
- Better CPU cache utilization

### 3. Faster Factorization (Medium-Hard)

The `factorise_no_trial()` function could be optimized:
- Use Pollard's rho algorithm for small composites
- ECM for slightly larger numbers
- Better trial division with wheel factorization

### 4. Parallel Processing (Hard, Very High Impact)

The loop in `check_for_remaining_relations()` is embarrassingly parallel:

```cpp
for (PotentiallySmoothPoint* smooth_iter = head_psp_; smooth_iter; ...) {
    // Each iteration is independent!
    if (is_actually_smooth(...)) { ... }
}
```

Could use:
- OpenMP for simple parallelization
- Thread pool for better control
- SIMD for multiple primality tests simultaneously

### 5. Better Smoothness Bounds (Algorithmic)

The current approach checks all candidates that pass the sieve. Consider:
- Tighter initial bounds to reduce candidates
- Two-stage checking (quick test first, expensive verification second)
- Probabilistic filtering before full check

## Recommended Next Steps

1. **Profile `check_for_remaining_relations` specifically**:
   ```bash
   perf record -g -e cpu-clock --call-graph dwarf \
     ./gbin/lsieve 1000000 1000100
   perf report
   ```

2. **Focus on these functions**:
   - `is_actually_smooth()`
   - `is_probable_prime()`
   - `factorise_no_trial()`

3. **Implement early-exit optimizations** (low-hanging fruit):
   - Add quick rejection tests before expensive operations
   - Skip obviously non-smooth candidates

4. **Consider parallelization**:
   - The 49% spent in serial checking could become 5-10% with 8-core parallelization
   - Potential 40% overall speedup!

## Why Previous Optimizations Didn't Help

The micro-optimizations (precomputing powers, caching primes, optimizing ceil/floor) targeted phases that take <2% of total time each. Even a 50% improvement in these phases would only yield 1% overall speedup.

Meanwhile, a 50% improvement in `check_for_remaining_relations` would yield **25% overall speedup**.

## Expected Impact

| Optimization | Difficulty | Expected Speedup |
|--------------|-----------|------------------|
| Early exit tests | Easy | 10-20% in check phase = 5-10% overall |
| Batch processing | Medium | 20-30% in check phase = 10-15% overall |
| Parallel checking | Hard | 80-90% in check phase = 40-45% overall |
| **Combined** | | **60-70% overall speedup possible** |

This would take the current 760 relations/sec to **1200-1300 relations/sec**!
