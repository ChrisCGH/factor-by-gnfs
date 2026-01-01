# Deep Profiling Results with perf

## Executive Summary

**perf profiling reveals the actual CPU-level bottlenecks**, which differ from the high-level phase timers:

```
Function                                      CPU Time
----------------------------------------------------------
QS (Quadratic Sieve factorization)            23.7%  ← #1 BOTTLENECK
SieveCache::add (memory operations)            17.2%  ← #2 BOTTLENECK  
std::__introsort_loop (sorting)                10.2%  ← #3 BOTTLENECK
check_interval1                                 6.3%
GMP primality testing                           4.7%
is_probable_prime wrapper                       4.6%
```

## Key Insight: Timer vs CPU Profiling

The built-in timers showed "check relations" taking 49% of wall-clock time, but perf reveals:
- **QS factorization** (23.7%) is the biggest CPU consumer
- **Memory operations** (17.2%) in sieve cache are second
- **Sorting** (10.2%) is third

The discrepancy occurs because:
- Timers measure **wall-clock time** including I/O, memory stalls, etc.
- perf measures **actual CPU cycles** spent executing instructions

## Detailed Analysis

### 1. QS (Quadratic Sieve) - 23.7% of CPU time

**What it does:** Factors composite numbers using quadratic sieve algorithm.

**Why it's expensive:**
- Called from `is_actually_smooth()` via `factorise_no_trial()`
- Complex algorithm with many modular arithmetic operations
- Working with large numbers (VeryLong types)

**Location in code:**
```cpp
// LatticeSiever.cpp:873
if (!remaining_quotient.factorise_no_trial(&fac)) return false;
```

**Optimization opportunities:**
1. **Early rejection** - Check if number is obviously too large before factoring
2. **Pollard's rho first** - Try faster algorithms before QS
3. **Caching** - Cache factorizations of numbers we see repeatedly
4. **Batch factoring** - Accumulate numbers and factor together

**Expected impact:** 15-20% overall speedup if optimized

### 2. SieveCache::add() - 17.2% of CPU time

**What it does:** Adds sieve contributions to memory arrays.

**Why it's expensive:**
- Memory-bound operation touching large arrays
- Called millions of times during sieving
- Poor cache locality

**Location in code:**
```cpp
// Called from sieve1(), sieve2(), sieve1_again()
sieveCache_.add(ptr, f_span, e22, iter);
```

**Optimization opportunities:**
1. **Cache blocking** - Process data in L3-cache-sized chunks
2. **Prefetching** - Software prefetch for better memory pipeline
3. **SIMD** - Vectorize the additions if possible
4. **Reduce array size** - Tune sieve region to fit in cache

**Expected impact:** 10-15% overall speedup if optimized

### 3. std::__introsort_loop - 10.2% of CPU time

**What it does:** Sorting operation (likely factor lists).

**Why it's expensive:**
- Called frequently on factor lists
- Branch-heavy algorithm with unpredictable patterns

**Location in code:**
```cpp
// LatticeSiever.cpp:933
std::sort(factors1.begin(), factors1.end());
```

**Optimization opportunities:**
1. **Insertion sort for small lists** - Faster than introsort for n < 10
2. **Avoid re-sorting** - Keep lists sorted as we add elements
3. **Count sort** - If factors are in bounded range
4. **Parallel sort** - If lists are large enough

**Expected impact:** 5-8% overall speedup if optimized

### 4. check_interval1 - 6.3% of CPU time

**What it does:** Evaluates sieve locations to find smooth candidates.

**Why it's expensive:**
- Evaluates polynomial at many lattice points
- Memory access to sieve array

**Optimization opportunities:**
1. **Early exit** - Skip obviously non-smooth locations
2. **SIMD** - Vectorize polynomial evaluation
3. **Better thresholds** - Tune cutoff values

**Expected impact:** 3-5% overall speedup if optimized

### 5. GMP Primality Testing - 4.7% of CPU time

**What it does:** Miller-Rabin primality test from GMP library.

**Why it's expensive:**
- Modular exponentiation is slow
- Called for every remaining quotient

**Optimization opportunities:**
1. **Reduce iterations** - Use fewer Miller-Rabin rounds
2. **Trial division first** - Eliminate composites cheaply
3. **Deterministic test** - For numbers in range, use deterministic variant
4. **Batch testing** - Not easily parallelizable but could pipeline

**Expected impact:** 2-4% overall speedup if optimized

## Revised Optimization Strategy

Based on perf data, the priority order is:

### Priority 1: Optimize QS Factorization (23.7% CPU)

**Quick wins:**
```cpp
// Before factoring, check if too large
if (remaining_quotient > L * L) return false;

// Try Pollard's rho first (faster for small composites)
if (remaining_quotient < 1000000) {
    // Try pollard_rho with small number of iterations
    if (pollard_rho(remaining_quotient, factor)) {
        // Success with fast method
        return true;
    }
}

// Fall back to QS only if necessary
```

**Expected:** 15-20% overall speedup

### Priority 2: Optimize SieveCache Memory Access (17.2% CPU)

**Quick wins:**
```cpp
// Cache-blocking to fit in L3 (2 MiB)
const size_t CACHE_BLOCK = 256 * 1024; // 256KB blocks

// Process sieve in blocks
for (size_t block = 0; block < total_size; block += CACHE_BLOCK) {
    // Sieve this block with all primes
    for (auto prime : primes) {
        sieve_block(block, CACHE_BLOCK, prime);
    }
}
```

**Expected:** 10-15% overall speedup

### Priority 3: Optimize Sorting (10.2% CPU)

**Quick wins:**
```cpp
// Use insertion sort for small factor lists
if (factors.size() < 10) {
    insertion_sort(factors.begin(), factors.end());
} else {
    std::sort(factors.begin(), factors.end());
}

// Or keep factors sorted as we insert
void add_factor_sorted(vector<long>& factors, long p) {
    auto pos = std::lower_bound(factors.begin(), factors.end(), p);
    factors.insert(pos, p);
}
```

**Expected:** 5-8% overall speedup

## Combined Potential

| Optimization | CPU Reduction | Overall Speedup |
|--------------|---------------|-----------------|
| QS optimization | 23.7% → 12% | ~15% faster |
| Memory optimization | 17.2% → 10% | ~10% faster |
| Sorting optimization | 10.2% → 5% | ~7% faster |
| **Combined** | | **35-40% faster** |

This would take performance from **760 relations/sec to 1025-1065 relations/sec**.

## Implementation Roadmap

### Phase 1: QS Factorization (1-2 days)
1. Add early rejection tests
2. Implement Pollard's rho for small composites
3. Add factorization caching

### Phase 2: Memory Optimization (2-3 days)
1. Implement cache-blocking in sieve operations
2. Add software prefetching hints
3. Tune sieve region size

### Phase 3: Sorting Optimization (1 day)
1. Use insertion sort for small lists
2. Keep factors sorted during insertion
3. Profile to verify improvement

## How to Use This Data

1. **Focus on QS first** - Biggest impact (23.7% of CPU)
2. **Measure after each change** - Re-run perf to verify
3. **Don't optimize check_interval prematurely** - Only 6.3% of CPU

The profiling tools correctly identified that "check relations" was expensive (49% wall time), but perf reveals that within that phase, **QS factorization** is what's actually burning CPU cycles.

## Running More Detailed Profiling

To profile specific functions:
```bash
# Profile check_for_remaining_relations specifically
perf record -g -e cpu-clock --call-graph dwarf \
  -F 997 ./gbin/lsieve 1000000 1000100

# Generate annotated source
perf annotate --stdio QS > qs_annotate.txt

# Check cache misses specifically
perf stat -e cache-misses,cache-references,L1-dcache-load-misses \
  ./gbin/lsieve 1000000 1000100
```

## Conclusion

The built-in timers were correct that "check relations" is expensive, but perf reveals the **root cause**:
- 24% is QS factorization (algebraic operations)
- 17% is memory operations (cache-bound)
- 10% is sorting overhead

Optimizing these three areas has potential for **35-40% overall speedup**, compared to the 8-14% we expected from the original micro-optimizations.

This is why **profiling is essential** - assumptions about bottlenecks are often wrong!
