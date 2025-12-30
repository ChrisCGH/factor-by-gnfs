# Cache Performance Analysis

## Critical Finding: Severe Cache Miss Problem

**Cache performance counters reveal a critical memory bottleneck:**

```
Metric                          Value          Analysis
------------------------------------------------------------------------
Cache misses                    921,535,136    68.77% miss rate ← CRITICAL!
Cache references              1,339,971,894    
L1 dcache load misses           981,832,002    Most loads miss L1!
Total time                      47.0 seconds   
```

**68.77% cache miss rate is SEVERE** - this means 2 out of 3 memory accesses miss the cache and go to main memory, which is 100-200x slower than cache access.

## Impact Analysis

### What This Means

With a 68.77% cache miss rate:
- **Most memory accesses stall the CPU** waiting for data from RAM
- **SieveCache::add (17.2% CPU)** is likely much worse in wall-clock time due to cache misses
- Even CPU-bound operations like QS are slowed by waiting for operands

### Why Timer vs perf vs Cache Stats All Show Different Pictures

1. **Built-in timers (49% in check_relations)**: Wall-clock time including all stalls
2. **perf CPU profiling (24% in QS)**: Time CPU spends executing instructions
3. **Cache stats (69% miss rate)**: Time CPU spends **waiting** for memory

The real story: **CPU spends more time waiting for memory than executing instructions!**

## Root Cause Analysis

### Where Cache Misses Occur

Based on perf showing `SieveCache::add` at 17.2% CPU and cache miss rate of 69%, the sieve operations are the primary source:

1. **Sieve array access patterns** - Large sieve arrays touched in non-contiguous patterns
2. **Factor base iteration** - Jumping through factor base with irregular strides
3. **Candidate checking** - Random access to potentially smooth points

### Memory Hierarchy

```
Location        Size        Latency      Miss Penalty
-------------------------------------------------------
L1 cache        64 KiB      ~4 cycles    → L2 or RAM
L2 cache        512 KiB     ~12 cycles   → L3 or RAM
L3 cache        2 MiB       ~40 cycles   → RAM
RAM             16+ GB      ~200 cycles  ← 69% of accesses!
```

With 69% cache miss rate and 981M L1 misses, most accesses go to RAM at ~200 cycles each. This adds:
- 981M × 200 cycles ≈ 196 billion cycles of stalls
- At 3 GHz CPU: ~65 seconds of stall time
- Actual runtime: 47 seconds (includes overlapped computation)

**Memory stalls dominate the runtime!**

## Optimization Priorities - REVISED

The cache data completely changes the optimization priorities:

### Priority 1: Fix Cache Miss Rate (CRITICAL - 35-50% speedup potential)

**Target: Reduce 69% miss rate to <20%**

This is now the #1 priority, superseding QS optimization.

#### Strategy A: Cache-Blocking (Immediate, High Impact)

Current code likely sieves entire region with each prime, thrashing cache:

```cpp
// BAD: Sieves entire region with each prime (poor locality)
for (auto prime : factor_base) {
    for (location in entire_sieve_region) {
        sieve[location] += log(prime);
    }
}
```

**Fix with cache-blocking:**

```cpp
// GOOD: Process small blocks that fit in cache
const size_t L3_SIZE = 2 * 1024 * 1024;  // 2 MiB L3 cache
const size_t BLOCK_SIZE = L3_SIZE / 4;    // Leave room for other data

for (size_t block_start = 0; block_start < sieve_size; block_start += BLOCK_SIZE) {
    size_t block_end = std::min(block_start + BLOCK_SIZE, sieve_size);
    
    // Sieve this block with ALL primes before moving to next block
    for (auto prime : factor_base) {
        for (size_t loc = block_start; loc < block_end; loc += prime) {
            sieve[loc] += log(prime);
        }
    }
}
```

**Expected impact:**
- Reduce cache miss rate from 69% to 15-25%
- **35-50% overall speedup** (47s → 23-30s)
- Transforms memory-bound into CPU-bound workload

#### Strategy B: Reduce Sieve Region Size (Quick Win)

If sieve region > 2 MiB (L3 size), it won't fit in cache:

```cpp
// Check current sieve size
size_t sieve_bytes = (max_c - min_c + 1) * (max_d - min_d + 1) * sizeof(SIEVE_TYPE);

// If > 2 MiB, reduce region
if (sieve_bytes > 2 * 1024 * 1024) {
    // Split into multiple smaller sieves
    // OR reduce region bounds in config
}
```

**Expected impact:**
- If sieve > 2 MiB: 20-30% speedup from this alone
- If sieve < 2 MiB: Root cause is access pattern, use Strategy A

#### Strategy C: Prefetching (Medium Difficulty)

Add software prefetch hints for predictable access patterns:

```cpp
for (size_t i = 0; i < sieve_size; i += prime) {
    // Prefetch next few locations
    __builtin_prefetch(&sieve[i + prime], 1, 3);
    __builtin_prefetch(&sieve[i + prime*2], 1, 3);
    
    sieve[i] += log_prime;
}
```

**Expected impact:**
- 5-10% additional speedup
- Best combined with cache-blocking

### Priority 2: Optimize QS Factorization (15-20% speedup)

Still important but now secondary to cache issues:

1. **Early rejection** - Check if number too large before factoring
2. **Pollard's rho first** - Try fast algorithm for small composites
3. **Cache factorizations** - Avoid re-factoring same numbers

### Priority 3: Optimize Sorting (5-8% speedup)

Lower priority given cache problems:

1. **Insertion sort for small lists**
2. **Keep factors sorted during insertion**

## Implementation Roadmap - REVISED

### Phase 1: Cache Optimization (CRITICAL - 2-3 days)

**Day 1: Measure and understand**
```bash
# Find current sieve size
# Check if it fits in L3 cache (2 MiB)
# Profile memory access patterns
```

**Day 2: Implement cache-blocking**
- Modify `sieve_by_vectors1()` and `sieve_by_vectors2()`
- Process sieve in L3-cache-sized blocks
- Verify cache miss rate drops to <20%

**Day 3: Tune and validate**
- Optimize block size for best performance
- Add prefetching if helpful
- Measure: expect 35-50% speedup

### Phase 2: QS Factorization (1-2 days)

Only after cache issues are resolved:
- Implement early rejection
- Add Pollard's rho for small composites
- Expect: 15-20% additional speedup on already-improved baseline

### Phase 3: Sorting (1 day)

Final polish:
- Insertion sort for small lists
- Expect: 5-8% additional speedup

## Expected Results - UPDATED

```
Optimization               Speedup      Cumulative Time
---------------------------------------------------------
Baseline                   1.00×        47.0 seconds
Cache-blocking             1.50×        31.3 seconds  ← PHASE 1
QS optimization            1.18×        26.5 seconds  ← PHASE 2  
Sorting optimization       1.07×        24.8 seconds  ← PHASE 3
---------------------------------------------------------
Total                      1.90×        24.8 seconds

Relations/sec: 760 → 1440 (90% faster!)
```

This is **much better** than the 35-40% we estimated from CPU profiling alone, because we now understand that memory is the primary bottleneck.

## How to Validate

### Before Changes
```bash
perf stat -e cache-misses,cache-references,L1-dcache-load-misses \
  ./gbin/lsieve 1000000 1001000

# Current: 69% miss rate, 47 seconds
```

### After Cache-Blocking
```bash
perf stat -e cache-misses,cache-references,L1-dcache-load-misses \
  ./gbin/lsieve 1000000 1001000

# Target: <20% miss rate, <30 seconds
```

### Detailed Cache Analysis
```bash
# Profile specific cache events
perf stat -e L1-dcache-loads,L1-dcache-load-misses,\
           L1-dcache-stores,L1-dcache-store-misses,\
           LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses \
  ./gbin/lsieve 1000000 1001000

# Check which functions cause cache misses
perf record -e cache-misses -g ./gbin/lsieve 1000000 1001000
perf report
```

## Code Locations to Modify

### sieve_by_vectors1() - LatticeSiever.cpp:~1170

Current structure causes cache thrashing. Need to add blocking:

```cpp
void LatticeSiever::sieve_by_vectors1()
{
    // Determine sieve size
    size_t sieve_size = fixed_sieve_array_size;
    const size_t BLOCK_SIZE = 512 * 1024;  // 512 KB blocks
    
    // Process in cache-friendly blocks
    for (size_t block = 0; block < sieve_size; block += BLOCK_SIZE) {
        size_t block_end = std::min(block + BLOCK_SIZE, sieve_size);
        
        // Sieve this block with all primes
        auto iter = alg_factor_base_->begin();
        auto enditer = alg_factor_base_->end();
        for (; iter != enditer; ++iter) {
            // ... existing sieve logic but only for [block, block_end)
        }
    }
}
```

### sieve_by_vectors2() - LatticeSiever.cpp:~1277

Similar blocking needed here.

### SieveCache::add() - LatticeSiever.cpp or LatticeSiever.h

This function at 17.2% CPU with 69% miss rate is the main culprit. Cache-blocking will help this most.

## Key Insight

**The 69% cache miss rate explains everything:**

- Why "check relations" takes 49% wall time: Waiting for memory during factorization
- Why SieveCache::add takes 17% CPU time: Memory-bound operations
- Why QS takes 24% CPU but seems to take longer: Waiting for operands from RAM

**Fixing cache misses will speed up ALL phases**, not just sieving. The QS factorization will run faster with better cache behavior for its working data.

## Conclusion

The cache profiling data reveals that **memory access is the real bottleneck**, not computation. The 69% cache miss rate is critical and must be addressed first.

**Revised optimization strategy:**
1. **Cache-blocking (critical)**: 35-50% speedup
2. **QS optimization**: 15-20% additional speedup
3. **Sorting optimization**: 5-8% additional speedup

**Total potential: ~90% speedup (760 → 1440 relations/sec)**

This completely changes the optimization approach from "optimize expensive functions" to "fix memory access patterns first, then optimize functions."
