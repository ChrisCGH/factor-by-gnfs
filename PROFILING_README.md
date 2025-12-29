# Performance Profiling Tools for LatticeSiever

This directory contains tools to help profile and optimize the LatticeSiever performance.

## Quick Start

### Step 1: Run Performance Analysis

```bash
./analyze_performance.sh [min_q] [max_q]
```

Example:
```bash
./analyze_performance.sh 1000000 1000010
```

This will:
- Build lsieve if needed
- Run a test sieve
- Display built-in timer results
- Show performance summary

### Step 2: Analyze Timing Data

```bash
./analyze_timings.py [path/to/sieve.tim]
```

Example:
```bash
./analyze_timings.py gnfs/sieve.tim
```

This provides:
- Visual breakdown of time spent in each phase
- Specific optimization recommendations for bottlenecks
- System cache information
- Next steps for deep profiling

## Understanding the Results

The LatticeSiever has built-in timers that track:

1. **sieve by vectors 1** - Algebraic sieving (scanning memory with prime patterns)
2. **check interval 1** - Checking which sieve locations might be smooth
3. **sieve by vectors 2** - Rational sieving
4. **check interval 2** - Checking rational smoothness
5. **remove factors** - Trial division to find actual factors
6. **eliminate** - Filtering out non-smooth candidates
7. **check relations** - Final validation

### What Different Patterns Mean

**Pattern 1: Sieve by vectors dominates (>40% total time)**
- This is memory-bound work
- Optimization strategies:
  - Reduce sieve region size to fit in cache
  - Improve memory access patterns
  - Consider SIMD vectorization
- Not fixed by compiler optimization!

**Pattern 2: Check interval dominates (>30% total time)**
- Evaluating which candidates might be smooth
- Optimization strategies:
  - Better early-exit strategies
  - Optimize polynomial evaluation
  - Reduce branch mispredictions
- May benefit from algorithmic changes

**Pattern 3: Remove factors dominates (>25% total time)**
- Trial division is expensive
- Our small prime caching should help
- Further optimizations:
  - Wheel factorization
  - Montgomery multiplication
  - Better divisibility tests

**Pattern 4: Balanced time across phases**
- No single bottleneck - already well-optimized
- Small improvements in multiple areas needed
- Focus on cache efficiency and algorithmic improvements

## Deep Profiling with perf

For Linux systems, use `perf` for detailed analysis:

```bash
cd gnfs

# Record profile data
perf record -g --call-graph dwarf ./gbin/lsieve 1000000 1000100

# View interactive report
perf report

# Or generate annotated source
perf annotate --stdio > perf_annotate.txt
```

Key metrics:

```bash
# Check cache performance
perf stat -e cache-misses,cache-references,L1-dcache-load-misses,LLC-load-misses \
  ./gbin/lsieve 1000000 1000100

# Check branch prediction
perf stat -e branch-misses,branch-instructions,branch-load-misses \
  ./gbin/lsieve 1000000 1000100
```

### Interpreting perf Results

**Cache miss rate > 10%**: Memory-bound, optimize data locality
**Branch miss rate > 5%**: Too many unpredictable branches
**High cycles per instruction**: CPU stalls, likely memory/cache issues

## What the Compiler Already Optimizes

With `-O3 -march=native -flto`, the compiler automatically does:

✅ Loop unrolling and vectorization
✅ Function inlining (including across files with LTO)
✅ Constant propagation and folding
✅ Dead code elimination
✅ Common subexpression elimination
✅ Loop invariant code motion

These are **NOT** effective optimizations anymore:
- ❌ Manually hoisting loop invariants → compiler already does this
- ❌ Small inline functions → LTO inlines everything
- ❌ Local micro-optimizations → compiler optimizes better

These **ARE** effective optimizations:
- ✅ Algorithmic improvements (better algorithms)
- ✅ Cache-conscious data structures
- ✅ Reducing total work (skip unnecessary computation)
- ✅ Memory access pattern improvements
- ✅ SIMD explicit vectorization (for complex patterns)

## Real Optimization Examples

Based on profiling results, here are proven optimization strategies:

### 1. Cache-Blocking for Sieving

If "sieve by vectors" takes >40% of time:

```cpp
// Instead of sieving entire region at once
for (all primes p) {
    for (all locations in region) {
        sieve[location] += log(p);
    }
}

// Use cache-blocking
const int BLOCK_SIZE = 8192; // Fits in L1 cache
for (block_start = 0; block_start < region_size; block_start += BLOCK_SIZE) {
    for (all primes p) {
        for (locations in [block_start, block_start+BLOCK_SIZE)) {
            sieve[location] += log(p);
        }
    }
}
```

### 2. Early Exit in Smoothness Checking

If "check interval" takes >30% of time:

```cpp
// Add early exit for obvious non-smooth cases
if (abs(value) > threshold * 10) {
    continue; // Skip expensive evaluation
}
```

### 3. Batch Processing

Process multiple candidates together for better cache utilization:

```cpp
// Instead of: for each candidate { check, factor, eliminate }
// Do: for each candidate { check }, for each candidate { factor }, ...
```

## Tools Summary

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `analyze_performance.sh` | Quick performance snapshot | First step, always |
| `analyze_timings.py` | Detailed phase analysis | After each run |
| `perf record/report` | CPU profiling | Identify hot functions |
| `perf stat` | Cache/branch analysis | Memory-bound issues |
| `valgrind --tool=callgrind` | Detailed call graph | Complex bottlenecks |
| `gprof` | Traditional profiling | If perf unavailable |

## Expected Performance

Typical LatticeSiever performance on modern hardware:
- Small numbers (100-150 digits): 100-500 relations/sec
- Medium numbers (150-200 digits): 50-200 relations/sec  
- Large numbers (200+ digits): 10-100 relations/sec

Performance depends on:
- CPU speed and cache size
- Factor base size
- Smoothness bounds
- Sieve region size
- Memory bandwidth

## Getting Help

1. Check `PROFILING_GUIDE.md` for detailed instructions
2. Run `./analyze_timings.py` to get specific recommendations
3. Use `perf` to identify actual hot spots
4. Focus on phases taking >20% of total time

Remember: **Profile first, optimize second!**
