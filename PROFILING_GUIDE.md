# LatticeSiever Performance Profiling Guide

This guide explains how to profile the LatticeSiever to identify real performance bottlenecks.

## Quick Start: Built-in Timing

The code already has built-in timers! The `Timer` class in `LatticeSiever.cpp` tracks time spent in each phase:

1. **sieve by vectors 1** - Algebraic factor base sieving
2. **check interval 1** - Checking algebraic smoothness candidates
3. **sieve by vectors 2** - Rational factor base sieving  
4. **check interval 2** - Checking rational smoothness candidates
5. **remove factors for rational** - Trial division by rational sieved factors
6. **eliminate rational** - Eliminating non-smooth rational candidates
7. **remove factors for algebraic** - Trial division by algebraic sieved factors
8. **eliminate algebraic** - Eliminating non-smooth algebraic candidates
9. **check relations** - Final relation validation

### View Timing Output

The timing data is written to `sieve.tim`. To see where time is being spent:

```bash
cd gnfs
cat sieve.tim
```

This will show cumulative time for each phase. The phases with the highest time are the bottlenecks.

## Method 1: Linux `perf` Profiling

For detailed CPU profiling on Linux:

```bash
# Build with profiling-friendly flags (debug symbols retained)
cd gnfs
make clean
make lsieve

# Run with perf to collect samples
perf record -g ./gbin/lsieve <min_q> <max_q>

# Generate report
perf report

# Or annotate source code
perf annotate
```

Key metrics to look for:
- **High sample count functions** - Where CPU time is spent
- **Cache misses** - `perf stat -e cache-misses ./gbin/lsieve ...`
- **Branch mispredictions** - `perf stat -e branch-misses ./gbin/lsieve ...`

## Method 2: gprof Profiling

To use traditional gprof profiling:

1. Edit `makefile` and uncomment the PROFILE line:
   ```makefile
   PROFILE = -pg
   ```

2. Rebuild:
   ```bash
   make clean
   make lsieve
   ```

3. Run the program:
   ```bash
   ./gbin/lsieve <min_q> <max_q>
   ```

4. Analyze results:
   ```bash
   gprof ./gbin/lsieve gmon.out > profile.txt
   less profile.txt
   ```

## Method 3: Valgrind Callgrind

For detailed call graph and cache analysis:

```bash
# Build normally (with -g for debug symbols)
make lsieve

# Run with callgrind
valgrind --tool=callgrind --dump-instr=yes --collect-jumps=yes \
  ./gbin/lsieve <min_q> <max_q>

# Visualize with kcachegrind
kcachegrind callgrind.out.*
```

Or use `callgrind_annotate` for text output:
```bash
callgrind_annotate --auto=yes callgrind.out.*
```

## Understanding the Results

### Common Bottlenecks in Sieving

Based on typical GNFS siever performance:

1. **Memory bandwidth** - Sieving touches lots of memory
   - Look for cache miss rates > 5%
   - L3 cache misses are expensive (>100 cycles)

2. **Integer modular arithmetic** - `modasm`, `mulmodasm2`, `inverse`
   - These are in hot loops in `sieve_by_vectors1/2`
   - May benefit from Montgomery multiplication

3. **Trial division** - `divide_by_small_primes1/2`
   - Dividing by many small primes
   - Our optimization should help here, but verify

4. **Vector generation** - `generate_ef_lattice`, lattice basis reduction
   - LLL-style reductions in tight loops
   - May benefit from better algorithms

5. **Bit array operations** - `sieve_bit_array_.set/isSet`
   - Frequent in `check_interval1/2`
   - Cache-friendly access patterns help

### What to Look For

After profiling, identify:

1. **Functions with >10% total time** - Primary optimization targets
2. **Functions called millions of times** - Small savings multiply
3. **Cache miss hotspots** - Memory access patterns to optimize
4. **Branch misprediction hotspots** - Unpredictable branches to reduce

## Interpreting Built-in Timers

Run a test sieve and check `sieve.tim`:

```bash
cd gnfs
./gbin/lsieve 1000000 1000100
cat sieve.tim
```

Look for phases taking >20% of total time. Common patterns:

- **sieve by vectors dominates** → Optimize `sieve1/sieve2` loops
- **check interval dominates** → Optimize smooth candidate checking
- **remove factors dominates** → Trial division is expensive
- **eliminate dominates** → Factor checking logic needs work

## Next Steps After Profiling

Once you've identified the real bottlenecks:

1. **Memory-bound** → Improve cache locality, reduce memory footprint
2. **CPU-bound in arithmetic** → Use faster algorithms (Montgomery, etc.)
3. **Branch-heavy** → Reduce unpredictable branches, add likely/unlikely hints
4. **I/O-bound** → Buffer writes, optimize file I/O patterns

The compiler with `-O3 -march=native -flto` already does:
- Loop unrolling and vectorization
- Function inlining
- Constant propagation
- Dead code elimination

Real improvements require:
- Better algorithms
- Improved memory access patterns
- Reduced work (skip unnecessary computation)
- Hardware-specific optimizations (SIMD, etc.)
