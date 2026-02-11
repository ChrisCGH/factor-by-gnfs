# Relation Rate Tracking and Auto-Tuning Guide

## Overview

The lattice siever now includes automatic tracking of relation generation rates and can dynamically adjust sieving parameters to maximize throughput.

## Features

### Rate Tracking
- Tracks relation generation rate (relations per second) across sieving operations
- Calculates moving averages and trends
- Identifies best-performing parameter configurations
- Detects degrading performance

### Parameter Auto-Tuning
When enabled, the system automatically adjusts:
- **B1, B2**: Factor base bounds (algebraic and rational)
- **SIEVE_BOUND_ADJUSTMENT1/2**: Sieve threshold adjustments
- **INITIAL_CUTOFF**: Initial sieve cutoff value

### Safety Bounds
Parameters are constrained to safe ranges:
- B1: 100,000 to 5,000,000
- B2: 50,000 to 3,000,000
- SIEVE_BOUND_ADJUSTMENT: 0 to 50
- INITIAL_CUTOFF: 5 to 100

## Configuration

### Enabling Auto-Tuning

Add to your sieve configuration file (e.g., `sieve.cfg`):

```
ENABLE_AUTO_TUNING = true
```

To disable (default):

```
ENABLE_AUTO_TUNING = false
```

Or simply omit the line.

### Example Configuration

```
# Configuration file for GNFS siever
SIEVE_ID = test_problem
N = 123456789...
f1 = Polynomial
...
B1 = 1000000
B2 = 700000
SIEVE_BOUND_ADJUSTMENT1 = 11
SIEVE_BOUND_ADJUSTMENT2 = 0
INITIAL_CUTOFF = 10

# Enable automatic parameter tuning
ENABLE_AUTO_TUNING = true

# Enable verbose output to see tuning in action
```

## How It Works

### Rate Tracking Algorithm

1. **Recording**: After each special-q sieving, records:
   - Current relation rate (relations/second)
   - Total elapsed time
   - Number of relations found
   - Current parameter values

2. **Analysis**: Maintains a sliding window of recent samples (default: 20)
   - Calculates moving average over last 5 samples
   - Computes rate trend using linear regression over last 10 samples
   - Identifies best rate achieved and corresponding parameters

3. **Trend Detection**:
   - Positive trend: Rate improving over time
   - Negative trend: Rate degrading
   - Stable: Low variance in recent rates (coefficient of variation < 10%)

### Auto-Tuning Algorithm

When auto-tuning is enabled:

1. **Trigger Conditions**: Adjustments suggested when:
   - At least 10 samples collected (need baseline data)
   - At least 5 samples since last adjustment (avoid thrashing)
   - Rate is degrading (negative trend < -0.5) OR
   - Current rate is significantly below best (< 80% of best)

2. **Adjustment Strategy**:
   - If best configuration known: Move gradually (10% steps) towards it
   - If no clear best yet: Small exploratory adjustments
   - Example exploratory changes:
     - B1: ±1000
     - SIEVE_BOUND_ADJUSTMENT: ±1
     - INITIAL_CUTOFF: ±5

3. **Application**: Parameters updated with safety bounds enforced

4. **Convergence**: If rate stable and improving, no changes made (let it run)

## Output

### Verbose Mode

When running with verbose output (set `LATTICE_SIEVER_VERBOSE=1` environment variable), you'll see:

```
5 -> 123 -> 8 (850.5,3062,73489), (978.2,3521,84514), [trend:2.3, avg:975.5, best:1005.2]
```

Breaking down:
- `5 -> 123 -> 8`: Potentially smooth counts at each stage
- `(850.5,3062,73489)`: Current special-q rate (per sec, hour, day)
- `(978.2,3521,84514)`: Running average rate
- `[trend:2.3, avg:975.5, best:1005.2]`: Optimization metrics
  - `trend`: Rate trend (positive = improving)
  - `avg`: Moving average
  - `best`: Best rate achieved

### When Adjusting

```
### Auto-tuning: Adjusting parameters ###
Rate trend: -0.8
Current rate: 850.2 rel/sec
Best rate: 1005.2 rel/sec
Adjustments: B1=500, B2=200, SBA1=1, SBA2=0, IC=2
New parameters: B1=1000500, B2=700200, SBA1=12, SBA2=0, IC=12
###################################
```

## Performance Considerations

### Factor Base Rebuilding

**Important**: Changes to B1 and B2 require rebuilding the factor bases to take full effect. Currently, the implementation:
- Records the parameter change
- Applies it for logging/tracking
- Full effect requires restarting the sieving session

Future enhancement: Dynamic factor base rebuilding.

### Other Parameters

SIEVE_BOUND_ADJUSTMENT and INITIAL_CUTOFF changes take effect immediately as they only affect sieving thresholds.

## Testing

A test suite is included in `RateOptimizerTest.cpp`:

```bash
cd gnfs
g++ -std=c++14 -I. -o RateOptimizerTest RateOptimizerTest.cpp
./RateOptimizerTest
```

Tests verify:
- Basic rate recording
- Trend calculation (increasing/decreasing)
- Stability detection
- Adjustment suggestions
- Enable/disable functionality

## Tuning the Auto-Tuner

If you want to modify the auto-tuning behavior, see `RateOptimizer.h`:

Key parameters:
- `history_size`: Number of samples to retain (default: 20)
- Trend window: Number of samples for trend calculation (default: 10)
- Moving average window: Number of samples for average (default: 5)
- Adjustment frequency: Minimum samples between adjustments (default: 5)
- Stability threshold: Coefficient of variation for stability (default: 0.15)

## Examples

### Conservative Auto-Tuning

For problems where you want to be cautious:

```
ENABLE_AUTO_TUNING = true
# Start with well-tested parameters
B1 = 1000000
B2 = 700000
```

The system will only adjust if performance degrades significantly.

### Aggressive Exploration

For exploratory runs where you want to find optimal parameters:

```
ENABLE_AUTO_TUNING = true
# Start with smaller values, let it explore
B1 = 500000
B2 = 350000
```

Monitor the output to see which parameters give best results.

### Manual Optimization

Run with auto-tuning disabled but watch the rate statistics:

```
ENABLE_AUTO_TUNING = false
```

Note the rate outputs, then manually adjust parameters between runs based on what works best.

## Troubleshooting

### Rate Not Improving

If auto-tuning is enabled but rates aren't improving:

1. **Check sample count**: Need at least 10 samples before adjustments start
2. **Look for stability**: If rate is stable and good, it won't adjust
3. **Parameter bounds**: May have hit min/max bounds
4. **Verbose output**: Enable to see what's happening

### Unstable Behavior

If parameters keep changing:

1. **Increase adjustment frequency**: Modify `sample_count_ - adjustment_count_ < 5` to larger value
2. **Tighten stability threshold**: Reduce from 0.15 to 0.10
3. **Disable during critical runs**: Turn off auto-tuning for production factorization

### Performance Degradation

If performance gets worse with auto-tuning:

1. **Disable it**: `ENABLE_AUTO_TUNING = false`
2. **Check verbose output**: See what parameters it's trying
3. **Use known-good parameters**: Set them explicitly in config

## Future Enhancements

Potential improvements:
- Dynamic factor base rebuilding for B1/B2 changes
- Multi-objective optimization (relations + quality)
- Machine learning-based parameter prediction
- Cross-problem parameter transfer learning
- Adaptive exploration vs. exploitation strategies
- Save/load parameter history for resume

## References

The auto-tuning algorithm uses a simple hill-climbing approach with:
- Linear regression for trend detection
- Coefficient of variation for stability
- Gradual parameter adjustment (10% steps)
- Safety bounds to prevent extremes

For more sophisticated parameter optimization, consider:
- Bayesian optimization
- Genetic algorithms
- Reinforcement learning

But the simple approach often works well for this problem domain.
