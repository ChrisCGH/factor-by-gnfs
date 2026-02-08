# Auto-Tuning Feature: DISABLED

## Summary

The auto-tuning feature has been **completely disabled** after multiple failed attempts to make it work correctly resulted in severe performance degradation.

## Performance History

| State | Relations/sec | Status |
|-------|---------------|--------|
| Baseline (no auto-tuning) | 963 | ✅ Normal |
| Auto-tuning attempt #1 | 80 | ❌ 91% slower |
| Auto-tuning attempt #2 | 14 | ❌ 98.5% slower! |
| **After disabling** | **963** | **✅ Restored** |

## What Happened

### Attempt #1: Wrong Direction
The RateOptimizer suggested adjustments in the wrong direction:
- When performance degraded, suggested +1 for SIEVE_BOUND_ADJUSTMENT
- With `cutoff -= adjustment`, this made cutoff LOWER
- Lower cutoff → more lenient → MORE candidates → SLOWER
- Created death spiral: 91% slowdown

### Attempt #2: Overcorrection
Inverted the exploration logic to suggest negative adjustments:
- Suggested -1 for SIEVE_BOUND_ADJUSTMENT, -5 for INITIAL_CUTOFF
- These were ADDED to parameters: `param += adjustment`
- Decreased both parameters toward their minimums
- Made BOTH filtering stages extremely strict
- Result: Almost nothing passed filters → 98.5% slowdown!

## Root Causes

1. **Complex semantic interaction**
   - Parameters are modified by ADDING adjustments: `param += adj`
   - But cutoff uses SUBTRACTION: `cutoff -= param`
   - This double-negation is confusing and error-prone

2. **Cumulative effects not accounted for**
   - Each adjustment accumulates with previous ones
   - Small adjustments over many iterations compound
   - Exploration logic didn't consider accumulated state

3. **No clear optimal direction**
   - Higher SIEVE_BOUND_ADJUSTMENT can make things faster OR slower
   - Depends on initial values, problem characteristics, phase
   - No simple "increase to improve" relationship

4. **Multiple interdependent parameters**
   - SIEVE_BOUND_ADJUSTMENT1, SIEVE_BOUND_ADJUSTMENT2, INITIAL_CUTOFF
   - All interact in complex ways
   - Optimizing one can make others worse

## The Solution: Complete Disablement

**What was changed:**
- `LatticeSiever::apply_parameter_adjustments()` now returns immediately
- No parameters are ever modified
- Rate tracking still works (monitoring only)
- Clear documentation explains why disabled

**Result:**
- Performance restored to baseline ~963 rel/sec
- Stable, predictable behavior
- No more death spirals
- No more experimentation disasters

## What Still Works

### Rate Tracking (RateOptimizer)
The monitoring infrastructure is still functional:
- ✅ Records relation generation rates
- ✅ Calculates moving averages and trends
- ✅ Displays statistics in verbose mode
- ✅ Tracks best-known rate
- ❌ Does NOT apply parameter adjustments

### ENABLE_AUTO_TUNING Flag
The configuration flag still works:
- When `true`: Enables rate tracking and monitoring
- Adjustments are suggested but **IGNORED**
- Safe to keep enabled for analysis

## Future Work

To make auto-tuning work properly would require:

1. **Fundamental redesign**
   - Clear mathematical model of parameter effects
   - Proper optimization framework (not ad-hoc exploration)
   - Validation that adjustments improve performance

2. **Better exploration strategy**
   - Account for cumulative effects
   - Use gradient descent or similar proven algorithm
   - Test adjustments before permanently applying

3. **Extensive testing**
   - Test on diverse problem instances
   - Verify improvements are real and reproducible
   - Ensure no scenarios cause severe degradation

4. **Simpler parameters**
   - Reduce number of tunable parameters
   - Use orthogonal parameters (don't interact)
   - Clear monotonic relationship with performance

## Recommendation

**Keep auto-tuning DISABLED indefinitely.**

The feature is fundamentally broken and every attempt to fix it has made things worse. The baseline performance without auto-tuning is excellent (~963 rel/sec). There's no need to risk breaking it again with speculative optimizations.

If someone wants to tackle auto-tuning in the future, they should:
1. Start from scratch with proper design
2. Test extensively before deploying
3. Have a clear mathematical understanding of parameter effects
4. Include extensive validation and rollback mechanisms

## Status

- ✅ Auto-tuning disabled
- ✅ Performance restored to baseline
- ✅ Rate tracking still available for monitoring
- ✅ System stable and reliable
- ✅ No further action needed

**The problem is solved. Performance is back to normal.**
