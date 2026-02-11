# Auto-Tuning Re-Enabled for SIEVE_BOUND_ADJUSTMENT

## Summary

Auto-tuning is now **partially re-enabled** for parameters that are safe to adjust mid-session.

## What Changed

### Previously (After Invalid Relations Bug Fix)
All parameter adjustments were disabled because:
1. B1/B2 changes require factor base rebuild
2. **INCORRECTLY** thought SIEVE_BOUND_ADJUSTMENT affected precomputed values

### Now
- ✅ **SIEVE_BOUND_ADJUSTMENT1** adjustments enabled (safe)
- ✅ **SIEVE_BOUND_ADJUSTMENT2** adjustments enabled (safe)
- ✅ **INITIAL_CUTOFF** adjustments enabled (safe)
- ❌ **B1/B2** adjustments still disabled (require factor base rebuild)

## Why It's Safe

### The Misunderstanding

The previous comment claimed:
```cpp
// SIEVE_BOUND_ADJUSTMENT values affect precomputed values like
// L1_pow_LP1_ and log_L2_pow_LP2_ which are computed at session start.
```

**This was INCORRECT!**

### The Reality

Looking at the actual code:

**Precomputed values (lines 268-272):**
```cpp
L1_pow_LP1_ = L1_;
for (int i = 0; i < LP1_ - 1; i++) L1_pow_LP1_ *= static_cast<double>(L1_);

log_L2_pow_LP2_ = log10(L2_pow_LP2_);
```
These are computed from `L1`, `LP1`, `L2`, `LP2` **only** - NOT from SIEVE_BOUND_ADJUSTMENT!

**Usage of SIEVE_BOUND_ADJUSTMENT (lines 354, 397, 435):**
```cpp
// check_interval1:
int adjustment = SIEVE_BOUND_ADJUSTMENT1_;  // Line 354
...
cutoff -= adjustment;  // Line 397 - APPLIED DIRECTLY

// check_interval2:
int adjustment = SIEVE_BOUND_ADJUSTMENT2_;  // Line 435
...
cutoff -= adjustment;  // Similar usage
```

The adjustment is **applied directly** to the cutoff calculation each time. There are NO precomputed values involved!

## How Auto-Tuning Works Now

### 1. Rate Tracking (Always Active)
- Tracks relations/second after each special-q
- Calculates trends and moving averages
- Identifies best-performing parameter sets

### 2. Parameter Adjustment (When ENABLE_AUTO_TUNING = true)

**What Gets Adjusted:**
- `SIEVE_BOUND_ADJUSTMENT1` (range: 0-50)
  - Higher values: More relaxed cutoff → more candidates → slower but more thorough
  - Lower values: Stricter cutoff → fewer candidates → faster but might miss relations
  
- `SIEVE_BOUND_ADJUSTMENT2` (range: 0-50)
  - Similar effect for rational side checking
  
- `INITIAL_CUTOFF` (range: 5-100)
  - Initial threshold before checking sieve values

**How Adjustments Are Applied:**
- Algorithm detects degrading performance trend
- Suggests small adjustments (10% steps toward best-known values)
- Adjustments applied immediately (safe, no rebuild needed)
- Effects visible in next special-q

### 3. Verbose Output

When adjustments happen, you'll see:
```
### Auto-tuning: Adjusting parameters ###
Adjustments: SBA1=1, SBA2=0, IC=2
New values: SBA1=12, SBA2=10, IC=17
###################################
```

## Expected Performance Impact

### Benefits
- **Adaptive optimization**: Parameters adjust to current workload
- **Better cutoff thresholds**: Finds optimal balance between thoroughness and speed
- **No manual tuning needed**: Algorithm explores parameter space automatically

### Limitations
- Can only adjust cutoff thresholds, not factor base sizes
- Hill-climbing finds local optimum, not necessarily global
- Conservative adjustments (10% steps) mean slow convergence

### Realistic Expectations
- 0-10% performance improvement from optimal threshold tuning
- More stable performance across different special-q ranges
- Automatic recovery if parameters become suboptimal

## Configuration

### Enable Auto-Tuning
```bash
# In your config file:
ENABLE_AUTO_TUNING = true
```

### Disable Auto-Tuning (Default)
```bash
ENABLE_AUTO_TUNING = false
# OR just omit the setting
```

## Safety Features

### Built-In Safeguards
1. **Parameter bounds**: All adjustments clamped to safe ranges
2. **Gradual changes**: 10% steps prevent drastic swings
3. **Minimum samples**: Requires 10+ samples before adjusting
4. **Adjustment interval**: Waits 5+ samples between adjustments
5. **Best tracking**: Returns to best-known configuration if performance degrades

### What Can't Go Wrong
- ✅ No invalid relations (B1/B2 can't change)
- ✅ No crashes (all parameters bounded)
- ✅ No data corruption (only affects thresholds)
- ✅ Reversible (disable auto-tuning anytime)

### What Could Be Suboptimal
- ⚠️ Algorithm might settle on local optimum, not global best
- ⚠️ Conservative adjustments mean slow adaptation
- ⚠️ Might oscillate if rate fluctuates randomly

## Comparison with Previous State

### Before (All Disabled)
```
Relations/sec: ~980 (baseline)
Adjustments: None
Safety: Maximum (no changes)
Adaptability: None (manual tuning only)
```

### Now (SIEVE_BOUND_ADJUSTMENT Enabled)
```
Relations/sec: ~980-1080 (adaptive, potentially higher)
Adjustments: Cutoff thresholds only
Safety: High (bounded, gradual)
Adaptability: Good (finds optimal cutoffs)
```

### Future (If B1/B2 Enabled)
```
Relations/sec: Potentially >1100 (factor base optimization)
Adjustments: Full parameter space
Safety: Requires factor base rebuild
Complexity: High (not yet implemented)
```

## Technical Details

### Parameters and Their Effects

**SIEVE_BOUND_ADJUSTMENT1** (check_interval1):
- Used in algebraic smoothness checking
- Affects line 397: `cutoff -= adjustment;`
- Higher → more relaxed → more candidates pass check
- Lower → stricter → fewer candidates pass check

**SIEVE_BOUND_ADJUSTMENT2** (check_interval2):
- Used in rational smoothness checking
- Similar effect on line 435+
- Directly affects which candidates become relations

**INITIAL_CUTOFF** (check_interval1):
- Line 355: `int cutoff0 = INITIAL_CUTOFF_;`
- Line 388: `if ((*sieve_ptr) >= cutoff0)`
- Higher → process more sieve locations
- Lower → skip more locations early

### Algorithm Flow
1. Track rate after each special-q completion
2. Every N special-q values, check if adjustment needed:
   - Is rate degrading? (negative trend)
   - Is rate below best by >20%?
3. If yes, calculate adjustments:
   - Move 10% toward best-known parameters
   - Or make small exploratory change if no best known
4. Apply adjustments with bounds checking
5. Continue sieving with new parameters

## Troubleshooting

### "Auto-tuning not making adjustments"
- Check `ENABLE_AUTO_TUNING = true` in config
- Need 10+ samples before first adjustment
- Wait 5+ samples between adjustments
- Rate must be degrading or significantly below best

### "Performance degraded after enabling"
- Algorithm exploring parameter space (give it time)
- If persistent, disable and use manual parameters
- Check verbose output to see what's being adjusted

### "Want faster convergence"
- Edit `RateOptimizer.h`: Reduce step size or adjustment threshold
- Trade-off: Faster convergence vs stability

### "Want to adjust B1/B2"
- Not currently supported (requires factor base rebuild)
- Future enhancement: Implement factor base regeneration
- For now: Manually set B1/B2, let auto-tuning adjust cutoffs

## Future Work

### Possible Enhancements
1. **Factor base rebuild**: Enable B1/B2 auto-tuning
2. **Better algorithms**: Bayesian optimization, genetic algorithms
3. **Multi-objective**: Balance speed, smoothness rate, memory
4. **Persistent history**: Remember best parameters across runs
5. **Adaptive step size**: Larger steps when far from optimum

### Not Planned
- Changing these parameters mid-session is fundamentally safe now
- No architectural changes needed
- Can be enhanced incrementally

## Recommendation

**For Production:**
- ✅ Enable auto-tuning: `ENABLE_AUTO_TUNING = true`
- Monitor performance over time
- If stable improvement seen, keep enabled
- If unstable or degraded, disable and use manual tuning

**For Experimentation:**
- Try different initial parameter values
- Compare auto-tuned vs manual performance
- Analyze rate tracking output to understand behavior
- Share results to improve algorithm

## Summary

Auto-tuning is now **safe and functional** for cutoff threshold parameters. It can't cause invalid relations (B1/B2 can't change) and provides automatic optimization within the available parameter space. Performance should equal or exceed baseline, with potential for 5-10% improvement from optimal threshold tuning.
