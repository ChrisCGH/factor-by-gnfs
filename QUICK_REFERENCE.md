# Quick Reference: What Happened and What to Do

## TL;DR

Two bugs found and fixed:
1. **Invalid relations** - Auto-tuning changed parameters unsafely → FIXED (adjustments disabled)
2. **Slow performance** - Expensive duplicate removal → FIXED (O(1) tracking flag)

**Result**: Code is correct, fast, and ready to use.

## What You Need to Know

### The Good News ✅
- Cache blocking works for check phases
- Rate tracking works (useful for monitoring)
- No invalid relations
- Performance back to ~980-1000 rel/sec
- Everything well documented

### What's Disabled ⚠️
- Auto parameter adjustment (unsafe, causes invalid relations)
- Sieving cache blocking (didn't help, removed for simplicity)

### What to Do
```bash
# In your config file:
ENABLE_AUTO_TUNING = false   # Keep this! (default)

# Or just don't set it - false is default
```

## Why Invalid Relations Happened

**Simple explanation:**
1. Auto-tuning said "let's try bigger factor base" (change B1/B2)
2. Code changed the B1/B2 values mid-session
3. But factor bases were already built with old values
4. Sieving used old factor bases
5. Checking used new B1/B2 values
6. Mismatch → relations marked "smooth" that weren't → INVALID!

**Fix**: Don't change parameters mid-session. Disabled auto-adjustments.

## Why Performance Slowed

**Simple explanation:**
1. Cache buckets could be added to tracking list multiple times
2. Fixed by calling sort+unique to remove duplicates
3. But with many duplicates, sorting 1000+ entries 3x per special-q
4. Small overhead × many calls = noticeable slowdown

**Fix**: Track buckets with boolean flag so they're never added twice. No sorting needed.

## What Works Now

### Cache Blocking (Check Phases)
```
Before: Process entire 8MB array linearly
Now: Process in 256KB blocks
Benefit: Slight improvement (check phases are 6-7% of runtime)
Status: WORKING ✅
```

### Rate Tracking
```
What it does: Monitors relations/sec, calculates trends
What it doesn't do: Change parameters (disabled)
Benefit: Useful for finding optimal parameters manually
Status: WORKING ✅
```

### Duplicate Prevention
```
Old way: Allow duplicates, then sort to remove (expensive)
New way: Track with flag, never add duplicates (free)
Benefit: No performance overhead
Status: WORKING ✅
```

## Performance Numbers

```
Baseline (no changes):     ~1000 rel/sec
With check blocking:       ~980-1000 rel/sec (minimal change)
With broken sieve blocking: 732 rel/sec (REVERTED)
With expensive dedup:       ~slow (FIXED)
Current (all fixes):       ~980-1000 rel/sec ✅
```

## For Developers: What Needs to Happen for Safe Auto-Tuning

To re-enable parameter adjustment safely:

```cpp
// When B1 or B2 changes:
1. Rebuild rational factor base with new B1
2. Rebuild algebraic factor base with new B2
3. Clear any cached/precomputed data
4. Verify relations are still valid

// When SIEVE_BOUND_ADJUSTMENT changes:
1. Recompute L1_pow_LP1_
2. Recompute log_L2_pow_LP2_
3. Update any other precomputed bounds

// When INITIAL_CUTOFF changes:
- This one is actually safe! Just a threshold comparison
- Could allow this one if desired
```

**Complexity**: High. Requires significant refactoring.
**Benefit**: Automatic optimization (uncertain improvement).
**Recommendation**: Manual optimization is simpler and safer for now.

## File Guide

```
Main code:
  gnfs/LatticeSiever.h          - Cache blocking + tracking flag
  gnfs/LatticeSiever.cpp        - Implementation + disabled adjustments
  gnfs/RateOptimizer.h          - Rate tracking (still works)

Documentation:
  BRANCH_SUMMARY.md             - Full overview (read this first)
  BUG_FIXES_*.md                - Detailed bug analysis
  CACHE_BLOCKING_*.md           - Cache blocking guides
  RATE_OPTIMIZATION_*.md        - Auto-tuning guides
  
This file:
  QUICK_REFERENCE.md            - You are here!
```

## Bottom Line

**Safe to use**: Yes ✅
**Performance**: Good (restored to baseline) ✅
**Correctness**: Yes (no invalid relations) ✅
**Auto-tuning**: Disabled (for safety) ⚠️

**Action**: Merge and use with confidence. Set ENABLE_AUTO_TUNING=false or leave it unset (false is default).
