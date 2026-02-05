# Quick Start: Relation Rate Optimization

## What It Does
Automatically tracks relation generation rates and adjusts sieving parameters to maximize throughput.

## Enable It
Add to your `sieve.cfg`:
```
ENABLE_AUTO_TUNING = true
```

## What You'll See
```
5 -> 123 -> 8 (850,3062,73489), (978,3521,84514), [trend:2.3, avg:975, best:1005]
```
- First numbers: Sieving stages
- First parentheses: Current rate (per sec/hour/day)
- Second parentheses: Running average rate
- Brackets: Trend, moving average, best rate

## When Parameters Adjust
```
### Auto-tuning: Adjusting parameters ###
Rate trend: -0.8
Current rate: 850 rel/sec
Best rate: 1005 rel/sec
Adjustments: B1=500, B2=200, SBA1=1, IC=2
New parameters: B1=1000500, B2=700200, SBA1=12, IC=12
###################################
```

## Adjustable Parameters
- **B1, B2**: Factor base bounds
- **SIEVE_BOUND_ADJUSTMENT1/2**: Sieve thresholds
- **INITIAL_CUTOFF**: Initial cutoff value

All have safety bounds to prevent extreme values.

## Disable It
Remove the line or set:
```
ENABLE_AUTO_TUNING = false
```

## More Info
- **Complete Guide**: See `RATE_OPTIMIZATION_GUIDE.md`
- **Implementation Details**: See `RATE_OPTIMIZATION_IMPLEMENTATION.md`
- **Test Suite**: Run `gnfs/RateOptimizerTest`

## Key Points
✅ Disabled by default (backward compatible)  
✅ Safe (bounded parameters)  
✅ Adaptive (responds to performance changes)  
✅ Informative (detailed metrics in verbose mode)  
✅ Tested (comprehensive unit tests)  

That's it! The system will automatically try to maximize your relation generation rate.
