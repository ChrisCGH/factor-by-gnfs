# Quick Start: Statistics Output Feature

## TL;DR

Add these 2 lines to your `sieve.cfg` file to get real-time sieving statistics:

```
STATS_FILE = sieve_stats.csv
STATS_INTERVAL = 100
```

Then monitor progress with:
```bash
tail -f sieve_stats.csv
```

## What You Get

A CSV file that updates every 100 relations with:
- Elapsed time
- Total relations found
- Current and average relations/second
- Projected relations/hour and /day

## Example Output

```csv
elapsed_time,total_relations,current_rel_per_sec,running_avg_rel_per_sec,rel_per_hour,rel_per_day
10.5,100,9.52,9.52,34272,822528
25.8,200,7.75,7.75,27900,669600
41.2,300,6.55,7.28,26208,628992
```

## Quick Commands

```bash
# Watch statistics update in real-time
tail -f sieve_stats.csv

# Get latest statistics
tail -n 1 sieve_stats.csv

# Extract current average rate
tail -n 1 sieve_stats.csv | cut -d',' -f4

# Count total relations found
tail -n 1 sieve_stats.csv | cut -d',' -f2
```

## Configuration Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| STATS_FILE | sieve_stats.csv | Output filename |
| STATS_INTERVAL | 100 | Relations between writes |

## Recommended Settings

- **Testing**: `STATS_INTERVAL = 10` (frequent updates)
- **Normal**: `STATS_INTERVAL = 100` (balanced)
- **Long runs**: `STATS_INTERVAL = 500` (less overhead)
- **Disabled**: `STATS_INTERVAL = 0` (no output)

## Visualization

### Python
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('sieve_stats.csv')
df.plot(x='elapsed_time', y='running_avg_rel_per_sec')
plt.title('Sieving Rate Over Time')
plt.xlabel('Time (seconds)')
plt.ylabel('Relations/Second')
plt.show()
```

### Excel/LibreOffice
Just open `sieve_stats.csv` and create charts from the data.

## Troubleshooting

**File not created?**
- Check `STATS_INTERVAL > 0`
- Verify directory is writable

**No data?**
- Wait for N relations (where N = STATS_INTERVAL)
- Check siever is finding relations

**Disable output?**
- Set `STATS_INTERVAL = 0` or omit the parameters

## More Information

See complete documentation:
- **STATS_OUTPUT_GUIDE.md** - Full user guide
- **STATS_CONFIG_EXAMPLE.txt** - Configuration reference
- **STATS_FEATURE_SUMMARY.md** - Implementation details

## That's It!

The feature is fully integrated and ready to use. Just add the configuration and start sieving!
