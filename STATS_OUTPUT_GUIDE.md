# Statistics Output Guide for Lattice Siever

## Overview

The lattice siever can now output statistics on relation generation rate to a configurable file. This allows real-time monitoring of sieving performance and data collection for analysis.

## Configuration

Add these parameters to your `sieve.cfg` file:

```
# Statistics output file (default: sieve_stats.csv)
STATS_FILE = sieve_stats.csv

# Write statistics every N relations (default: 100)
# Set to 0 to disable statistics output
STATS_INTERVAL = 100
```

### Configuration Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `STATS_FILE` | string | "sieve_stats.csv" | Path to output CSV file |
| `STATS_INTERVAL` | integer | 100 | Number of relations between statistics writes |

## Output Format

The statistics file is a CSV (Comma-Separated Values) file with the following columns:

| Column | Description | Units |
|--------|-------------|-------|
| `elapsed_time` | Total sieving time since start | seconds |
| `total_relations` | Cumulative number of relations found | count |
| `current_rel_per_sec` | Relations per second for current special-q | relations/sec |
| `running_avg_rel_per_sec` | Overall average relations per second | relations/sec |
| `rel_per_hour` | Running average × 3600 | relations/hour |
| `rel_per_day` | Running average × 86400 | relations/day |

### Example Output

```csv
elapsed_time,total_relations,current_rel_per_sec,running_avg_rel_per_sec,rel_per_hour,rel_per_day
10.5,100,9.52,9.52,34272,822528
25.8,200,7.75,7.75,27900,669600
41.2,300,6.55,7.28,26208,628992
57.0,400,6.33,7.02,25272,606528
```

## Usage Examples

### Basic Usage

1. Add configuration to `sieve.cfg`:
   ```
   STATS_FILE = my_stats.csv
   STATS_INTERVAL = 50
   ```

2. Run the siever:
   ```bash
   ./lsieve sieve.cfg
   ```

3. Monitor progress:
   ```bash
   tail -f my_stats.csv
   ```

### Disabling Statistics Output

To disable statistics output, either:

**Option 1:** Set interval to 0
```
STATS_INTERVAL = 0
```

**Option 2:** Omit the configuration parameters entirely

### Different Output Frequencies

**Every 10 relations (frequent updates):**
```
STATS_INTERVAL = 10
```

**Every 500 relations (less frequent):**
```
STATS_INTERVAL = 500
```

**Every 1000 relations (minimal overhead):**
```
STATS_INTERVAL = 1000
```

## Analysis Examples

### Using Python/pandas

```python
import pandas as pd
import matplotlib.pyplot as plt

# Read statistics
df = pd.read_csv('sieve_stats.csv')

# Plot relations over time
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(df['elapsed_time'], df['total_relations'])
plt.xlabel('Elapsed Time (seconds)')
plt.ylabel('Total Relations')
plt.title('Relations Found Over Time')

# Plot rate over time
plt.subplot(1, 2, 2)
plt.plot(df['elapsed_time'], df['running_avg_rel_per_sec'])
plt.xlabel('Elapsed Time (seconds)')
plt.ylabel('Relations per Second')
plt.title('Sieving Rate Over Time')

plt.tight_layout()
plt.savefig('sieving_performance.png')
plt.show()
```

### Using Unix Tools

**Calculate average rate:**
```bash
tail -n 1 sieve_stats.csv | cut -d',' -f4
```

**Count total relations:**
```bash
tail -n 1 sieve_stats.csv | cut -d',' -f2
```

**Extract recent performance:**
```bash
tail -n 10 sieve_stats.csv
```

## Implementation Details

### When Statistics Are Written

- Statistics are written after every N relations (where N = `STATS_INTERVAL`)
- The counter resets after each write
- Example with `STATS_INTERVAL = 100`:
  - Write after 100 relations (cumulative: 100)
  - Write after next 100 relations (cumulative: 200)
  - Write after next 100 relations (cumulative: 300)
  - And so on...

### File Handling

- File is created when siever starts (if `STATS_INTERVAL > 0`)
- Header row is written immediately
- Data is appended incrementally
- File is flushed after each write for real-time monitoring
- File is closed when siever exits

### Performance Impact

The performance impact of statistics output is minimal:
- File I/O only occurs every N relations (not continuously)
- Simple CSV write with 6 numeric values
- File is flushed but not re-opened on each write
- Typical overhead: < 0.1% with `STATS_INTERVAL = 100`

### Recommended Settings

| Scenario | Recommended Interval | Rationale |
|----------|---------------------|-----------|
| Quick testing | 10-50 | Frequent updates for monitoring |
| Normal sieving | 100-200 | Good balance of detail and overhead |
| Long runs | 500-1000 | Minimize file size and overhead |
| No monitoring | 0 | Disable to eliminate overhead |

## Troubleshooting

### File Not Created

**Issue:** Statistics file doesn't appear

**Solutions:**
1. Check that `STATS_INTERVAL > 0` in config
2. Verify file path is valid and writable
3. Check permissions in output directory
4. Ensure siever has started and found relations

### No Data Written

**Issue:** File created but empty (except header)

**Solutions:**
1. Check that relations are being found (console output)
2. Verify `STATS_INTERVAL` is reasonable (not too large)
3. Wait for N relations to be found before first write
4. Check console for error messages

### Unexpected Values

**Issue:** Statistics values seem incorrect

**Solutions:**
1. Verify relations are being found correctly
2. Check that timer is working (elapsed_time increasing)
3. Compare with console verbose output
4. Ensure siever wasn't restarted mid-run

## Advanced Usage

### Multiple Runs

To track statistics across multiple siever runs:

1. Use unique filenames:
   ```
   STATS_FILE = stats_run1.csv
   ```

2. Or append timestamp:
   ```bash
   STATS_FILE=sieve_stats_$(date +%Y%m%d_%H%M%S).csv ./lsieve sieve.cfg
   ```

### Real-Time Monitoring Script

```bash
#!/bin/bash
# monitor_sieving.sh - Display live statistics

STATS_FILE="sieve_stats.csv"

echo "Monitoring sieving performance..."
echo "Press Ctrl+C to stop"
echo

while true; do
    if [ -f "$STATS_FILE" ]; then
        clear
        echo "=== Latest Sieving Statistics ==="
        echo
        tail -n 1 "$STATS_FILE" | awk -F',' '{
            printf "Elapsed Time: %.1f seconds\n", $1
            printf "Total Relations: %d\n", $2
            printf "Current Rate: %.2f rel/sec\n", $3
            printf "Average Rate: %.2f rel/sec\n", $4
            printf "Per Hour: %d relations\n", int($5)
            printf "Per Day: %d relations\n", int($6)
        }'
    else
        echo "Waiting for statistics file..."
    fi
    sleep 2
done
```

### Gnuplot Visualization

```gnuplot
# plot_stats.gp - Generate performance graphs

set datafile separator ","
set terminal png size 1200,800
set output 'sieving_performance.png'

set multiplot layout 2,2

# Plot 1: Total relations over time
set title "Total Relations Found"
set xlabel "Elapsed Time (seconds)"
set ylabel "Total Relations"
plot 'sieve_stats.csv' using 1:2 with lines title "Relations"

# Plot 2: Current rate
set title "Current Sieving Rate"
set xlabel "Elapsed Time (seconds)"
set ylabel "Relations per Second"
plot 'sieve_stats.csv' using 1:3 with lines title "Current Rate"

# Plot 3: Running average rate
set title "Running Average Rate"
set xlabel "Elapsed Time (seconds)"
set ylabel "Relations per Second"
plot 'sieve_stats.csv' using 1:4 with lines title "Avg Rate"

# Plot 4: Projected daily rate
set title "Projected Daily Yield"
set xlabel "Elapsed Time (seconds)"
set ylabel "Relations per Day"
plot 'sieve_stats.csv' using 1:6 with lines title "Per Day"

unset multiplot
```

Run with: `gnuplot plot_stats.gp`

## See Also

- Main lattice siever documentation
- Configuration file format (`sieve.cfg`)
- Performance tuning guide
- Rate optimization features
