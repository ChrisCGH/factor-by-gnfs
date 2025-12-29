#!/bin/bash
# Performance analysis helper for LatticeSiever

set -e

GNFS_DIR="/home/runner/work/factor-by-gnfs/factor-by-gnfs/gnfs"
cd "$GNFS_DIR"

echo "=== LatticeSiever Performance Analysis ==="
echo ""

# Check if lsieve exists
if [ ! -f "gbin/lsieve" ]; then
    echo "Building lsieve..."
    make lsieve 2>&1 | tail -5
    echo ""
fi

# Default test range - adjust based on your needs
MIN_Q=${1:-1000000}
MAX_Q=${2:-1000010}

echo "Running test sieve from q=$MIN_Q to q=$MAX_Q..."
echo ""

# Clear old timing data
rm -f sieve.tim

# Run the siever
./gbin/lsieve $MIN_Q $MAX_Q 2>&1 | tee /tmp/sieve_output.txt

echo ""
echo "=== Built-in Timer Results ==="
if [ -f "sieve.tim" ]; then
    cat sieve.tim
    echo ""
    echo "=== Top Time Consumers ==="
    grep -v "^$" sieve.tim | sort -k2 -nr | head -10
else
    echo "Warning: sieve.tim not found"
fi

echo ""
echo "=== Performance Summary from Output ==="
grep "relations.*second\|relations/sec" /tmp/sieve_output.txt || echo "No performance data in output"

echo ""
echo "=== Analysis Complete ==="
echo ""
echo "To run detailed profiling:"
echo "  1. perf:      perf record -g ./gbin/lsieve $MIN_Q $MAX_Q && perf report"
echo "  2. callgrind: valgrind --tool=callgrind ./gbin/lsieve $MIN_Q $MAX_Q"
echo "  3. gprof:     Uncomment PROFILE=-pg in makefile, rebuild, run, then: gprof ./gbin/lsieve"
echo ""
echo "See PROFILING_GUIDE.md for detailed instructions"
