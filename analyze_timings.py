#!/usr/bin/env python3
"""
Analyze LatticeSiever performance data and suggest optimizations
"""

import sys
import re
import subprocess
from pathlib import Path

def parse_timer_file(timer_file):
    """Parse sieve.tim file and return timing data"""
    timings = {}
    try:
        with open(timer_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or ':' not in line:
                    continue
                parts = line.split(':')
                if len(parts) >= 2:
                    name = parts[0].strip()
                    time_str = parts[1].strip().split()[0]
                    try:
                        time_val = float(time_str)
                        timings[name] = time_val
                    except ValueError:
                        pass
    except FileNotFoundError:
        print(f"Warning: {timer_file} not found", file=sys.stderr)
    return timings

def analyze_timings(timings):
    """Analyze timing data and provide recommendations"""
    if not timings:
        return []
    
    total_time = sum(timings.values())
    if total_time == 0:
        return []
    
    # Calculate percentages
    timing_pcts = {name: (time/total_time)*100 for name, time in timings.items()}
    
    # Sort by time
    sorted_timings = sorted(timing_pcts.items(), key=lambda x: x[1], reverse=True)
    
    recommendations = []
    
    for name, pct in sorted_timings[:5]:  # Top 5 consumers
        time = timings[name]
        rec = {
            'phase': name,
            'time': time,
            'percent': pct,
            'suggestions': []
        }
        
        # Phase-specific optimization suggestions
        if 'sieve by vectors' in name:
            if pct > 30:
                rec['suggestions'].extend([
                    "Memory-bound: Consider cache-blocking techniques",
                    "Try reducing the sieve region size to fit in L3 cache",
                    "Consider SIMD vectorization for array operations",
                    "Profile memory access patterns with 'perf mem record'"
                ])
        
        elif 'check interval' in name:
            if pct > 20:
                rec['suggestions'].extend([
                    "Evaluate speedup: Try different smoothness bounds",
                    "Consider early-exit strategies for non-smooth candidates",
                    "Profile branch mispredictions: perf stat -e branch-misses",
                    "Optimize evaluate_on_lattice() function"
                ])
        
        elif 'remove factors' in name:
            if pct > 15:
                rec['suggestions'].extend([
                    "Trial division bottleneck: Our cached primes optimization should help",
                    "Consider wheel factorization for small primes",
                    "Try Montgomery multiplication for modular arithmetic",
                    "Use GMP's mpz_divisible_ui_p for faster divisibility tests"
                ])
        
        elif 'eliminate' in name:
            if pct > 10:
                rec['suggestions'].extend([
                    "Optimize factor checking logic",
                    "Consider parallel processing of candidates",
                    "Reduce unnecessary memory allocations"
                ])
        
        if pct > 10 and not rec['suggestions']:
            rec['suggestions'].append(
                f"This phase takes {pct:.1f}% of time - profile with 'perf record -g'"
            )
        
        recommendations.append(rec)
    
    return recommendations

def get_cache_info():
    """Try to get cache size information"""
    try:
        result = subprocess.run(['lscpu'], capture_output=True, text=True)
        cache_info = []
        for line in result.stdout.split('\n'):
            if 'cache' in line.lower():
                cache_info.append(line.strip())
        return cache_info
    except:
        return []

def main():
    timer_file = sys.argv[1] if len(sys.argv) > 1 else 'gnfs/sieve.tim'
    
    print("=" * 70)
    print("LatticeSiever Performance Analysis")
    print("=" * 70)
    print()
    
    # Parse timings
    timings = parse_timer_file(timer_file)
    
    if not timings:
        print("No timing data found. Run the siever first:")
        print("  cd gnfs && ./gbin/lsieve <min_q> <max_q>")
        return 1
    
    # Display raw timings
    total_time = sum(timings.values())
    print(f"Total Time: {total_time:.3f}s")
    print()
    print("Time Distribution:")
    print("-" * 70)
    
    for name, time in sorted(timings.items(), key=lambda x: x[1], reverse=True):
        pct = (time/total_time)*100
        bar = '█' * int(pct / 2)
        print(f"{name:30s} {time:8.3f}s  {pct:5.1f}%  {bar}")
    print()
    
    # Analyze and recommend
    recommendations = analyze_timings(timings)
    
    if recommendations:
        print("=" * 70)
        print("Optimization Recommendations")
        print("=" * 70)
        print()
        
        for i, rec in enumerate(recommendations, 1):
            print(f"{i}. {rec['phase']} - {rec['percent']:.1f}% of total time")
            if rec['suggestions']:
                for suggestion in rec['suggestions']:
                    print(f"   • {suggestion}")
            print()
    
    # Cache information
    cache_info = get_cache_info()
    if cache_info:
        print("=" * 70)
        print("System Cache Information")
        print("=" * 70)
        for info in cache_info:
            print(f"  {info}")
        print()
        print("Note: Sieve array size should fit in L3 cache for best performance")
        print()
    
    print("=" * 70)
    print("Next Steps")
    print("=" * 70)
    print("1. Focus on phases taking >20% of total time")
    print("2. Run: perf record -g ./gbin/lsieve <args> && perf report")
    print("3. Check cache misses: perf stat -e cache-misses,cache-references ./gbin/lsieve <args>")
    print("4. See PROFILING_GUIDE.md for detailed profiling instructions")
    print()
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
