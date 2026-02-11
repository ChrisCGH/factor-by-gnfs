#!/usr/bin/env python3
"""
Hill-climbing parameter tuner for the GNFS lattice siever.

Reads an existing sieve.cfg, runs lsieve in sampling mode with different
parameter combinations, and uses hill-climbing to find parameters that
maximise the relation generation rate (relations per second).

Usage:
    python3 tune_sieve.py [OPTIONS]

Options:
    --config FILE       Path to sieve.cfg (default: sieve.cfg)
    --lsieve PATH       Path to lsieve binary (default: gnfs/gbin/lsieve)
    --min-q Q           Minimum q for sampling (default: 1000000)
    --max-q Q           Maximum q for sampling (default: 1000500)
    --max-iterations N  Maximum hill-climbing iterations (default: 50)
    --output FILE       Write best config to FILE (default: sieve.cfg)
    --dry-run           Show what would be done without running lsieve
    --verbose           Show detailed output
"""

import argparse
import copy
import os
import re
import subprocess
import sys
import tempfile


# ---------------------------------------------------------------------------
# Configuration parsing
# ---------------------------------------------------------------------------

def parse_sieve_config(path):
    """Parse a sieve.cfg file and return (lines, params).

    *lines*  – the raw lines of the file (for faithful round-tripping).
    *params* – dict mapping parameter names to their string values.
    """
    lines = []
    params = {}
    with open(path, "r") as fh:
        for raw in fh:
            lines.append(raw)
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            eqpos = stripped.find("=")
            if eqpos < 0:
                continue
            key = stripped[:eqpos].strip()
            val = stripped[eqpos + 1:].strip()
            params[key] = val
    return lines, params


def write_sieve_config(lines, params, dest):
    """Write *lines* back to *dest*, replacing any parameter present in *params*."""
    with open(dest, "w") as fh:
        in_polynomial = False
        for raw in lines:
            stripped = raw.strip()
            # Track multi-line polynomial blocks (DEGREE + coefficient lines)
            if in_polynomial:
                if stripped and not stripped.startswith("#"):
                    eqpos = stripped.find("=")
                    # A line with '=' that is a known key ends the polynomial
                    if eqpos >= 0:
                        key = stripped[:eqpos].strip()
                        if key != "DEGREE" and key in params:
                            in_polynomial = False
                            # fall through to normal handling
                        else:
                            fh.write(raw)
                            continue
                    else:
                        fh.write(raw)
                        continue
                else:
                    fh.write(raw)
                    continue

            if not stripped or stripped.startswith("#"):
                fh.write(raw)
                continue

            eqpos = stripped.find("=")
            if eqpos < 0:
                fh.write(raw)
                continue

            key = stripped[:eqpos].strip()
            if key in ("f1", "f2"):
                # Keep polynomial lines as-is; enter polynomial mode if multi-line
                fh.write(raw)
                val = stripped[eqpos + 1:].strip()
                if val == "Polynomial":
                    in_polynomial = True
                continue

            if key in params:
                fh.write(f"{key} = {params[key]}\n")
            else:
                fh.write(raw)


# ---------------------------------------------------------------------------
# Running the siever and measuring the relation rate
# ---------------------------------------------------------------------------

def run_lsieve(lsieve_path, config_path, min_q, max_q):
    """Run lsieve in sampling mode and return the relation rate (rel/sec).

    Returns *None* if the run fails or no rate can be parsed.
    """
    # lsieve reads sieve.cfg from the current working directory, so we run
    # it in the directory containing the config file.
    cfg_dir = os.path.dirname(os.path.abspath(config_path)) or "."
    lsieve_abs = os.path.abspath(lsieve_path)

    env = dict(os.environ)
    env["LATTICE_SIEVER_VERBOSE"] = "1"

    try:
        proc = subprocess.run(
            [lsieve_abs, "-s", str(min_q), str(max_q)],
            cwd=cfg_dir,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=600,
        )
    except FileNotFoundError:
        print(f"Error: lsieve binary not found at {lsieve_abs}", file=sys.stderr)
        return None
    except subprocess.TimeoutExpired:
        print("Warning: lsieve timed out", file=sys.stderr)
        return None

    output = proc.stdout.decode("utf-8", errors="replace")
    stderr_output = proc.stderr.decode("utf-8", errors="replace")

    # Parse the running-average relation rate from the verbose output.
    # Output format per q:
    #   <n1> -> <n2> -> <rels> (<inst_rate>,<hr>,<day>), (<run_rate>,<hr>,<day>),
    # We want the *last* running-average rate printed.
    rate = _parse_relation_rate(output)

    # Fallback: try to parse the stats CSV if present.
    if rate is None:
        stats_path = os.path.join(cfg_dir, "sieve_stats.csv")
        rate = _parse_stats_csv(stats_path)

    if rate is None and stderr_output:
        # Check stderr for rate information
        rate = _parse_relation_rate(stderr_output)

    return rate


def _parse_relation_rate(text):
    """Extract the last running-average relation rate from lsieve output.

    The verbose output contains entries like:
        ... (123.4,444960,10679040), (456.7,1644120,39458880),
    The second parenthesised group is the running average.  We want the
    first number (relations per second) from the *last* such group.
    """
    # Match pairs of parenthesised rate groups on a line.
    # The running average is the second group.
    pattern = re.compile(
        r"\([\d.e+-]+,[\d.e+-]+,[\d.e+-]+\)"
        r",?\s*"
        r"\(([\d.e+-]+),[\d.e+-]+,[\d.e+-]+\)"
    )
    last_rate = None
    for m in pattern.finditer(text):
        try:
            last_rate = float(m.group(1))
        except ValueError:
            pass
    return last_rate


def _parse_stats_csv(path):
    """Parse the last line of a sieve_stats.csv and return running_avg rel/sec."""
    try:
        with open(path, "r") as fh:
            lines = fh.readlines()
        for line in reversed(lines):
            parts = line.strip().split(",")
            if len(parts) >= 4:
                try:
                    return float(parts[3])  # running_avg_rel_per_sec
                except ValueError:
                    continue
    except FileNotFoundError:
        pass
    return None


# ---------------------------------------------------------------------------
# Hill-climbing optimiser
# ---------------------------------------------------------------------------

# The parameters we tune and their step sizes.
TUNABLE_PARAMS = [
    # (config key, step size, min value, max value)
    ("SIEVE_BOUND_ADJUSTMENT1", 1, 0, 40),
    ("SIEVE_BOUND_ADJUSTMENT2", 1, 0, 40),
    ("SMALL_PRIME_BOUND1", 100, 0, 5000),
    ("SMALL_PRIME_BOUND2", 100, 0, 5000),
    ("INITIAL_CUTOFF", 5, 0, 400),
]


def _set_param(params, key, value):
    """Set an integer parameter in *params*, clamping to allowed range."""
    for pkey, _step, lo, hi in TUNABLE_PARAMS:
        if pkey == key:
            value = max(lo, min(hi, value))
            break
    params[key] = str(int(value))


def _get_param(params, key):
    """Get an integer parameter from *params*."""
    return int(params.get(key, "0"))


def evaluate(params, lines, config_path, lsieve_path, min_q, max_q):
    """Write *params* to a temporary config and measure the relation rate."""
    cfg_dir = os.path.dirname(os.path.abspath(config_path)) or "."
    tmp_cfg = os.path.join(cfg_dir, "sieve.cfg")

    # Write config
    write_sieve_config(lines, params, tmp_cfg)

    rate = run_lsieve(lsieve_path, tmp_cfg, min_q, max_q)
    return rate


def hill_climb(lines, params, config_path, lsieve_path, min_q, max_q,
               max_iterations, verbose=False):
    """Perform hill-climbing over the tunable sieve parameters.

    Returns (best_params, best_rate).
    """
    current_params = copy.deepcopy(params)

    print("=== Measuring baseline rate ===")
    best_rate = evaluate(current_params, lines, config_path,
                         lsieve_path, min_q, max_q)
    if best_rate is None:
        print("Error: could not measure baseline rate. "
              "Check that lsieve runs correctly.", file=sys.stderr)
        return current_params, None

    best_params = copy.deepcopy(current_params)
    print(f"Baseline rate: {best_rate:.2f} rel/sec")
    print(f"Baseline parameters:")
    for key, _step, _lo, _hi in TUNABLE_PARAMS:
        print(f"  {key} = {_get_param(current_params, key)}")
    print()

    for iteration in range(1, max_iterations + 1):
        improved = False
        for key, step, lo, hi in TUNABLE_PARAMS:
            current_val = _get_param(current_params, key)
            for direction in (+1, -1):
                candidate_params = copy.deepcopy(current_params)
                new_val = current_val + direction * step
                if new_val < lo or new_val > hi:
                    continue
                _set_param(candidate_params, key, new_val)

                if verbose:
                    print(f"  Iteration {iteration}: trying {key} = {new_val} "
                          f"(was {current_val}) ...", end=" ", flush=True)

                rate = evaluate(candidate_params, lines, config_path,
                                lsieve_path, min_q, max_q)

                if rate is None:
                    if verbose:
                        print("FAILED")
                    continue

                if verbose:
                    print(f"{rate:.2f} rel/sec", end="")

                if rate > best_rate:
                    best_rate = rate
                    best_params = copy.deepcopy(candidate_params)
                    current_params = copy.deepcopy(candidate_params)
                    improved = True
                    if verbose:
                        print(f"  ** NEW BEST **")
                    else:
                        print(f"  Iteration {iteration}: {key} = {new_val} -> "
                              f"{rate:.2f} rel/sec  ** NEW BEST **")
                    break  # restart inner loop with new current params
                else:
                    if verbose:
                        print()

            if improved:
                break  # restart parameter scan from first param

        if not improved:
            print(f"  Iteration {iteration}: no improvement found – stopping.")
            break

    return best_params, best_rate


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Hill-climbing parameter tuner for the GNFS lattice siever."
    )
    parser.add_argument(
        "--config", default="sieve.cfg",
        help="Path to sieve.cfg (default: sieve.cfg)"
    )
    parser.add_argument(
        "--lsieve", default="gnfs/gbin/lsieve",
        help="Path to lsieve binary (default: gnfs/gbin/lsieve)"
    )
    parser.add_argument(
        "--min-q", type=int, default=1000000,
        help="Minimum q for sampling (default: 1000000)"
    )
    parser.add_argument(
        "--max-q", type=int, default=1000500,
        help="Maximum q for sampling (default: 1000500)"
    )
    parser.add_argument(
        "--max-iterations", type=int, default=50,
        help="Maximum hill-climbing iterations (default: 50)"
    )
    parser.add_argument(
        "--output", default=None,
        help="Write best config to FILE (default: same as --config)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show what would be done without running lsieve"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Show detailed output"
    )
    args = parser.parse_args()

    config_path = args.config
    output_path = args.output or config_path

    if not os.path.isfile(config_path):
        print(f"Error: config file not found: {config_path}", file=sys.stderr)
        return 1

    if not args.dry_run and not os.path.isfile(args.lsieve):
        print(f"Error: lsieve binary not found: {args.lsieve}", file=sys.stderr)
        print("Build it first:  cd gnfs && make lsieve", file=sys.stderr)
        return 1

    lines, params = parse_sieve_config(config_path)

    print("=" * 60)
    print("GNFS Lattice Siever Parameter Tuner")
    print("=" * 60)
    print(f"Config file : {config_path}")
    print(f"lsieve      : {args.lsieve}")
    print(f"q range     : {args.min_q} – {args.max_q}")
    print(f"Max iters   : {args.max_iterations}")
    print()

    if args.dry_run:
        print("Dry-run mode: showing tunable parameters and their ranges.")
        print()
        for key, step, lo, hi in TUNABLE_PARAMS:
            cur = params.get(key, "(not set)")
            print(f"  {key:30s} = {cur:>6s}  (step={step}, range={lo}..{hi})")
        print()
        print("Re-run without --dry-run to start tuning.")
        return 0

    best_params, best_rate = hill_climb(
        lines, params, config_path, args.lsieve,
        args.min_q, args.max_q, args.max_iterations,
        verbose=args.verbose,
    )

    print()
    print("=" * 60)
    if best_rate is not None:
        print(f"Best rate achieved: {best_rate:.2f} rel/sec")
    else:
        print("Could not determine a relation rate.")
    print("Best parameters:")
    for key, _step, _lo, _hi in TUNABLE_PARAMS:
        print(f"  {key} = {_get_param(best_params, key)}")
    print("=" * 60)

    # Write the best configuration
    write_sieve_config(lines, best_params, output_path)
    print(f"Best configuration written to {output_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
