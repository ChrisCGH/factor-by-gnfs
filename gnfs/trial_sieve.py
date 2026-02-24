#!/usr/bin/env python3
"""
Trial sieving utility for GNFS polynomial selection.

Reads output from the skewed polynomial selection program, generates
sieve.cfg files from a template, and runs trial sieving with lsieve
to find the polynomials that give the best yield.

The skewed program produces output blocks in this format:

    f1 = <polynomial>
    m = <m>
    a = <a>
    b = <b>
    s = <s>
    alpha = <alpha>
    E(F) = <E_F>
    N = <N>

    ==============================================================================

Each block is separated by a line of '=' characters.

A template file is used to generate sieve.cfg files.  Placeholders in
the template are replaced with values from the skewed output:

    <N>      - the number being factored
    <m>      - the common root
    <f1>     - the algebraic polynomial
    <a>      - the 'a' value (used in f2)
    <b>      - the 'b' value (used in f2)
    <s>      - the skewedness
    <alpha>  - the alpha value
    <E_F>    - the E(F) value

Usage:
    python3 trial_sieve.py [OPTIONS]

Options:
    --skewed FILE       Path to skewed output file (default: skewed.out)
    --template FILE     Path to template sieve.cfg (default: sieve.cfg.template)
    --lsieve PATH       Path to lsieve binary (default: gbin/lsieve)
    --min-q Q           Minimum q for trial sieving (default: 1000000)
    --max-q Q           Maximum q for trial sieving (default: 1010000)
    --timeout SECS      Timeout per lsieve run in seconds (default: 600)
    --workdir DIR       Working directory for sieve runs (default: trial_sieve_work)
    --dry-run           Parse and generate configs without running lsieve
    --verbose           Show detailed output
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile


# ---------------------------------------------------------------------------
# Parsing skewed output
# ---------------------------------------------------------------------------

def parse_skewed_output(path):
    """Parse a skewed polynomial selection output file.

    Returns a list of dicts, each containing the fields from one
    polynomial block: f1, m, a, b, s, alpha, E_F, N.
    Not all fields are required to be present in every block.
    """
    blocks = []
    current = {}

    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()

            # Block separator
            if line and all(c == '=' for c in line):
                if current:
                    blocks.append(current)
                    current = {}
                continue

            # Skip empty lines
            if not line:
                continue

            # Parse "key = value" lines
            if line.startswith("f1 = "):
                current["f1"] = line[len("f1 = "):]
            elif line.startswith("m = "):
                current["m"] = line[len("m = "):]
            elif line.startswith("a = "):
                current["a"] = line[len("a = "):]
            elif line.startswith("b = "):
                current["b"] = line[len("b = "):]
            elif line.startswith("s = "):
                current["s"] = line[len("s = "):]
            elif line.startswith("alpha = "):
                current["alpha"] = line[len("alpha = "):]
            elif line.startswith("E(F) = "):
                current["E_F"] = line[len("E(F) = "):]
            elif line.startswith("N = "):
                current["N"] = line[len("N = "):]
            elif "f1" not in current and not line.startswith("#"):
                # First non-key line before f1 is set: could be the
                # polynomial without the "f1 = " prefix (older format)
                if re.search(r'X(\^\d+)?', line) or re.match(r'^-?\d+\s', line):
                    current["f1"] = line

    # Don't forget the last block if not terminated by separator
    if current:
        blocks.append(current)

    return blocks


# ---------------------------------------------------------------------------
# Template processing
# ---------------------------------------------------------------------------

def read_template(path):
    """Read a template file and return its contents as a string."""
    with open(path, "r") as fh:
        return fh.read()


def generate_config(template, block):
    """Generate a sieve.cfg from a template by substituting placeholders.

    Supported placeholders: <N>, <m>, <f1>, <a>, <b>, <s>, <alpha>, <E_F>.
    """
    result = template
    replacements = {
        "<N>": block.get("N", ""),
        "<m>": block.get("m", ""),
        "<f1>": block.get("f1", ""),
        "<a>": block.get("a", ""),
        "<b>": block.get("b", ""),
        "<s>": block.get("s", ""),
        "<alpha>": block.get("alpha", ""),
        "<E_F>": block.get("E_F", ""),
    }
    for placeholder, value in replacements.items():
        result = result.replace(placeholder, value)
    return result


# ---------------------------------------------------------------------------
# Running trial sieving
# ---------------------------------------------------------------------------

def run_trial_sieve(lsieve_path, config_content, min_q, max_q,
                    workdir, block_index, timeout=600,
                    sample_count=None):
    """Run lsieve in sampling mode for one polynomial.

    Creates a temporary sieve.cfg in a subdirectory of workdir,
    runs lsieve, and returns the number of relations found.
    """
    block_dir = os.path.join(workdir, f"poly_{block_index}")
    os.makedirs(block_dir, exist_ok=True)

    config_path = os.path.join(block_dir, "sieve.cfg")
    with open(config_path, "w") as fh:
        fh.write(config_content)

    lsieve_abs = os.path.abspath(lsieve_path)

    cmd = [lsieve_abs, "-s"]
    if sample_count is not None:
        cmd.extend(["-n", str(sample_count)])
    cmd.extend([str(min_q), str(max_q)])

    try:
        proc = subprocess.run(
            cmd,
            cwd=block_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
    except FileNotFoundError:
        print(f"Error: lsieve binary not found at {lsieve_abs}",
              file=sys.stderr)
        return None
    except subprocess.TimeoutExpired:
        print(f"Warning: lsieve timed out for polynomial {block_index}",
              file=sys.stderr)
        return None

    output = proc.stdout.decode("utf-8", errors="replace")
    stderr_output = proc.stderr.decode("utf-8", errors="replace")

    # Count relations from output and relation file
    relations = count_relations(block_dir, output, stderr_output)
    return relations


def count_relations(block_dir, stdout_text, stderr_text):
    """Count the number of relations found during a trial sieve run.

    Tries multiple strategies:
    1. Count lines in the relation file
    2. Parse relation counts from lsieve output
    """
    # Strategy 1: look for relation files in the block directory
    for fname in os.listdir(block_dir):
        if fname.endswith(".relations"):
            rel_path = os.path.join(block_dir, fname)
            try:
                with open(rel_path, "r") as fh:
                    count = sum(1 for line in fh
                                if line.strip() and not line.startswith("#"))
                if count > 0:
                    return count
            except (IOError, OSError):
                pass

    # Strategy 2: parse the output for relation counts
    # lsieve output format: "<n1> -> <n2> -> <rels> ..."
    total = 0
    pattern = re.compile(r'(\d+)\s*->\s*(\d+)\s*->\s*(\d+)')
    for text in (stdout_text, stderr_text):
        for m in pattern.finditer(text):
            try:
                total += int(m.group(3))
            except ValueError:
                pass
    if total > 0:
        return total

    return 0


# ---------------------------------------------------------------------------
# Summary and reporting
# ---------------------------------------------------------------------------

def format_block_summary(index, block, relations):
    """Format a one-line summary for a polynomial block."""
    e_f = block.get("E_F", "?")
    alpha = block.get("alpha", "?")
    s = block.get("s", "?")
    rels = relations if relations is not None else "ERROR"
    return (f"Poly {index:3d}: relations={rels:>6}, "
            f"E(F)={e_f:>10}, alpha={alpha:>10}, s={s}")


def print_results(results, verbose=False):
    """Print trial sieving results sorted by yield (relations found)."""
    # Sort by relations found, descending
    ranked = sorted(results,
                    key=lambda r: r["relations"] if r["relations"] is not None else -1,
                    reverse=True)

    print()
    print("=" * 70)
    print("Trial Sieving Results (ranked by yield)")
    print("=" * 70)
    for rank, r in enumerate(ranked, 1):
        idx = r["index"]
        block = r["block"]
        rels = r["relations"]
        e_f = block.get("E_F", "?")
        alpha = block.get("alpha", "?")
        s = block.get("s", "?")
        rels_str = str(rels) if rels is not None else "ERROR"
        print(f"  #{rank:2d}  Poly {idx:3d}: relations={rels_str:>6}, "
              f"E(F)={e_f:>10}, alpha={alpha:>10}, s={s}")
        if verbose and "f1" in block:
            print(f"        f1 = {block['f1']}")
    print("=" * 70)

    if ranked and ranked[0]["relations"] is not None:
        best = ranked[0]
        print(f"\nBest polynomial: #{best['index']}")
        print(f"  Relations found: {best['relations']}")
        for key in ("f1", "m", "a", "b", "s", "alpha", "E_F", "N"):
            if key in best["block"]:
                label = "E(F)" if key == "E_F" else key
                print(f"  {label} = {best['block'][key]}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Trial sieving utility for GNFS polynomial selection."
    )
    parser.add_argument(
        "--skewed", default="skewed.out",
        help="Path to skewed output file (default: skewed.out)"
    )
    parser.add_argument(
        "--template", default="sieve.cfg.template",
        help="Path to template sieve.cfg (default: sieve.cfg.template)"
    )
    parser.add_argument(
        "--lsieve", default="gbin/lsieve",
        help="Path to lsieve binary (default: gbin/lsieve)"
    )
    parser.add_argument(
        "--min-q", type=int, default=1000000,
        help="Minimum q for trial sieving (default: 1000000)"
    )
    parser.add_argument(
        "--max-q", type=int, default=1010000,
        help="Maximum q for trial sieving (default: 1010000)"
    )
    parser.add_argument(
        "--timeout", type=int, default=600,
        help="Timeout per lsieve run in seconds (default: 600)"
    )
    parser.add_argument(
        "--workdir", default="trial_sieve_work",
        help="Working directory for sieve runs (default: trial_sieve_work)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Parse and generate configs without running lsieve"
    )
    parser.add_argument(
        "--sample-count", type=int, default=None,
        help="Number of sample points for lsieve -s (default: 100)"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Show detailed output"
    )
    args = parser.parse_args()

    if not os.path.isfile(args.skewed):
        print(f"Error: skewed output file not found: {args.skewed}",
              file=sys.stderr)
        return 1

    if not os.path.isfile(args.template):
        print(f"Error: template file not found: {args.template}",
              file=sys.stderr)
        return 1

    if not args.dry_run and not os.path.isfile(args.lsieve):
        print(f"Error: lsieve binary not found: {args.lsieve}",
              file=sys.stderr)
        print("Build it first:  cd gnfs && make lsieve", file=sys.stderr)
        return 1

    # Parse skewed output
    blocks = parse_skewed_output(args.skewed)
    if not blocks:
        print("Error: no polynomial blocks found in skewed output",
              file=sys.stderr)
        return 1

    template = read_template(args.template)

    print("=" * 70)
    print("GNFS Trial Sieving")
    print("=" * 70)
    print(f"Skewed output : {args.skewed}")
    print(f"Template      : {args.template}")
    print(f"Polynomials   : {len(blocks)}")
    if not args.dry_run:
        print(f"lsieve        : {args.lsieve}")
        print(f"q range       : {args.min_q} - {args.max_q}")
        print(f"Working dir   : {args.workdir}")
    print()

    os.makedirs(args.workdir, exist_ok=True)

    results = []
    for i, block in enumerate(blocks, 1):
        config = generate_config(template, block)

        if args.dry_run:
            config_path = os.path.join(args.workdir, f"poly_{i}", "sieve.cfg")
            os.makedirs(os.path.dirname(config_path), exist_ok=True)
            with open(config_path, "w") as fh:
                fh.write(config)
            print(f"Poly {i}: generated {config_path}")
            if args.verbose:
                for key in ("f1", "m", "a", "b", "s", "alpha", "E_F", "N"):
                    if key in block:
                        label = "E(F)" if key == "E_F" else key
                        print(f"  {label} = {block[key]}")
            continue

        print(f"Trial sieving polynomial {i}/{len(blocks)} ...", flush=True)
        relations = run_trial_sieve(
            args.lsieve, config, args.min_q, args.max_q,
            args.workdir, i, timeout=args.timeout,
            sample_count=args.sample_count
        )
        result = {"index": i, "block": block, "relations": relations}
        results.append(result)
        print(format_block_summary(i, block, relations))

    if not args.dry_run and results:
        print_results(results, verbose=args.verbose)

    return 0


if __name__ == "__main__":
    sys.exit(main())
