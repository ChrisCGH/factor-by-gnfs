#!/usr/bin/env python3
"""Validate sieve.cfg files used by the GNFS lattice siever."""

import argparse
import math
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass
class PolynomialSpec:
    coeffs: List[int]  # constant term first
    source: str

    @property
    def degree(self) -> int:
        return len(self.coeffs) - 1


def parse_sieve_cfg(path: str) -> Tuple[Dict[str, str], Dict[str, PolynomialSpec]]:
    params: Dict[str, str] = {}
    polys: Dict[str, PolynomialSpec] = {}

    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        i += 1
        if not stripped or stripped.startswith("#"):
            continue
        if "=" not in stripped:
            continue

        key, val = stripped.split("=", 1)
        key = key.strip()
        val = val.strip()

        if key in ("f1", "f2") and val == "Polynomial":
            degree, coeffs, i = _parse_multiline_polynomial(lines, i)
            if degree != len(coeffs) - 1:
                raise ValueError(f"{key}: declared DEGREE={degree}, read {len(coeffs)} coefficients")
            polys[key] = PolynomialSpec(coeffs=coeffs, source="multiline")
            params[key] = val
            continue

        if key in ("f1", "f2"):
            polys[key] = PolynomialSpec(coeffs=_parse_inline_polynomial(val), source="inline")

        params[key] = val

    return params, polys


def _parse_multiline_polynomial(lines: List[str], start_index: int) -> Tuple[int, List[int], int]:
    i = start_index

    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError("Expected DEGREE = <n> after polynomial declaration")
        key, val = line.split("=", 1)
        if key.strip() != "DEGREE":
            raise ValueError("Expected DEGREE = <n> after polynomial declaration")
        try:
            degree = int(val.strip())
        except ValueError as exc:
            raise ValueError(f"Invalid DEGREE value: {val.strip()}") from exc
        if degree < 0:
            raise ValueError("DEGREE must be non-negative")
        break
    else:
        raise ValueError("Missing DEGREE after polynomial declaration")

    coeffs: List[int] = []
    while i < len(lines) and len(coeffs) < degree + 1:
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            raise ValueError("Encountered next key before reading all polynomial coefficients")
        try:
            coeffs.append(int(line))
        except ValueError as exc:
            raise ValueError(f"Invalid polynomial coefficient: {line}") from exc

    if len(coeffs) != degree + 1:
        raise ValueError(
            f"Expected {degree + 1} polynomial coefficients, got {len(coeffs)}"
        )

    return degree, coeffs, i


def _parse_inline_polynomial(expr: str) -> List[int]:
    compact = re.sub(r"\s+", "", expr).replace("*", "").replace("x", "X")
    if not compact:
        raise ValueError("Polynomial expression is empty")

    term_re = re.compile(r"[+-]?\d+(?:X(?:\^\d+)?)?")
    matches = list(term_re.finditer(compact))
    if not matches:
        raise ValueError(f"Unable to parse polynomial: {expr}")

    consumed = "".join(m.group(0) for m in matches)
    if consumed != compact:
        raise ValueError(f"Invalid polynomial syntax: {expr}")

    coeff_map: Dict[int, int] = {}
    for m in matches:
        term = m.group(0)
        if "X" not in term:
            coeff = int(term)
            power = 0
        else:
            coeff_part, power_part = term.split("X", 1)
            coeff = int(coeff_part)
            power = 1
            if power_part:
                if not power_part.startswith("^"):
                    raise ValueError(f"Invalid power syntax in term: {term}")
                power = int(power_part[1:])
            if power < 0:
                raise ValueError("Negative exponents are not supported")
        coeff_map[power] = coeff_map.get(power, 0) + coeff

    max_power = max(coeff_map)
    coeffs = [0] * (max_power + 1)
    for power, coeff in coeff_map.items():
        coeffs[power] = coeff
    return coeffs


def eval_poly_mod(coeffs: List[int], m: int, modulus: int) -> int:
    acc = 0
    for coeff in reversed(coeffs):
        acc = (acc * m + coeff) % modulus
    return acc


def _content(coeffs: List[int]) -> int:
    g = 0
    for c in coeffs:
        g = math.gcd(g, abs(c))
    return g


def validate_config(path: str) -> List[str]:
    errors: List[str] = []

    try:
        params, polys = parse_sieve_cfg(path)
    except ValueError as exc:
        return [f"Parse error: {exc}"]

    for required in ("N", "m", "f1", "f2"):
        if required not in params:
            errors.append(f"Missing required setting: {required}")

    if errors:
        return errors

    try:
        n = int(params["N"])
    except ValueError:
        errors.append("N must be an integer")
        return errors

    try:
        m = int(params["m"])
    except ValueError:
        errors.append("m must be an integer")
        return errors

    if n <= 1:
        errors.append("N must be greater than 1")
    if n % 2 == 0:
        errors.append("N must be odd for GNFS")
    if m <= 0:
        errors.append("m must be positive")
    if m >= n:
        errors.append("m should be smaller than N")

    for name in ("f1", "f2"):
        if name not in polys:
            errors.append(f"Unable to parse {name} polynomial")
            continue

        poly = polys[name]
        if not poly.coeffs:
            errors.append(f"{name} has no coefficients")
            continue

        if poly.coeffs[-1] == 0:
            errors.append(f"{name} has zero leading coefficient")

        content = _content(poly.coeffs)
        if content != 1:
            errors.append(f"{name} is not primitive (content gcd = {content})")

        residue = eval_poly_mod(poly.coeffs, m % n if n else m, n) if n > 1 else 1
        if residue != 0:
            errors.append(f"{name}(m) mod N = {residue}, expected 0")

    if "f1" in polys and polys["f1"].degree < 2:
        errors.append("f1 must have degree >= 2 for GNFS")

    if "f2" in polys and polys["f2"].degree != 1:
        errors.append("f2 must be linear (degree 1) for GNFS lattice sieving")

    return errors


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", default="sieve.cfg", help="Path to sieve.cfg")
    args = parser.parse_args(argv)

    errors = validate_config(args.config)
    if errors:
        print(f"{args.config}: INVALID")
        for err in errors:
            print(f"  - {err}")
        return 1

    print(f"{args.config}: OK")
    print("  f1/f2 satisfy f(m) ≡ 0 (mod N) and pass GNFS suitability checks.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
