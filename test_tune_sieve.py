#!/usr/bin/env python3
"""Unit tests for tune_sieve.py"""

import os
import random
import sys
import tempfile
import textwrap
import unittest

# Ensure the repo root is on the path so we can import tune_sieve
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tune_sieve


class TestParseSieveConfig(unittest.TestCase):
    """Tests for parse_sieve_config / write_sieve_config round-tripping."""

    SAMPLE_CONFIG = textwrap.dedent("""\
        # Configuration file for GNFS siever
        SIEVE_ID = tst250
        N = 3675041894739039405533259197211548846143110109152323761665377505538520830273
        m = 1039055185163739067158686044210305011518848843996071107667782459117743366330
        SKEWEDNESS = 1448288
        f1 = Polynomial
        DEGREE = 4
        1368110313470481075218020
        -248951329984989566
        -3210990608931
        6631
        1
        f2 = Polynomial
        DEGREE = 1
        -7786020182927250339
        541880117
        B1 = 1000000
        L1 = 50000000
        LP1 = 2
        B2 = 700000
        L2 = 50000000
        LP2 = 2
        SIEVE_BOUND_ADJUSTMENT1 = 11
        SIEVE_BOUND_ADJUSTMENT2 = 0
        SMALL_PRIME_BOUND1 = 500
        SMALL_PRIME_BOUND2 = 500
        MIN_A = -2000000000
        MAX_A = 2000000000
        MIN_B = 1
        MAX_B = 1000
        INITIAL_CUTOFF = 10
        RELATION_FILE = tst250.relations
    """)

    def _write_sample(self):
        fd, path = tempfile.mkstemp(suffix=".cfg")
        os.close(fd)
        with open(path, "w") as fh:
            fh.write(self.SAMPLE_CONFIG)
        return path

    def test_parse_reads_all_params(self):
        path = self._write_sample()
        try:
            _lines, params = tune_sieve.parse_sieve_config(path)
            self.assertEqual(params["SIEVE_ID"], "tst250")
            self.assertEqual(params["B1"], "1000000")
            self.assertEqual(params["SIEVE_BOUND_ADJUSTMENT1"], "11")
            self.assertEqual(params["INITIAL_CUTOFF"], "10")
            self.assertEqual(params["SMALL_PRIME_BOUND1"], "500")
        finally:
            os.unlink(path)

    def test_round_trip_preserves_structure(self):
        path = self._write_sample()
        try:
            lines, params = tune_sieve.parse_sieve_config(path)
            # Modify one parameter
            params["SIEVE_BOUND_ADJUSTMENT1"] = "15"
            out_path = path + ".out"
            tune_sieve.write_sieve_config(lines, params, out_path)

            with open(out_path) as fh:
                content = fh.read()

            # The modified parameter should appear
            self.assertIn("SIEVE_BOUND_ADJUSTMENT1 = 15", content)
            # Other parameters should be unchanged
            self.assertIn("B1 = 1000000", content)
            self.assertIn("INITIAL_CUTOFF = 10", content)
            # Polynomials should still be present
            self.assertIn("f1 = Polynomial", content)
            self.assertIn("DEGREE = 4", content)
            self.assertIn("1368110313470481075218020", content)
            # Comments should be preserved
            self.assertIn("# Configuration file", content)
        finally:
            os.unlink(path)
            if os.path.exists(path + ".out"):
                os.unlink(path + ".out")

    def test_write_updates_multiple_params(self):
        path = self._write_sample()
        try:
            lines, params = tune_sieve.parse_sieve_config(path)
            params["SIEVE_BOUND_ADJUSTMENT1"] = "20"
            params["SIEVE_BOUND_ADJUSTMENT2"] = "5"
            params["INITIAL_CUTOFF"] = "25"
            out_path = path + ".out"
            tune_sieve.write_sieve_config(lines, params, out_path)

            _lines2, params2 = tune_sieve.parse_sieve_config(out_path)
            self.assertEqual(params2["SIEVE_BOUND_ADJUSTMENT1"], "20")
            self.assertEqual(params2["SIEVE_BOUND_ADJUSTMENT2"], "5")
            self.assertEqual(params2["INITIAL_CUTOFF"], "25")
            # Unchanged
            self.assertEqual(params2["B1"], "1000000")
        finally:
            os.unlink(path)
            if os.path.exists(path + ".out"):
                os.unlink(path + ".out")


class TestParseRelationRate(unittest.TestCase):
    """Tests for _parse_relation_rate."""

    def test_parses_verbose_output(self):
        output = (
            "42 -> 10 -> 5 (123.4,444240,10661760), (456.7,1644120,39458880),\n"
            "38 -> 8 -> 3 (100.0,360000,8640000), (789.1,2840760,68178240),\n"
        )
        rate = tune_sieve._parse_relation_rate(output)
        self.assertAlmostEqual(rate, 789.1, places=1)

    def test_parses_single_line(self):
        output = "5 (963.5,3468600,83246400), (963.5,3468600,83246400),\n"
        rate = tune_sieve._parse_relation_rate(output)
        self.assertAlmostEqual(rate, 963.5, places=1)

    def test_returns_none_on_empty(self):
        self.assertIsNone(tune_sieve._parse_relation_rate(""))

    def test_returns_none_on_garbage(self):
        self.assertIsNone(tune_sieve._parse_relation_rate("no rate here\n"))

    def test_handles_scientific_notation(self):
        output = "5 (1.23e+02,444240,10661760), (4.56e+02,1644120,39458880),\n"
        rate = tune_sieve._parse_relation_rate(output)
        self.assertAlmostEqual(rate, 456.0, places=0)


class TestParseStatsCsv(unittest.TestCase):
    """Tests for _parse_stats_csv."""

    def test_parses_csv(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        os.close(fd)
        try:
            with open(path, "w") as fh:
                fh.write("elapsed_time,total_relations,current_rel_per_sec,"
                         "running_avg_rel_per_sec,rel_per_hour,rel_per_day\n")
                fh.write("10.5,500,50.0,47.6,171360,4112640\n")
                fh.write("21.0,1000,48.0,47.8,172080,4129920\n")
            rate = tune_sieve._parse_stats_csv(path)
            self.assertAlmostEqual(rate, 47.8, places=1)
        finally:
            os.unlink(path)

    def test_returns_none_for_missing_file(self):
        self.assertIsNone(tune_sieve._parse_stats_csv("/nonexistent/path.csv"))


class TestParamHelpers(unittest.TestCase):
    """Tests for _get_param / _set_param."""

    def test_get_param_default(self):
        self.assertEqual(tune_sieve._get_param({}, "X"), 0)

    def test_get_param(self):
        self.assertEqual(tune_sieve._get_param({"X": "42"}, "X"), 42)

    def test_set_param_clamps_low(self):
        params = {}
        tune_sieve._set_param(params, "SIEVE_BOUND_ADJUSTMENT1", -5)
        self.assertEqual(int(params["SIEVE_BOUND_ADJUSTMENT1"]), 0)

    def test_set_param_clamps_high(self):
        params = {}
        tune_sieve._set_param(params, "SIEVE_BOUND_ADJUSTMENT1", 999)
        self.assertEqual(int(params["SIEVE_BOUND_ADJUSTMENT1"]), 40)

    def test_set_param_normal(self):
        params = {}
        tune_sieve._set_param(params, "INITIAL_CUTOFF", 25)
        self.assertEqual(int(params["INITIAL_CUTOFF"]), 25)

    def test_set_param_min_a_clamps_high(self):
        params = {}
        tune_sieve._set_param(params, "MIN_A", 999)
        self.assertEqual(int(params["MIN_A"]), 0)

    def test_set_param_max_a_clamps_low(self):
        params = {}
        tune_sieve._set_param(params, "MAX_A", -999)
        self.assertEqual(int(params["MAX_A"]), 0)

    def test_set_param_min_b(self):
        params = {}
        tune_sieve._set_param(params, "MIN_B", 0)
        self.assertEqual(int(params["MIN_B"]), 1)

    def test_set_param_max_b(self):
        params = {}
        tune_sieve._set_param(params, "MAX_B", 500)
        self.assertEqual(int(params["MAX_B"]), 500)

    def test_set_param_min_a_normal(self):
        params = {}
        tune_sieve._set_param(params, "MIN_A", -2000000000)
        self.assertEqual(int(params["MIN_A"]), -2000000000)


class TestInlinePoly(unittest.TestCase):
    """Test that configs with inline (single-line) polynomials are handled."""

    INLINE_CONFIG = textwrap.dedent("""\
        SIEVE_ID = inline_test
        N = 12345
        m = 678
        f1 = 1 + 2 X + 3 X^2
        f2 = 4 + 5 X
        B1 = 500
        SIEVE_BOUND_ADJUSTMENT1 = 5
        INITIAL_CUTOFF = 10
    """)

    def test_round_trip_inline_poly(self):
        fd, path = tempfile.mkstemp(suffix=".cfg")
        os.close(fd)
        try:
            with open(path, "w") as fh:
                fh.write(self.INLINE_CONFIG)
            lines, params = tune_sieve.parse_sieve_config(path)
            params["SIEVE_BOUND_ADJUSTMENT1"] = "8"
            out_path = path + ".out"
            tune_sieve.write_sieve_config(lines, params, out_path)

            _lines2, params2 = tune_sieve.parse_sieve_config(out_path)
            self.assertEqual(params2["SIEVE_BOUND_ADJUSTMENT1"], "8")
            self.assertEqual(params2["INITIAL_CUTOFF"], "10")
        finally:
            os.unlink(path)
            if os.path.exists(path + ".out"):
                os.unlink(path + ".out")


class TestRandomPerturbation(unittest.TestCase):
    """Tests for random_perturbation."""

    def test_returns_valid_params(self):
        """Perturbed params should still be within allowed ranges."""
        params = {
            "SIEVE_BOUND_ADJUSTMENT1": "20",
            "SIEVE_BOUND_ADJUSTMENT2": "10",
            "SMALL_PRIME_BOUND1": "500",
            "SMALL_PRIME_BOUND2": "500",
            "INITIAL_CUTOFF": "50",
            "MIN_A": "-2000000000",
            "MAX_A": "2000000000",
            "MIN_B": "1",
            "MAX_B": "1000",
        }
        random.seed(42)
        result = tune_sieve.random_perturbation(params)
        for key, _step, lo, hi in tune_sieve.TUNABLE_PARAMS:
            val = int(result[key])
            self.assertGreaterEqual(val, lo, f"{key} below minimum")
            self.assertLessEqual(val, hi, f"{key} above maximum")

    def test_does_not_mutate_original(self):
        """Original params dict should be unchanged."""
        params = {"SIEVE_BOUND_ADJUSTMENT1": "20", "INITIAL_CUTOFF": "50"}
        original_copy = dict(params)
        tune_sieve.random_perturbation(params)
        self.assertEqual(params, original_copy)

    def test_clamps_at_boundaries(self):
        """Params at boundaries should stay within range after perturbation."""
        params = {
            "SIEVE_BOUND_ADJUSTMENT1": "0",
            "SIEVE_BOUND_ADJUSTMENT2": "40",
            "SMALL_PRIME_BOUND1": "0",
            "SMALL_PRIME_BOUND2": "5000",
            "INITIAL_CUTOFF": "0",
            "MIN_A": "0",
            "MAX_A": "0",
            "MIN_B": "1",
            "MAX_B": "1",
        }
        random.seed(123)
        for _ in range(20):
            result = tune_sieve.random_perturbation(params)
            for key, _step, lo, hi in tune_sieve.TUNABLE_PARAMS:
                val = int(result[key])
                self.assertGreaterEqual(val, lo, f"{key} below minimum")
                self.assertLessEqual(val, hi, f"{key} above maximum")


class TestHillClimbRandomRestarts(unittest.TestCase):
    """Tests for the random restart behaviour of hill_climb."""

    def test_random_restarts_zero_stops_immediately(self):
        """With random_restarts=0, hill_climb should stop at first local max."""
        call_count = [0]

        def fake_evaluate(params, lines, config_path, lsieve_path,
                          min_q, max_q, timeout=600, **kwargs):
            call_count[0] += 1
            # Baseline call returns a rate; all others return worse.
            if call_count[0] == 1:
                return 100.0
            return 50.0

        import unittest.mock
        with unittest.mock.patch("tune_sieve.evaluate", side_effect=fake_evaluate):
            params = {
                "SIEVE_BOUND_ADJUSTMENT1": "20",
                "SIEVE_BOUND_ADJUSTMENT2": "10",
                "SMALL_PRIME_BOUND1": "500",
                "SMALL_PRIME_BOUND2": "500",
                "INITIAL_CUTOFF": "50",
                "MIN_A": "-2000000000",
                "MAX_A": "2000000000",
                "MIN_B": "1",
                "MAX_B": "1000",
            }
            calls_before = call_count[0]
            _best_params, best_rate = tune_sieve.hill_climb(
                [], params, "/fake/cfg", "/fake/lsieve",
                1000000, 1000500, max_iterations=3,
                random_restarts=0,
            )
        self.assertAlmostEqual(best_rate, 100.0)
        # With random_restarts=0 there should be no random-restart evaluate
        # calls after the greedy phase of iteration 1 fails.  The total
        # calls are: 1 (baseline) + N (greedy tries) + 0 (restarts).
        # We just verify the count doesn't include extra restart calls.
        total = call_count[0] - calls_before
        # Greedy tries = up to 2 per param (some may be skipped at boundary)
        max_greedy = 2 * len(tune_sieve.TUNABLE_PARAMS)
        self.assertLessEqual(total, 1 + max_greedy)

    def test_random_restart_finds_improvement(self):
        """If a random perturbation beats best_rate, hill-climbing continues."""
        call_count = [0]
        greedy_per_iter = [0]  # track how many greedy calls iteration 1 made

        def fake_evaluate(params, lines, config_path, lsieve_path,
                          min_q, max_q, timeout=600, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                return 100.0  # baseline
            return 50.0  # everything else fails

        import unittest.mock

        # First run: figure out how many greedy calls one iteration makes.
        with unittest.mock.patch("tune_sieve.evaluate", side_effect=fake_evaluate):
            params = {
                "SIEVE_BOUND_ADJUSTMENT1": "20",
                "SIEVE_BOUND_ADJUSTMENT2": "10",
                "SMALL_PRIME_BOUND1": "500",
                "SMALL_PRIME_BOUND2": "500",
                "INITIAL_CUTOFF": "50",
                "MIN_A": "-2000000000",
                "MAX_A": "2000000000",
                "MIN_B": "1",
                "MAX_B": "1000",
            }
            tune_sieve.hill_climb(
                [], dict(params), "/fake/cfg", "/fake/lsieve",
                1000000, 1000500, max_iterations=1,
                random_restarts=0,
            )
        # calls: 1 baseline + N greedy
        greedy_count = call_count[0] - 1

        # Second run: make the first random restart after greedy-fail succeed.
        call_count[0] = 0

        def fake_evaluate2(params, lines, config_path, lsieve_path,
                           min_q, max_q, timeout=600, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                return 100.0  # baseline
            # First random restart attempt of iteration 1:
            if call_count[0] == 1 + greedy_count + 1:
                return 200.0
            return 50.0

        with unittest.mock.patch("tune_sieve.evaluate", side_effect=fake_evaluate2):
            _best_params, best_rate = tune_sieve.hill_climb(
                [], dict(params), "/fake/cfg", "/fake/lsieve",
                1000000, 1000500, max_iterations=2,
                random_restarts=3,
            )
        self.assertAlmostEqual(best_rate, 200.0)

    def test_random_restart_verbose_shows_param_changes(self):
        """Verbose output for random restarts should show changed parameters."""
        call_count = [0]

        def fake_evaluate(params, lines, config_path, lsieve_path,
                          min_q, max_q, timeout=600, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                return 100.0  # baseline
            return 50.0

        import io
        import unittest.mock

        params = {
            "SIEVE_BOUND_ADJUSTMENT1": "20",
            "SIEVE_BOUND_ADJUSTMENT2": "10",
            "SMALL_PRIME_BOUND1": "500",
            "SMALL_PRIME_BOUND2": "500",
            "INITIAL_CUTOFF": "50",
            "MIN_A": "-2000000000",
            "MAX_A": "2000000000",
            "MIN_B": "1",
            "MAX_B": "1000",
        }
        captured = io.StringIO()
        random.seed(42)
        with unittest.mock.patch("tune_sieve.evaluate", side_effect=fake_evaluate), \
             unittest.mock.patch("sys.stdout", captured):
            tune_sieve.hill_climb(
                [], dict(params), "/fake/cfg", "/fake/lsieve",
                1000000, 1000500, max_iterations=1,
                verbose=True, random_restarts=3,
            )
        output = captured.getvalue()
        # The random restart lines should contain "[" with parameter changes
        restart_lines = [l for l in output.splitlines()
                         if "random restart attempt" in l]
        self.assertTrue(len(restart_lines) > 0,
                        "Expected random restart lines in verbose output")
        for line in restart_lines:
            self.assertIn("[", line,
                          "Random restart line should show parameter changes")
            self.assertIn("]", line)


if __name__ == "__main__":
    unittest.main()
