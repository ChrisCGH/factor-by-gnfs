#!/usr/bin/env python3
"""Unit tests for trial_sieve.py"""

import os
import sys
import tempfile
import textwrap
import unittest

# Ensure the gnfs directory is on the path so we can import trial_sieve
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "gnfs"))
import trial_sieve


class TestParseSkewedOutput(unittest.TestCase):
    """Tests for parse_skewed_output."""

    SAMPLE_OUTPUT = textwrap.dedent("""\
        f1 = 28401178922957826568481403296310 - 3199761622362244507168766869 X - 68308394838556620909345 X^2 + 2832627855183353597 X^3 + 25714112986408 X^4 + 158988960 X^5
        m = 265326213255482044308304372573181143657052889388773311344205525057953532130583924590773156033550113811312873866059392293321271779657292494557824
        a = 62626764709090441
        b = 1141593610769955636830559061
        s = 35743
        alpha = -4.50856
        E(F) = 40.9678
        N = 308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341

        ==============================================================================

        f1 = 12345678901234567890 + 9876543210 X + 111213 X^2 + 42 X^3
        m = 99999999999999999999
        a = 12345
        b = 67890
        s = 10000
        alpha = -3.21
        E(F) = 38.5
        N = 308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341

        ==============================================================================
    """)

    def _write_sample(self, content=None):
        fd, path = tempfile.mkstemp(suffix=".out")
        os.close(fd)
        with open(path, "w") as fh:
            fh.write(self.SAMPLE_OUTPUT if content is None else content)
        return path

    def test_parses_two_blocks(self):
        path = self._write_sample()
        try:
            blocks = trial_sieve.parse_skewed_output(path)
            self.assertEqual(len(blocks), 2)
        finally:
            os.unlink(path)

    def test_first_block_values(self):
        path = self._write_sample()
        try:
            blocks = trial_sieve.parse_skewed_output(path)
            b = blocks[0]
            self.assertIn("f1", b)
            self.assertIn("158988960 X^5", b["f1"])
            self.assertEqual(b["a"], "62626764709090441")
            self.assertEqual(b["b"], "1141593610769955636830559061")
            self.assertEqual(b["s"], "35743")
            self.assertEqual(b["alpha"], "-4.50856")
            self.assertEqual(b["E_F"], "40.9678")
            self.assertTrue(b["m"].startswith("265326213"))
            self.assertTrue(b["N"].startswith("308268017"))
        finally:
            os.unlink(path)

    def test_second_block_values(self):
        path = self._write_sample()
        try:
            blocks = trial_sieve.parse_skewed_output(path)
            b = blocks[1]
            self.assertEqual(b["a"], "12345")
            self.assertEqual(b["b"], "67890")
            self.assertEqual(b["s"], "10000")
            self.assertEqual(b["alpha"], "-3.21")
            self.assertEqual(b["E_F"], "38.5")
        finally:
            os.unlink(path)

    def test_parses_block_without_trailing_separator(self):
        content = textwrap.dedent("""\
            f1 = 100 + 200 X
            m = 999
            a = 10
            b = 20
            s = 100
            alpha = -1.0
            E(F) = 30.0
            N = 12345
        """)
        path = self._write_sample(content)
        try:
            blocks = trial_sieve.parse_skewed_output(path)
            self.assertEqual(len(blocks), 1)
            self.assertEqual(blocks[0]["a"], "10")
        finally:
            os.unlink(path)

    def test_parses_empty_file(self):
        path = self._write_sample("")
        try:
            blocks = trial_sieve.parse_skewed_output(path)
            self.assertEqual(len(blocks), 0)
        finally:
            os.unlink(path)

    def test_handles_extra_lines(self):
        """Extra lines (like f1 evaluation result) should not corrupt parsing."""
        content = textwrap.dedent("""\
            f1 = 100 + 200 X
            m = 999
            a = 10
            b = 20
            s = 100
            alpha = -1.0
            E(F) = 30.0
            N = 12345
            0

            ==============================================================================
        """)
        path = self._write_sample(content)
        try:
            blocks = trial_sieve.parse_skewed_output(path)
            self.assertEqual(len(blocks), 1)
            self.assertEqual(blocks[0]["f1"], "100 + 200 X")
            self.assertEqual(blocks[0]["N"], "12345")
        finally:
            os.unlink(path)


class TestGenerateConfig(unittest.TestCase):
    """Tests for generate_config."""

    TEMPLATE = textwrap.dedent("""\
        # Configuration file for GNFS siever
        SIEVE_ID = C144_104_45
        N = <N>
        m = <m>
        # alpha = <alpha>
        # E(F) = <E_F>
        SKEWEDNESS = <s>
        f1 = <f1>
        f2 = -<b> + <a> X
        B1 = 16000000
        L1 = 300000000
        LP1 = 2
        B2 = 16000000
        L2 = 300000000
        LP2 = 2
        SIEVE_BOUND_ADJUSTMENT1 = 26
        SIEVE_BOUND_ADJUSTMENT2 = 6
        SMALL_PRIME_BOUND1 = 1000
        SMALL_PRIME_BOUND2 = 1000
        MIN_A = -2600000000
        MAX_A = 2600000000
        MIN_B = 1
        MAX_B = 80000
        INITIAL_CUTOFF = 15
        RELATION_FILE = C144_104_45.relations
    """)

    def test_substitution(self):
        block = {
            "f1": "100 + 200 X + 300 X^2",
            "m": "999888777",
            "a": "12345",
            "b": "67890",
            "s": "5000",
            "alpha": "-2.5",
            "E_F": "35.0",
            "N": "123456789",
        }
        result = trial_sieve.generate_config(self.TEMPLATE, block)
        self.assertIn("N = 123456789", result)
        self.assertIn("m = 999888777", result)
        self.assertIn("f1 = 100 + 200 X + 300 X^2", result)
        self.assertIn("f2 = -67890 + 12345 X", result)
        self.assertIn("SKEWEDNESS = 5000", result)
        self.assertIn("# alpha = -2.5", result)
        self.assertIn("# E(F) = 35.0", result)
        # Fixed values should remain
        self.assertIn("B1 = 16000000", result)
        self.assertIn("INITIAL_CUTOFF = 15", result)

    def test_missing_fields_become_empty(self):
        block = {"f1": "100 + 200 X", "m": "999"}
        result = trial_sieve.generate_config(self.TEMPLATE, block)
        self.assertIn("N = \n", result)
        self.assertIn("f1 = 100 + 200 X", result)

    def test_preserves_non_placeholder_content(self):
        block = {
            "f1": "1", "m": "2", "a": "3", "b": "4",
            "s": "5", "alpha": "6", "E_F": "7", "N": "8",
        }
        result = trial_sieve.generate_config(self.TEMPLATE, block)
        self.assertIn("RELATION_FILE = C144_104_45.relations", result)
        self.assertIn("# Configuration file for GNFS siever", result)


class TestReadTemplate(unittest.TestCase):
    """Tests for read_template."""

    def test_reads_file_content(self):
        fd, path = tempfile.mkstemp(suffix=".template")
        os.close(fd)
        try:
            with open(path, "w") as fh:
                fh.write("N = <N>\nf1 = <f1>\n")
            content = trial_sieve.read_template(path)
            self.assertEqual(content, "N = <N>\nf1 = <f1>\n")
        finally:
            os.unlink(path)


class TestCountRelations(unittest.TestCase):
    """Tests for count_relations."""

    def test_counts_from_relation_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rel_path = os.path.join(tmpdir, "test.relations")
            with open(rel_path, "w") as fh:
                fh.write("# header\n")
                fh.write("relation1\n")
                fh.write("relation2\n")
                fh.write("relation3\n")
            count = trial_sieve.count_relations(tmpdir, "", "")
            self.assertEqual(count, 3)

    def test_counts_from_stdout(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            stdout = "42 -> 10 -> 5\n38 -> 8 -> 3\n"
            count = trial_sieve.count_relations(tmpdir, stdout, "")
            self.assertEqual(count, 8)

    def test_returns_zero_with_no_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            count = trial_sieve.count_relations(tmpdir, "", "")
            self.assertEqual(count, 0)

    def test_ignores_comment_lines_in_relations(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rel_path = os.path.join(tmpdir, "test.relations")
            with open(rel_path, "w") as fh:
                fh.write("# comment\n")
                fh.write("\n")
                fh.write("relation1\n")
            count = trial_sieve.count_relations(tmpdir, "", "")
            self.assertEqual(count, 1)


class TestFormatBlockSummary(unittest.TestCase):
    """Tests for format_block_summary."""

    def test_format_with_all_fields(self):
        block = {"E_F": "40.0", "alpha": "-3.5", "s": "5000"}
        result = trial_sieve.format_block_summary(1, block, 42)
        self.assertIn("Poly   1", result)
        self.assertIn("relations=    42", result)
        self.assertIn("E(F)=", result)
        self.assertIn("alpha=", result)

    def test_format_with_none_relations(self):
        block = {"E_F": "40.0", "alpha": "-3.5", "s": "5000"}
        result = trial_sieve.format_block_summary(1, block, None)
        self.assertIn("ERROR", result)


class TestEndToEndDryRun(unittest.TestCase):
    """End-to-end test with --dry-run."""

    def test_dry_run_generates_configs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create skewed output
            skewed_path = os.path.join(tmpdir, "skewed.out")
            with open(skewed_path, "w") as fh:
                fh.write(textwrap.dedent("""\
                    f1 = 100 + 200 X + 300 X^2
                    m = 999
                    a = 10
                    b = 20
                    s = 100
                    alpha = -1.0
                    E(F) = 30.0
                    N = 12345

                    ==============================================================================
                """))

            # Create template
            template_path = os.path.join(tmpdir, "sieve.cfg.template")
            with open(template_path, "w") as fh:
                fh.write(textwrap.dedent("""\
                    N = <N>
                    m = <m>
                    SKEWEDNESS = <s>
                    f1 = <f1>
                    f2 = -<b> + <a> X
                """))

            workdir = os.path.join(tmpdir, "work")

            # Run in dry-run mode using the function directly
            blocks = trial_sieve.parse_skewed_output(skewed_path)
            template = trial_sieve.read_template(template_path)

            os.makedirs(workdir, exist_ok=True)
            for i, block in enumerate(blocks, 1):
                config = trial_sieve.generate_config(template, block)
                config_path = os.path.join(workdir, f"poly_{i}", "sieve.cfg")
                os.makedirs(os.path.dirname(config_path), exist_ok=True)
                with open(config_path, "w") as fh:
                    fh.write(config)

            # Verify the generated config
            config_path = os.path.join(workdir, "poly_1", "sieve.cfg")
            self.assertTrue(os.path.isfile(config_path))
            with open(config_path, "r") as fh:
                content = fh.read()
            self.assertIn("N = 12345", content)
            self.assertIn("m = 999", content)
            self.assertIn("SKEWEDNESS = 100", content)
            self.assertIn("f1 = 100 + 200 X + 300 X^2", content)
            self.assertIn("f2 = -20 + 10 X", content)


class TestMultipleBlocks(unittest.TestCase):
    """Test parsing multiple polynomial blocks."""

    def test_three_blocks(self):
        content = ""
        for i in range(1, 4):
            content += textwrap.dedent(f"""\
                f1 = {i*100} + {i*10} X
                m = {i*1000}
                a = {i}
                b = {i*10}
                s = {i*100}
                alpha = -{i}.0
                E(F) = {30+i}.0
                N = {99999-i}

                ==============================================================================

            """)

        fd, path = tempfile.mkstemp(suffix=".out")
        os.close(fd)
        try:
            with open(path, "w") as fh:
                fh.write(content)
            blocks = trial_sieve.parse_skewed_output(path)
            self.assertEqual(len(blocks), 3)
            self.assertEqual(blocks[0]["a"], "1")
            self.assertEqual(blocks[1]["a"], "2")
            self.assertEqual(blocks[2]["a"], "3")
            self.assertEqual(blocks[0]["E_F"], "31.0")
            self.assertEqual(blocks[1]["E_F"], "32.0")
            self.assertEqual(blocks[2]["E_F"], "33.0")
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
