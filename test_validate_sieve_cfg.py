#!/usr/bin/env python3

import os
import tempfile
import textwrap
import unittest

import validate_sieve_cfg


class ValidateSieveCfgTests(unittest.TestCase):
    VALID_MULTILINE = textwrap.dedent(
        """\
        SIEVE_ID = tst250
        N = 3675041894739039405533259197211548846143110109152323761665377505538520830273
        m = 1039055185163739067158686044210305011518848843996071107667782459117743366330
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
        """
    )

    VALID_INLINE = textwrap.dedent(
        """\
        N = 308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341
        m = 270408018939143932409845114498368658401589778718225878946223461855504808192356812169962599904520634544716909920361370733421226675947092513281426
        f1 = -145239798163761036421276279674832 - 1214966661472423650426575068 X + 345797467977205840635500 X^2 + 576326033846910003 X^3 - 77861900597420 X^4 + 61461120 X^5
        f2 = -1380597316096461543125283587 + 17878358832686441 X
        """
    )

    def _write_cfg(self, content: str) -> str:
        fd, path = tempfile.mkstemp(suffix=".cfg")
        os.close(fd)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(content)
        return path

    def test_valid_multiline_polynomials(self):
        path = self._write_cfg(self.VALID_MULTILINE)
        try:
            errors = validate_sieve_cfg.validate_config(path)
            self.assertEqual(errors, [])
        finally:
            os.unlink(path)

    def test_valid_inline_polynomials(self):
        path = self._write_cfg(self.VALID_INLINE)
        try:
            errors = validate_sieve_cfg.validate_config(path)
            self.assertEqual(errors, [])
        finally:
            os.unlink(path)

    def test_detects_bad_mod_condition(self):
        bad = self.VALID_MULTILINE.replace("541880117", "541880118")
        path = self._write_cfg(bad)
        try:
            errors = validate_sieve_cfg.validate_config(path)
            self.assertTrue(any("f2(m) mod N" in e for e in errors))
        finally:
            os.unlink(path)

    def test_detects_unsuitable_degrees(self):
        bad = textwrap.dedent(
            """\
            N = 101
            m = 10
            f1 = 1 + 1 X
            f2 = 1 + 2 X + 3 X^2
            """
        )
        path = self._write_cfg(bad)
        try:
            errors = validate_sieve_cfg.validate_config(path)
            self.assertIn("f1 must have degree >= 2 for GNFS", errors)
            self.assertIn("f2 must be linear (degree 1) for GNFS lattice sieving", errors)
        finally:
            os.unlink(path)

    def test_detects_non_primitive(self):
        bad = textwrap.dedent(
            """\
            N = 101
            m = 2
            f1 = 2 + 4 X + 6 X^2
            f2 = 4 + 6 X
            """
        )
        path = self._write_cfg(bad)
        try:
            errors = validate_sieve_cfg.validate_config(path)
            self.assertTrue(any("not primitive" in e for e in errors))
        finally:
            os.unlink(path)

    def test_inline_parser_supports_implicit_x_coefficients(self):
        coeffs = validate_sieve_cfg._parse_inline_polynomial("1 - X + X^2")
        self.assertEqual(coeffs, [1, -1, 1])


if __name__ == "__main__":
    unittest.main()
