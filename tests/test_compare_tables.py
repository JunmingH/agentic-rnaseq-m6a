"""Synthetic CLI regression tests; no research data or third-party packages.

Set COMPARE_TABLES_SCRIPT to run the same checks against another script version.
"""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(os.environ.get(
    "COMPARE_TABLES_SCRIPT",
    str(Path(__file__).resolve().parents[1] / "docs" / "compare_tables.py"),
)).resolve()


class CompareTablesTests(unittest.TestCase):
    def compare(self, original, new, key="id", renames=None, suffix=".csv"):
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / ("original" + suffix)
            second = Path(directory) / ("new" + suffix)
            first.write_text(original, encoding="utf-8")
            second.write_text(new, encoding="utf-8")
            command = [sys.executable, str(SCRIPT), str(first), str(second), key]
            if renames is not None:
                command.append(renames)
            return subprocess.run(command, capture_output=True, text=True, check=False)

    def assert_result(self, result, code, marker):
        self.assertEqual(result.returncode, code, result.stdout + result.stderr)
        self.assertIn(marker, result.stdout + result.stderr)
        if code:
            self.assertNotIn("RESULT: IDENTICAL", result.stdout)
            self.assertNotIn("RESULT: EQUAL", result.stdout)

    def test_missing_original_column_fails(self):
        result = self.compare("id,value\na,1\n", "id\na\n")
        self.assert_result(result, 1, "DIFFERENT")
        self.assertIn("missing in new", result.stdout)

    def test_one_sided_nan_fails(self):
        result = self.compare("id,value\na,1\n", "id,value\na,NaN\n")
        self.assert_result(result, 1, "DIFFERENT")

    def test_equal_tables_with_reordered_rows_and_columns(self):
        result = self.compare("id,value\na,1\nb,2\n", "value,id\n2,b\n1,a\n")
        self.assert_result(result, 0, "RESULT: IDENTICAL")

    def test_equal_numeric_spellings_and_missing_values(self):
        result = self.compare("id,value\na,1.0\nb,\n", "id,value\na,1\nb,NA\n")
        self.assert_result(result, 0, "RESULT: IDENTICAL")

    def test_extra_column_fails(self):
        result = self.compare("id,value\na,1\n", "id,value,extra\na,1,x\n")
        self.assert_result(result, 1, "DIFFERENT")

    def test_missing_and_extra_rows_fail(self):
        for original, new in (
            ("id,value\na,1\nb,2\n", "id,value\na,1\n"),
            ("id,value\na,1\n", "id,value\na,1\nb,2\n"),
            ("id,value\na,1\n", "id,value\nb,1\n"),
        ):
            with self.subTest(original=original, new=new):
                self.assert_result(self.compare(original, new), 1, "DIFFERENT")

    def test_normal_numeric_and_text_differences_fail(self):
        for value in ("2", "different", "", "NA"):
            with self.subTest(value=value):
                self.assert_result(self.compare("id,value\na,1\n", "id,value\na," + value + "\n"), 1, "DIFFERENT")

    def test_within_relative_tolerance_passes(self):
        result = self.compare("id,value\na,1\n", "id,value\na,1.0000005\n")
        self.assert_result(result, 0, "RESULT: EQUAL within 1e-6")

    def test_outside_relative_tolerance_fails(self):
        result = self.compare("id,value\na,1\n", "id,value\na,1.000002\n")
        self.assert_result(result, 1, "DIFFERENT")

    def test_tolerance_is_symmetric(self):
        for left, right in (("1.0000005", "1"), ("1.000002", "1")):
            with self.subTest(left=left):
                result = self.compare("id,value\na," + left + "\n", "id,value\na," + right + "\n")
                self.assert_result(result, 0 if left == "1.0000005" else 1,
                                   "RESULT: EQUAL" if left == "1.0000005" else "DIFFERENT")

    def test_composite_keys(self):
        original = "id,site,value\na,1,10\na,2,20\n"
        new = "id,site,value\na,2,20\na,1,10\n"
        self.assert_result(self.compare(original, new, key="id,site"), 0, "IDENTICAL")

    def test_original_column_renames_including_key(self):
        result = self.compare("old_id,old_value\na,1\n", "id,value\na,1\n",
                              renames="old_id:id,old_value:value")
        self.assert_result(result, 0, "IDENTICAL")

    def test_empty_rename_argument_is_compatible(self):
        self.assert_result(self.compare("id,value\na,1\n", "id,value\na,1\n", renames=""), 0, "IDENTICAL")

    def test_tsv_and_txt_input(self):
        for suffix in (".tsv", ".txt"):
            with self.subTest(suffix=suffix):
                self.assert_result(self.compare("id\tvalue\na\t1\n", "id\tvalue\na\t1\n", suffix=suffix), 0, "IDENTICAL")

    def test_quoted_csv_field(self):
        data = 'id,value\n"a,b","text, value"\n'
        self.assert_result(self.compare(data, data), 0, "IDENTICAL")

    def test_nonfinite_comparisons(self):
        cases = (
            ("NaN", "nan", 0), ("NaN", "1", 1), ("1", "NaN", 1),
            ("inf", "Infinity", 0), ("-inf", "-Infinity", 0),
            ("inf", "-inf", 1), ("inf", "1", 1), ("1", "inf", 1),
            ("NaN", "inf", 1), ("inf", "NaN", 1),
            ("1e308", "-1e308", 1), ("0", "-0", 0),
        )
        for left, right, code in cases:
            with self.subTest(left=left, right=right):
                result = self.compare("id,value\na," + left + "\n", "id,value\na," + right + "\n")
                self.assert_result(result, code, "IDENTICAL" if code == 0 else "DIFFERENT")

    def test_duplicate_keys_in_either_table_are_errors(self):
        duplicate = "id,value\na,1\na,1\n"
        valid = "id,value\na,1\n"
        for original, new in ((duplicate, valid), (valid, duplicate)):
            with self.subTest(original=original, new=new):
                self.assert_result(self.compare(original, new), 2, "duplicate key")

    def test_duplicate_composite_key_is_error(self):
        data = "id,site,value\na,1,10\na,1,10\n"
        self.assert_result(self.compare(data, data, key="id,site"), 2, "duplicate key")

    def test_duplicate_headers_in_either_table_are_errors(self):
        duplicate = "id,value,value\na,1,2\n"
        valid = "id,value\na,2\n"
        for original, new in ((duplicate, valid), (valid, duplicate)):
            with self.subTest(original=original, new=new):
                self.assert_result(self.compare(original, new), 2, "duplicate column")

    def test_rename_collision_is_error(self):
        result = self.compare("id,value\na,1\n", "id\na\n", renames="value:id")
        self.assert_result(result, 2, "duplicate column")

    def test_invalid_renames_are_errors(self):
        for renames in ("value", "value:", ":value", "value:v,value:w", "missing:value"):
            with self.subTest(renames=renames):
                result = self.compare("id,value\na,1\n", "id,value\na,1\n", renames=renames)
                self.assert_result(result, 2, "ERROR:")

    def test_empty_or_malformed_input_is_error(self):
        bad_tables = (
            "", "\n", "id,value\n", "id,\na,1\n",
            "id,value\na\n", "id,value\na,1,extra\n",
            'id,value\na,"unterminated\n', "id,value\n\na,1\n",
        )
        valid = "id,value\na,1\n"
        for malformed in bad_tables:
            for original, new in ((malformed, valid), (valid, malformed)):
                with self.subTest(original=original, new=new):
                    self.assert_result(self.compare(original, new), 2, "ERROR:")

    def test_missing_key_column_is_error(self):
        result = self.compare("id,value\na,1\n", "other,value\na,1\n")
        self.assert_result(result, 2, "missing key column")

    def test_empty_key_value_is_error(self):
        for key_value in ("", " "):
            with self.subTest(key_value=key_value):
                result = self.compare("id,value\na,1\n", "id,value\n" + key_value + ",1\n")
                self.assert_result(result, 2, "empty key")

    def test_invalid_key_argument_is_error(self):
        for key in ("", "id,", "id,id"):
            with self.subTest(key=key):
                self.assert_result(self.compare("id,value\na,1\n", "id,value\na,1\n", key=key), 2, "ERROR:")

    def test_missing_file_is_error(self):
        with tempfile.TemporaryDirectory() as directory:
            nonexistent = str(Path(directory) / "does-not-exist.csv")
            result = subprocess.run([sys.executable, str(SCRIPT), nonexistent, nonexistent, "id"],
                                    capture_output=True, text=True, check=False)
        self.assert_result(result, 2, "ERROR:")

    def test_missing_arguments_report_usage(self):
        result = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("usage:", result.stderr.lower())
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
