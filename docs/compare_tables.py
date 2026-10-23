#!/usr/bin/env python3
"""Compare CSV/TSV tables by unique keys (used for docs/validation.md).

Usage: compare_tables.py <original> <new> <key_col[,key_col2]> [old:new,...]
Renames apply to the original table; keys refer to the resulting column names.
Column/key sets must match, but column and row order may differ. Finite numbers
use abs(a-b) / max(abs(a), abs(b), 1e-300) < 1e-6. Both NaNs and same-sign
infinities compare equal; other nonfinite pairs differ. Empty strings and NA
compare equal, preserving the existing missing-value convention.

Exit 0: IDENTICAL or EQUAL within 1e-6; exit 1: DIFFERENT; exit 2: invalid input.
IDENTICAL means equality under these rules, not byte-for-byte file identity.
Empty/header-only tables, malformed rows, duplicate headers and keys are errors.
"""
import argparse
import csv
import math
import sys


RTOL = 1e-6


def parse_keys(value):
    keys = value.split(',')
    if any(not key.strip() for key in keys) or len(set(keys)) != len(keys):
        raise ValueError('key columns must be nonempty and unique')
    return keys


def parse_renames(value):
    renames = {}
    for entry in value.split(',') if value else []:
        parts = entry.split(':')
        if len(parts) != 2 or any(not part.strip() for part in parts):
            raise ValueError('column renames must have the form old:new,...')
        old, new = parts
        if old in renames:
            raise ValueError(f'duplicate rename source: {old!r}')
        renames[old] = new
    return renames


def validate_columns(columns, path):
    if not columns or any(not column.strip() for column in columns):
        raise ValueError(f'{path}: empty column name or missing header')
    if len(set(columns)) != len(columns):
        raise ValueError(f'{path}: duplicate column name (possibly after renaming)')


def load(path, rename, keys):
    delimiter = '\t' if path.endswith(('.tsv', '.txt')) else ','
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.reader(handle, delimiter=delimiter, strict=True)
        columns = next(reader, None)
        if columns is None:
            raise ValueError(f'{path}: empty input')
        validate_columns(columns, path)
        unknown = set(rename) - set(columns)
        if unknown:
            raise ValueError(f'{path}: rename source columns not found: {sorted(unknown)}')
        columns = [rename.get(column, column) for column in columns]
        validate_columns(columns, path)
        missing_keys = set(keys) - set(columns)
        if missing_keys:
            raise ValueError(f'{path}: missing key column(s): {sorted(missing_keys)}')

        rows = {}
        for values in reader:
            if len(values) != len(columns):
                raise ValueError(f'{path}:{reader.line_num}: expected {len(columns)} fields, got {len(values)}')
            row = dict(zip(columns, values))
            key = tuple(row[column] for column in keys)
            if any(not value.strip() for value in key):
                raise ValueError(f'{path}:{reader.line_num}: empty key value')
            if key in rows:
                raise ValueError(f'{path}:{reader.line_num}: duplicate key {key!r}')
            rows[key] = row
        if not rows:
            raise ValueError(f'{path}: no data rows')
    return rows, columns


def compare(original_path, new_path, keys, rename):
    original, original_columns = load(original_path, rename, keys)
    new, new_columns = load(new_path, {}, keys)
    original_keys, new_keys = set(original), set(new)
    print(f'{original_path}\n  vs {new_path}\n  rows {len(original)} / {len(new)}; '
          f'keys only in original {len(original_keys - new_keys)}, only in new {len(new_keys - original_keys)}')
    missing = [column for column in original_columns if column not in new_columns]
    extra = [column for column in new_columns if column not in original_columns]
    if missing:
        print('  columns of the original missing in new:', missing)
    if extra:
        print('  extra columns in new:', extra)

    worst = 0.0
    differences = 0
    for column in original_columns:
        if column not in new_columns:
            continue
        max_difference = 0.0
        text_differences = nonfinite_differences = 0
        for key in original_keys & new_keys:
            left, right = original[key][column], new[key][column]
            if left == right or (left in ('', 'NA') and right in ('', 'NA')):
                continue
            try:
                first, second = float(left), float(right)
            except ValueError:
                text_differences += 1
                continue
            if not math.isfinite(first) or not math.isfinite(second):
                if not ((math.isnan(first) and math.isnan(second)) or first == second):
                    nonfinite_differences += 1
                continue
            relative_difference = abs(first - second) / max(abs(first), abs(second), 1e-300)
            max_difference = max(max_difference, relative_difference)
        if max_difference > 0 or text_differences or nonfinite_differences:
            print(f'  {column}: max rel diff {max_difference:.3g}; '
                  f'text differences {text_differences}; nonfinite differences {nonfinite_differences}')
        worst = max(worst, max_difference)
        differences += text_differences + nonfinite_differences

    ok = not missing and not extra and original_keys == new_keys and differences == 0 and worst < RTOL
    result = 'DIFFERENT'
    if ok:
        result = 'IDENTICAL' if worst == 0 else f'EQUAL within 1e-6 (max rel diff {worst:.2g})'
    print('  RESULT:', result)
    return 0 if ok else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('original')
    parser.add_argument('new')
    parser.add_argument('key', help='key column name(s), separated by commas')
    parser.add_argument('renames', nargs='?', default='', help='original column renames: old:new,...')
    args = parser.parse_args()
    try:
        return compare(args.original, args.new, parse_keys(args.key), parse_renames(args.renames))
    except (OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
