# Contributing to Cerberus

## Development setup

```bash
git clone https://github.com/iowa69/cerberus
cd cerberus
conda env create -f environment.yml
conda activate cerberus
pip install -e ".[dev]"
```

## Before opening a pull request

```bash
ruff check .          # lint
pytest -q             # unit + regression tests (no external tools needed)
bash scripts/smoke_test.sh   # end-to-end against miniature references
```

CI runs all three on Python 3.10-3.13.

## What a change needs

**A test that fails without it.** `tests/test_regressions.py` is organised one
test per historical defect, each with a docstring naming the wrong behaviour it
prevents. Add to it in the same style.

**Honest documentation.** If a change alters what the tool guarantees — read
retention, host removal, output naming — update `README.md` and `CHANGELOG.md`
in the same commit. Cerberus is used to decide whether human sequence data is
safe to publish, so overstating a guarantee is a bug of the same severity as a
crash.

**Tuning parameters in the table.** Thresholds live in `_LENGTH_BUCKETS` and
`_BASE_PARAMS` in `cerberus/autotune.py` and nowhere else, so the behaviour
stays auditable.

## Areas that need care

- `cerberus/stages/align.py` — SAM flag arithmetic is pair-level logic
  expressed through per-record flags. Read the module docstring first, and
  verify any change against real BAMs where exactly one mate maps.
- `cerberus/utils/shell.py` — every stage of a pipe must have its exit status
  checked. A silently truncated dataset looks like a biological result.
- `cerberus/pipelines/gdpr.py` — changes here affect whether human sequence
  reaches a public release.
