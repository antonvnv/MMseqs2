# E-value Snowball Regression Reproduction

This directory contains scripts and data to reproduce the e-value snowball
regression introduced in the `new` branch.

## The Bug

The e-value formula in `src/alignment/EvalueComputation.h` was changed from:

    0.038 * pow(epa * a, 0.83)   // old: sublinear dampening

to:

    epa * a / 170.0              // new: linear, too permissive

This makes weak hits pass the e-value threshold (0.1), which contaminates
the PSSM profile in iterative search, causing a snowball effect:

- **old branch**: 9 hits for URS0003421724 (all genuine, median 99% identity)
- **new branch**: ~2700 hits (median ~10% identity, mostly noise)

## Quick Reproduction

```bash
# From this directory:
./test_regression
# Expected: ~2700 hits with median ~10% identity (the regression)

# To compare with the old (correct) behavior:
# 1. Check out the old branch: cd .. && git checkout old
# 2. Re-run: cd repro && ./test_regression
# Expected: 9 hits with median ~99% identity
```

## Prerequisites

- Pre-indexed `rnacentral_active` database (set `TARGET_DB` env var if not
  at the default location `../../../rnacentral_active/db`)
- CUDA-capable GPU (set `GPU=0` for CPU, but it will be very slow)
- CMake, make, CUDA toolkit

## Files

- `test_regression` - Main test script (builds mmseqs, runs search, reports stats)
- `seq_identity.py` - Calculates sequence identity statistics from a3m output
- `eval_branch.sh` - Evaluates a branch across all 25 queries in data/
- `URS0003421724.fasta` - Default query (590nt rRNA, the original regression case)
- `data/` - 25 diverse RNA query sequences for broader evaluation
- `regression.txt` - Detailed root cause analysis
- `regression_fix_final_summary.txt` - Evaluation of fix approaches

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `TARGET_DB` | `../../../rnacentral_active/db` | Path to pre-indexed target database |
| `QUERY_FASTA` | `URS0003421724.fasta` | Query sequence |
| `MAX_SEQS` | `10000` | Maximum sequences per search iteration |
| `NUM_ITER` | `3` | Number of profile search iterations |
| `GPU` | `1` | Use GPU (1) or CPU (0) |
| `THREADS` | `64` | Number of threads |
