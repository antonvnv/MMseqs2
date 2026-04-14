#!/bin/bash
# Evaluate a single branch across all repro/data/*.fasta sequences, skipping done ones.
# Usage: ./eval_branch.sh <branch_name>
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BRANCH=$1

if [ -z "$BRANCH" ]; then
    echo "Usage: $0 <branch_name>"
    exit 1
fi

MAX_SEQS=${MAX_SEQS:-10000}
NUM_ITER=${NUM_ITER:-3}

echo "=== Branch: $BRANCH, MAX_SEQS=$MAX_SEQS, NUM_ITER=$NUM_ITER ==="

pushd "$REPO_ROOT" > /dev/null
git checkout "$BRANCH" 2>&1 | head -1
COMMIT=$(git rev-parse --short HEAD)
popd > /dev/null

echo "Commit: $COMMIT"

# Build
mkdir -p "$REPO_ROOT/build"
pushd "$REPO_ROOT/build" > /dev/null
if [ ! -f Makefile ]; then
    cmake -DHAVE_MPI=0 \
          -DENABLE_CUDA=1 \
          -DHAVE_TESTS=1 \
          -DHAVE_AVX2=1 \
          -DCMAKE_CUDA_ARCHITECTURES="75;80;86;89;90" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX=. \
          ..
fi
make -j64 2>&1 | tail -1
popd > /dev/null

RESULTS_BASE="$SCRIPT_DIR/results/${COMMIT}_maxseqs_${MAX_SEQS}_numiter_${NUM_ITER}"

for fasta in "$SCRIPT_DIR"/data/*.fasta; do
    name=$(basename "$fasta" .fasta)
    if [ -f "$RESULTS_BASE/$name/stats.txt" ]; then
        echo "  SKIP $name (already done)"
        continue
    fi
    echo "  RUN  $name"
    if MAX_SEQS=$MAX_SEQS NUM_ITER=$NUM_ITER QUERY_FASTA="$fasta" \
        "$SCRIPT_DIR/test_regression" 2>&1 | grep -aE '^(Query:|Hits:|Identity:|  >=)'; then
        echo ""
    else
        echo "  FAIL $name (exit code $?)"
    fi
done

echo "=== Done: $BRANCH ($COMMIT) ==="
