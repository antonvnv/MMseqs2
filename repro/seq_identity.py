#!/usr/bin/env python3

"""Calculate sequence identity statistics for an a3m file.

Usage: python3 seq_identity.py result.a3m
"""

import sys
from dataclasses import dataclass


@dataclass
class Sequence:
    hdr: str
    seq: str


def parse_a3m(alignment_a3m: str) -> list[Sequence]:
    if not alignment_a3m or not alignment_a3m.strip():
        return []

    sequences = []
    lines = alignment_a3m.strip().split('\n')
    current_hdr = None
    current_seq = []

    for line in lines:
        if line.startswith('>'):
            if current_hdr is not None:
                sequences.append(Sequence(hdr=current_hdr, seq=''.join(current_seq)))
            current_hdr = line[1:]
            current_seq = []
        else:
            current_seq.append(line)

    if current_hdr is not None:
        sequences.append(Sequence(hdr=current_hdr, seq=''.join(current_seq)))

    return sequences


def remove_a3m_insertions(seq: str) -> str:
    return ''.join(c for c in seq if not c.islower())


def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    seq1 = remove_a3m_insertions(seq1)
    seq2 = remove_a3m_insertions(seq2)

    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1, seq2 = seq1[:min_len], seq2[:min_len]

    if not seq1:
        return 0.0

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')
    non_gap_positions = sum(1 for a, b in zip(seq1, seq2) if a != '-' or b != '-')

    return matches / non_gap_positions if non_gap_positions > 0 else 0.0


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <a3m_file>", file=sys.stderr)
        sys.exit(1)

    a3m_path = sys.argv[1]
    with open(a3m_path) as f:
        sequences = parse_a3m(f.read())

    n = len(sequences)
    if n < 2:
        print(f"Sequences: {n}")
        print("Not enough sequences to compute identity.")
        return

    query_seq = sequences[0].seq
    identities = sorted(calculate_sequence_identity(query_seq, s.seq) for s in sequences[1:])

    median = identities[len(identities) // 2]

    print(f"Query:    {sequences[0].hdr}")
    print(f"Hits:     {len(identities)}")
    print(f"Identity: mean={sum(identities)/len(identities):.2%}"
          f"  median={median:.2%}"
          f"  min={identities[0]:.2%}"
          f"  max={identities[-1]:.2%}")

    print("Hits above threshold:")
    for pct in range(10, 100, 10):
        t = pct / 100
        count = sum(1 for x in identities if x >= t)
        print(f"  >={pct:2d}%: {count}")



if __name__ == "__main__":
    main()
