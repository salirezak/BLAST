from itertools import product
from functools import lru_cache
import pandas as pd



@lru_cache(maxsize=None)
def generate_kmers(seq: str, k: int) -> list:
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]

def calculate_alignment_score(seq1: str, seq2: str, scoring_matrix: pd.DataFrame) -> int:
    return sum(scoring_matrix.loc[seq1[i], seq2[i]] for i in range(len(seq1)))

def generate_neighbors_map_filter(kmer: str, threshold: int, scoring_matrix: pd.DataFrame) -> set:
    neighbors = map(lambda x: "".join(x), product(scoring_matrix.columns, repeat=len(kmer)))
    neighbors = filter(lambda x: calculate_alignment_score(kmer, x, scoring_matrix) >= threshold, neighbors)
    return set(neighbors)

def generate_neighbors_recursive(kmer: str, threshold: int, scoring_matrix: pd.DataFrame) -> set:
    neighbors = set()
    def recurse(current, index, partial_score):
        if index == len(kmer):
            if partial_score >= threshold:
                neighbors.add(current)
            return
        for residue in scoring_matrix.columns:
            new_score = partial_score + scoring_matrix.loc[kmer[index], residue]
            if new_score + (len(kmer) - index - 1) * max(scoring_matrix.max()) >= threshold:
                recurse(current + residue, index + 1, new_score)
    recurse("", 0, 0)
    return neighbors

def generate_neighbors_stack(kmer: str, threshold: int, scoring_matrix: pd.DataFrame) -> set:
    neighbors = set()
    stack = [("", 0, 0)]  # (current, index, partial_score)
    while stack:
        current, index, partial_score = stack.pop()
        if index == len(kmer):
            if partial_score >= threshold:
                neighbors.add(current)
            continue
        for residue in scoring_matrix.columns:
            new_score = partial_score + scoring_matrix.loc[kmer[index], residue]
            if new_score + (len(kmer) - index - 1) * max(scoring_matrix.max()) >= threshold:
                stack.append((current + residue, index + 1, new_score))
    return neighbors

generate_neighbors = generate_neighbors_map_filter

def generate_neighborhood(seq: str, k: int, threshold: int, scoring_matrix: pd.DataFrame) -> list:
    return [generate_neighbors(kmer, threshold, scoring_matrix) for kmer in generate_kmers(seq, k)]
