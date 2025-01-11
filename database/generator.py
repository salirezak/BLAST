from pickle import dump
from fasta import read_fasta
from utils import generate_neighborhood
import pandas as pd



def load_fasta(file_path: str) -> dict:
    return read_fasta(file_path)

def generate_database(database_fasta: dict, k: int, threshold: int, scoring_matrix: pd.DataFrame) -> dict:
    neighborhoods = {}
    for name, seq in database_fasta.items():
        name = name[1:]
        print("generating neighborhood for", name)
        neighborhood = generate_neighborhood(seq, k, threshold, scoring_matrix)
        neighborhoods[name] = (seq, neighborhood)
    return neighborhoods

def save_database(file_path: str, neighborhoods: dict) -> None:
    with open(file_path, "wb") as f:
        dump(neighborhoods, f)