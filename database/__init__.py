import pandas as pd
from .loader import load_database
from .generator import generate_database, load_fasta, save_database



def get_database(fasta_file_path: str,k: int, threshold: int, scoring_matrix: pd.DataFrame, scoring_matrix_name: str):
    pickle_file_path = f"database/K{k}_T{threshold}_{scoring_matrix_name}.pkl"
    try:
        return load_database(pickle_file_path)
    except FileNotFoundError:
        database_fasta = load_fasta(fasta_file_path)
        neighborhoods = generate_database(database_fasta, k, threshold, scoring_matrix)
        save_database(pickle_file_path, neighborhoods)
        return neighborhoods