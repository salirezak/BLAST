from pickle import load



def load_database(file_path: str) -> dict:
    with open(file_path, "rb") as f:
        neighborhoods = load(f)
        return neighborhoods