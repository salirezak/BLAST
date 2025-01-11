def read_fasta(file_path: str) -> dict:
    sequences = {}
    with open(file_path, "r") as f:
        lines = filter(bool, map(lambda l: l.strip(), f.readlines()))

    for l in lines:
        if l.startswith(">"):
            sequences[key := l] = ""
        else:
            sequences[key] += l

    return sequences