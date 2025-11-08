def read_dna_seq(filename):
    """
    Read a DNA sequence from a file.
    """
    with open(filename, 'r') as f:
        content = f.read()
    return content.strip()
