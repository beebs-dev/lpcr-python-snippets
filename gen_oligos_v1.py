import random

oligo_length = 60
sticky_length = 14

def complement(seq):
    """Return the complement of a DNA sequence."""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(mapping[base] for base in seq)

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(mapping[base] for base in reversed(seq))

def gen_oligos(sequence, oligo_length, sticky_length):
    """Design linking PCR oligonucleotides from a given sequence."""

    # Normalize the sequence to upper case
    sequence = sequence.upper()

    # Check if it's even necessary to make oligos
    if len(sequence) <= oligo_length:
        return [sequence]

    # Check if the sequence is valid
    if not all(base in 'ACGT' for base in sequence):
        raise ValueError("Invalid sequence: only A, C, G, T are allowed.")

    # Check if the oligo length is valid
    if oligo_length <= 0 or oligo_length > len(sequence):
        raise ValueError("Invalid oligo length: must be between 1 and the length of the sequence.")

    # Check if the sticky length is valid
    if sticky_length < 0 or sticky_length >= oligo_length:
        raise ValueError("Invalid sticky length: must be between 0 and oligo length - 1.")
    
    # Keep track of remaining bases in the sequence.
    remaining = sequence

    # First oligo will be full length (special case)
    current_oligo = remaining[:oligo_length]
    
    # Consume the bases used for this first oligo.
    remaining = remaining[oligo_length:]

    # Keep track of the list of completed oligos.
    result = [current_oligo]

    while len(remaining) > 0:        
        # Compute the length of the current oligo's "body" (non-sticky part).
        # Try to use the full oligo length, but if there are not enough bases left,
        # use as many as possible.
        body_len = min(oligo_length - sticky_length, len(remaining))

        # Extract the body of the current oligo.
        current_oligo = remaining[:body_len]
        
        # Consume the bases that were just used for this oligo's body.
        remaining = remaining[body_len:]

        # For the sticky end: a reference to the previously built oligo.        
        last_oligo = result[-1]

        # Note: output oligos will ALWAYS be 5' to 3'
        # We have two cases here:
        #     (Case 1) for 5' to 3' oligos (i is even):
        #                              9876543210
        #     Prev: 3'-.........AAAAAAAAAAAAAAAAA-5' (sticky is 5' end of last oligo)
        #     Next:          5'-TTTTTTTTTTTTTTTTT... (sticky is 5' end of cur oligo)
        #                       0123456789    
        #
        #     (Case 2) for 3' to 5' oligos (i is odd):
        #                             0123456789N
        #     Prev: 5'-.........AAAAAAAAAAAAAAAAA-3' (sticky is 3' end of last oligo)
        #     Next:          3'-TTTTTTTTTTTTTTTTT... (sticky is 3' end of cur oligo)
        #                       N9876543210
        if len(result) % 2 == 0:
            # Current oligo is coding strand (Case 1).
            sticky = reverse_complement(last_oligo[:sticky_length])
            # Sticky end is 5' on current oligo.
            current_oligo = sticky + current_oligo
        else:
            # Current oligo is template strand (Case 2).
            sticky = reverse_complement(last_oligo[-sticky_length:])
            # Sticky end is 3' on current oligo.
            # Since this current strand is template side, we also
            # need to reverse complement the body of the oligo.
            current_oligo = reverse_complement(current_oligo) + sticky

        # Add the completed oligo to the result.
        result.append(current_oligo)

    # All done! Return the list of oligos.
    return result

def gen_random_seq(length):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def test_gen_oligos():
    """Test the gen_oligos function."""
    seq_length = 100 + oligo_length + random.randint(0, 1000)
    sequence = gen_random_seq(seq_length)
    oligos = gen_oligos(sequence, oligo_length, sticky_length)
    assert len(oligos) > 0, "No oligos generated."
    for i, oligo in enumerate(oligos):
        assert len(oligo) <= oligo_length, f"Invalid oligo {i} length: {len(oligo)} != {oligo_length}"
        assert all(base in 'ACGT' for base in oligo), f"Invalid base in oligo: {oligo}"
    for i in range(1, len(oligos)):
        # Note: all oligos are ALWAYS 5' to 3'
        if i % 2 == 0:
            # Current oligo is coding strand.
            # Current sticky end is at the 5' end.
            # Sticky end for last oligo is also at its 5' end.
            last_sticky = oligos[i-1][:sticky_length]
            current_sticky = oligos[i][:sticky_length]
        else:
            # Current oligo is template strand.
            # Current sticky end is at the 3' end.
            # Sticky end for last oligo is also at its 3' end.
            last_sticky = oligos[i-1][-sticky_length:]
            current_sticky = oligos[i][-sticky_length:]
        # Check that the sticky ends are reverse complements
        expected_sticky = reverse_complement(current_sticky)
        assert last_sticky == expected_sticky, \
            f"Sticky ends do not match: {last_sticky} != {expected_sticky}"

if __name__ == "__main__":
    num_tests = 10000
    for i in range(num_tests):
        test_gen_oligos()
    print(f"All {num_tests} tests passed.")
