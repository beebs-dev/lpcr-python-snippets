import random
from primer3.thermoanalysis import ThermoAnalysis

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

def calc_dg(seq1, seq2):
    analysis = ThermoAnalysis(mv_conc=50,
                              dv_conc=3,
                              dntp_conc=0.8,
                              dna_conc=200,
                              temp_c=70) # Note: should be 25
    result = analysis.calc_heterodimer(seq1, seq2)
    return result.dg / 1000  # Convert to kcal/mol

def gen_oligos_with_prelude(sequence, 
                            *args, 
                            prelude_length=24, 
                            prelude_sticky_length=14, 
                            **kwargs):
    # "Prelude" sequence is first N base pairs.
    prelude = sequence[:prelude_length]
    # Shorten the sequence passed to gen_oligos to only
    # include the sticky end of the prelude oligo.
    shortened = sequence[prelude_length-prelude_sticky_length:]
    oligos, sticky_lengths = gen_oligos(shortened, *args, **kwargs)
    # Add the prelude oligo to the beginning of the list.
    oligos = [prelude] + oligos
    sticky_lengths = [prelude_sticky_length] + sticky_lengths
    return oligos, sticky_lengths

def gen_oligos(sequence,
               oligo_length=60,
               sticky_length=14,
               sticky_lengths=None,
               reverse_oligo_length=None,
               reverse_oligo_lengths=None,
               max_dg=None):
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

    # Keep track of the sticky lengths for each overlap
    sticky_lengths = [] if sticky_lengths is None else sticky_lengths

    # Keep track of how long each reverse strand is.
    reverse_oligo_lengths = [] if reverse_oligo_lengths is None else reverse_oligo_lengths
    
    i = 0
    while len(remaining) > 0:
        # Determine sticky end lengths. If no value is provided, use the default.
        # If the delta G for the sticky end is not low enough, the length will
        # be increased and this loop will be retried.
        if i < len(sticky_lengths):
            stick_len = sticky_lengths[i]
        else:
            # Use the default sticky length.
            stick_len = sticky_length
            # Keep track of the sticky length for this oligo.
            sticky_lengths.append(stick_len)

        # Compute the length of the current oligo's "body" (non-sticky part).
        # Try to use the full oligo length, but if there are not enough bases left,
        # use as many as possible.
        if reverse_oligo_length is not None and len(result) % 2 == 1:
            # This is a reverse strand. Limit its length.
            # There is an issue where the second oligo binding to a reverse
            # strand will grow its sticky end to overlap with the reverse strand's
            # sticky end to the first oligo. In this case, the reverse strand needs
            # to be lengthened to accomodate both sticky ends without overlap.
            #reverse_oligo_lengths[len(reverse_oligo_lengths)-1] = reverse_oligo_length
            i_rv = i // 2
            if i_rv < len(reverse_oligo_lengths):
                # Use the length of the reverse oligo from the list.
                body_len = reverse_oligo_lengths[i_rv]
            else:
                # Use the default reverse oligo length.
                body_len = reverse_oligo_length
                reverse_oligo_lengths.append(body_len)
        else:
            body_len = oligo_length

        # Minimum body length is twice the minimum (default) sticky length.
        # This is to ensure that the oligo is long enough to accommodate
        # sticky ends on both sides.
        next_stick_len = sticky_lengths[i + 1] if i + 1 < len(sticky_lengths) else sticky_length
        body_len = max(body_len, stick_len + next_stick_len)

        # If there are not enough bases left, use as many as possible.
        body_len = min(body_len - stick_len, len(remaining))

        # Extract the body of the current oligo.
        current_oligo = remaining[:body_len]

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
            # Current oligo is forward strand (Case 1).
            source = last_oligo[:stick_len]
            sticky = reverse_complement(source)
            # Sticky end is 5' on current oligo.
            current_oligo = sticky + current_oligo
        else:
            # Current oligo is backward strand (Case 2).
            source = last_oligo[-stick_len:]
            sticky = reverse_complement(source)
            # Sticky end is 3' on current oligo.
            # Since this current strand is template side, we also
            # need to reverse complement the body of the oligo.
            current_oligo = reverse_complement(current_oligo) + sticky

        # Check if the oligo needs extension.
        dg = calc_dg(source, sticky)
        if max_dg is not None and dg > max_dg:
            # Try again with a longer sticky length.
            sticky_lengths[i] += 1

            # Is the current oligo a forward strand?
            if i % 2 == 1:
                rv_i = i // 2
                rev_strand_length = reverse_oligo_lengths[rv_i]
                prev_stick_len = sticky_lengths[i - 1]
                if rev_strand_length < stick_len + prev_stick_len:
                    # Increase the reverse strand length to accomodate
                    reverse_oligo_lengths[rv_i] += 1
                    # Recurse so we rebuild ALL oligos.
                    return gen_oligos(sequence,
                                      oligo_length=oligo_length,
                                      sticky_length=sticky_length,
                                      sticky_lengths=sticky_lengths,
                                      reverse_oligo_length=reverse_oligo_length,
                                      reverse_oligo_lengths=reverse_oligo_lengths,
                                      max_dg=max_dg)
            continue

        # Consume the bases that were just used for this oligo's body.
        remaining = remaining[body_len:]

        # Add the completed oligo to the result.
        result.append(current_oligo)
        i += 1

    # All done! Return the list of oligos.
    return result, sticky_lengths

def gen_random_seq(length):
    """Generate a random DNA sequence of given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def test_gen_oligos():
    """Test the gen_oligos function."""
    seq_length = 100 + oligo_length + random.randint(0, 1000)
    sequence = gen_random_seq(seq_length)
    sticky_lengths = [16, 14, 20, 15]
    oligos = gen_oligos(sequence,
                        oligo_length=oligo_length,
                        sticky_length=sticky_length,
                        sticky_lengths=sticky_lengths)
    assert len(oligos) > 0, "No oligos generated."
    for i, oligo in enumerate(oligos):
        assert len(oligo) <= oligo_length, f"Invalid oligo {i} length: {len(oligo)} != {oligo_length}"
        assert all(base in 'ACGT' for base in oligo), f"Invalid base in oligo: {oligo}"
    for i in range(1, len(oligos)):
        stick_length = sticky_lengths[i-1] if sticky_lengths is not None and i-1 < len(sticky_lengths) else sticky_length
        # Note: all oligos are ALWAYS 5' to 3'
        if i % 2 == 0:
            # Current oligo is coding strand.
            # Current sticky end is at the 5' end.
            # Sticky end for last oligo is also at its 5' end.
            last_sticky = oligos[i-1][:stick_length]
            current_sticky = oligos[i][:stick_length]
        else:
            # Current oligo is template strand.
            # Current sticky end is at the 3' end.
            # Sticky end for last oligo is also at its 3' end.
            last_sticky = oligos[i-1][-stick_length:]
            current_sticky = oligos[i][-stick_length:]
        # Check that the sticky ends are reverse complements
        expected_sticky = reverse_complement(current_sticky)
        assert last_sticky == expected_sticky, \
            f"Sticky ends do not match: {last_sticky} != {expected_sticky}"

if __name__ == "__main__":
    num_tests = 10000
    for i in range(num_tests):
        test_gen_oligos()
    print(f"All {num_tests} tests passed.")
