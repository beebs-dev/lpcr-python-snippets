import os

from primer3.thermoanalysis import ThermoAnalysis

from gen_oligos_v2 import (
    calc_dg,
    gen_oligos_with_prelude,
    gen_random_seq,
    reverse_complement,
)
from needle import needle_extra, needle_opts, run_needle
from utils import read_dna_seq

__VERSION__ = '0.1.0'
print('beebs.dev Linking PCR Design Tool')
print('Version: {}'.format(__VERSION__))
print('')

sticky_length = 12
reverse_oligo_length = 36
oligo_length = 60
max_dg = -4.0
prelude_length = 24
prelude_sticky_length = 14

#a = gen_random_seq(300)
a = read_dna_seq('genes/E130003G02Rik.txt')

print('Settings:')
print(f'  Default sticky end length: {sticky_length}')
print(f'  Default reverse oligo length: {reverse_oligo_length}')
print(f'  Oligo length limit: {oligo_length}')
print(f'  Max delta G: {max_dg} kcal/mol')
print(f'  Prelude length (overall): {prelude_length}')
print(f'  Prelude sticky length: {prelude_sticky_length}')
print('')
print(f'Input sequence ({len(a)} base pairs):')
print('5\'-{}-3\''.format(a))
print('')

oligos, sticky_lengths = gen_oligos_with_prelude(
    a,
    oligo_length=oligo_length,
    sticky_length=sticky_length,
    reverse_oligo_length=reverse_oligo_length,
    max_dg=max_dg,
    prelude_length=prelude_length,
    prelude_sticky_length=prelude_sticky_length)

print(f'Generated {len(oligos)} oligonucleotides:')
print('')
"""
5'------------3'     5'------------3'
          3'-------------5'
"""

# Extract the prelude oligo and sticky length.
# Print these out separately because both it and the next one are 5' to 3'.
prelude = oligos[0]
oligos = oligos[1:]
prelude_sticky_length = sticky_lengths[0]
sticky_lengths = sticky_lengths[1:]

prelude_sticky = prelude[-prelude_sticky_length:]
prelude_dg = calc_dg(prelude_sticky, reverse_complement(prelude_sticky))

spaces = ''.join(' ' for i in range(prelude_length - prelude_sticky_length))
bars = ''.join('|' for i in range(prelude_sticky_length))
print('Prelude (initial) oligo:')
print(f"[P] 5'-{prelude}-3'")
print(f"{spaces}       {bars}")
print(f"[0] {spaces}5'-{oligos[0]}-3'")
print(f'# sticky base pairs: {prelude_sticky_length}')
print('deltaG: {:.2f} kcal/mol'.format(prelude_dg))
print('')
print('')

for i in range(1, len(oligos)):
    sticky_length = sticky_lengths[i - 1]
    previous = oligos[i - 1]
    current = oligos[i]
    num_spaces = len(previous) - sticky_length
    spaces = ''.join(' ' for i in range(num_spaces))
    bars = ''.join('|' for i in range(sticky_length))
    if i % 2 == 0:
        # Previous oligo is template strand.
        print(f"[{i-1}] 3'-{previous[::-1]}-5'")
        print(f"{spaces}       {bars}")
        print(f"[{i}] {spaces}5'-{current}-3'")
        # Sticky end is at the 5' end of the previous oligo.
        source = previous[:sticky_length]
        target = current[:sticky_length]
    else:
        # Previous oligo is coding strand.
        print(f"[{i-1}] 5'-{previous}-3'")
        print(f"{spaces}       {bars}")
        print(f"[{i}] {spaces}3'-{current[::-1]}-5'")
        # Sticky end is at the 3' end of the previous oligo.
        source = previous[-sticky_length:]
        target = current[-sticky_length:]
    dg = calc_dg(source, target)
    print(f'# sticky base pairs: {sticky_length}')
    print('deltaG: {:.2f} kcal/mol'.format(dg))
    print("")
    print("")


def align_oligos(oligo1, oligo2):
    alignment = run_needle(oligo1, oligo2)
    print(alignment)
    print('')


def align_oligos_both(i, j, oligo1, oligo2):
    print('Forward alignments:')
    print(f'[{i}] 5\'-{oligo1}-3\'')
    print(f'[{j}] 5\'-{oligo2}-3\'')
    align_oligos(oligo1, oligo2)

    print('Backward alignments:')
    oligo2 = oligo2[::-1]
    print(f'[{i}] 5\'-{oligo1}-3\'')
    print(f'[{j}] 3\'-{oligo2}-5\'')
    align_oligos(oligo1, oligo2)


def run_seqalign(oligos):
    oligo1 = prelude
    for i in range(len(oligos)):
        oligo2 = oligos[i]
        align_oligos_both('P', i, oligo1, oligo2)
        print(''.join('=' for _ in range(80)))
        print('')

    for i in range(len(oligos)):
        oligo1 = oligos[i]
        for j in range(i + 1, len(oligos)):
            oligo2 = oligos[j]
            align_oligos_both(i, j, oligo1, oligo2)
            print(''.join('=' for _ in range(80)))
            print('')


def print_needle_args():
    for k, v in needle_opts.items():
        print(f'  -{k} {v}')
    for arg in needle_extra:
        print(f'  {arg}')


def write_oligos(oligos):
    with open('example_oligos.txt', 'w') as f:
        f.write(prelude + '\n')
        for oligo in oligos:
            f.write(oligo + '\n')


print('Running sequence alignment...')
print(
    'Each input strand is fed into EMBOSS Needle as it appears with the following arguments:'
)
print_needle_args()
print(
    'These match the default values at https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle'
)
print('')
run_seqalign(oligos)
write_oligos(oligos)

print('Done.')
print('')
