import os

# Key-value flags for EMBOSS Needle command line
needle_opts = dict(
    datafile='EDNAFULL',
    gapopen=10,
    gapextend=0.5,
    endopen=10.0,
    endextend=0.5,
    aformat3='pair',
)

# Extra flags for EMBOSS Needle command line
needle_extra = ['-auto', '-snucleotide1', '-snucleotide2']

# Runs EMBOSS Needle to perform an alignment between the two sequences.
def run_needle(aseq, bseq, remove_comments=True):
    afn = '/tmp/in.aseq'
    with open(afn, 'w') as f:
        f.write(aseq)
    bfn = '/tmp/in.bseq'
    with open(bfn, 'w') as f:
        f.write(bseq)
    ofn = '/tmp/out.needle'
    cmd = f'needle -asequence {afn} -bsequence {bfn} -outfile {ofn}'
    for k, v in needle_opts.items():
        cmd += f' -{k} {v}'
    for arg in needle_extra:
        cmd += f' {arg}'
    try:
        code = os.system(cmd)
        if code != 0:
            raise Exception(f'EMBOSS Needle error exit code: {code}')
        with open(ofn, 'r') as f:
            content = f.read()
    except:
        raise Exception('Error running EMBOSS Needle')
    finally:
        os.remove(afn)
        os.remove(bfn)
        os.remove(ofn)
    if remove_comments:
        lines = content.splitlines()
        content = []
        for line in lines:
            if line.startswith('#'):
                continue
            content.append(line)
        content = '\n'.join(content)
    return content
