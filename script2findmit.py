""" 
Script to find the smallest 30 scaffolds in the reference genome, which may be candidates for being the mitochondrial genome.
"""

from collections import OrderedDict

fa = "/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/reference/PitayaGenomic.fa"
lengths = OrderedDict()
name = None
with open(fa) as f:
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            lengths[name] = 0
        else:
            lengths[name] += len(line.strip())

small = sorted(lengths.items(), key=lambda x: x[1])[:30]
for k,v in small:
    print(k, v)
