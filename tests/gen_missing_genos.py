import random
import numpy as np

genolines = []
headers = []
with open('../data/snps.132k.clean.noX.ped', 'r') as f:
    for line in f:
        l = line.strip().split(' ')[6:]
        headers.append(line.strip().split(' ')[:6])
        line_pairs = zip(l[::2], l[1::2])
        k = random.randint(0, len(line_pairs)/20)
        indices = random.sample(range(len(line_pairs)), k)
        for idx in indices:
            line_pairs[idx] = ('0', '0')
        genolines.append(line_pairs)

with open('../data/snps.132k.clean.noX.missing_genos.ped', 'w') as w:
    for i in range(len(genolines)):
        for head in headers[i]:
            w.write(head + ' ')
        for geno in genolines[i]:
            w.write(geno[0] + ' ' + geno[1] + ' ')
        w.write('\n')
