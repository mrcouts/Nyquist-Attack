#!/usr/bin/env python2

import sys

from mpmath import mp
mp.dps = 50
    
N = long(sys.argv[1]) if len(sys.argv) > 1 else 1000000

inicio = -1*mp.pi
fim = 1*mp.pi
dx = (fim-inicio)/(N-1)

Idx = mp.mpf(1/dx)

with open('Cosseno.h', 'w') as f:

    f.write("\n".join([
        'using namespace std;',
        '',
        '#define COS_N {N}',
        '#define COS_inicio {inicio}',
        '#define COS_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const double Cosseno[] = {')

    for i in range(N):
        value = mp.cos(i*dx + inicio)
        v = "%.50e" % value
        f.write( v + ' ,\n')

    f.write('\n};')