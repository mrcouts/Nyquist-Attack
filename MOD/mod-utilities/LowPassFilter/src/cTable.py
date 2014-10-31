#!/usr/bin/env python2

import sys

from mpmath import mp
mp.dps = 50
    
N = long(sys.argv[1]) if len(sys.argv) > 1 else 100000

inicio = 0
fim = 1*mp.pi
dx = (fim-inicio)/(N-1)

Idx = mp.mpf(1/dx)

with open('C.h', 'w') as f:

    f.write("\n".join([
        'using namespace std;',
        '',
        '#define C_N {N}',
        '#define C_inicio {inicio}',
        '#define C_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const double _c_[] = {')

    for i in range(N):
        if i == 0:
            value = 2
        elif i== N-1:
            value = 0
        else:
            value = (i*dx + inicio)/mp.tan( (i*dx + inicio)/2 )
        v = "%.50e" % value
        f.write( v + ' ,\n')

    f.write('\n};')