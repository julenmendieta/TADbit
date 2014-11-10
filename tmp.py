from pytadbit.parsers.hic_parser import read_matrix
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import log2
import numpy as np

a = read_matrix('chrXX.tsv', scale=40000)[0]
b = a.get_section('chrT')
c = a.get_section('chrU')

a.write('lala.pik', gzip=False)
a.write('lala.pik.gzip', gzip=True)
a.write('lala.hmtx', format='mtx', gzip=False)
b.write('lala.mtx.gzip', format='mtx', gzip=True)
a.write('lala.mtx', format='mtx', headers=False, gzip=False)
c.write('lala.abc', format='abc', gzip=False)
a.write('lala.abc', format='abc', gzip=False)
a.write('lala.abc.gzip', format='abc', gzip=True)

d = a.get_sample((10,50))

e = a.get_scaled(90000)

log = True
fig, subaxes = plt.subplots(nrows=1, ncols=2)
axeslist = subaxes.flatten()
if log:
    axeslist[0].imshow(np.log2([[a[i * a._size + j] for j in xrange(a._size)]
                        for i in xrange(a._size)]), interpolation='none',
                       cmap='hot')
    im = axeslist[1].imshow(np.log2([[e[i * e._size + j] for j in xrange(e._size)]
                             for i in xrange(e._size)]), interpolation='none',
                            cmap='hot')
else:
    axeslist[0].imshow([[a[i * a._size + j] for j in xrange(a._size)]
                        for i in xrange(a._size)], interpolation='none',
                       cmap='hot')
    im = axeslist[1].imshow([[e[i * e._size + j] for j in xrange(e._size)]
                             for i in xrange(e._size)], interpolation='none',
                            cmap='hot')
fig.subplots_adjust(bottom=.2)
data = [a[i * a._size + j] for i in xrange(a._size) for j in xrange(a._size)]
gradient = np.linspace(0, 1, max(data))
gradient = np.vstack((gradient, gradient))
axe = fig.add_axes([0.1, 0.05, 0.5, 0.2])
if not log:
    cb  = axe.imshow(np.log2(gradient), aspect='auto', cmap='hot')
    h  = axe.hist(np.log2([i or 0.5 for i in data]), color='white',
                  normed=True, bins=20, histtype='step')
else:
    cb  = axe.imshow(gradient, aspect='auto', cmap='hot')
    h  = axe.hist(data, color='blue',
                  normed=True, bins=100, histtype='step')
axe.set_xlim((0, max(data)))
axe.set_ylim((0, max(h[0])))
plt.hist(data)
plt.show()
