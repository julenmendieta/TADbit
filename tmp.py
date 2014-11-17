from pytadbit.parsers.hic_parser import read_matrix
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import log2
from scipy.linalg import eigh
import numpy as np
from matplotlib.ticker import NullFormatter
from scipy.stats import norm as sc_norm, skew, kurtosis

def view_mat(a):
    cmap = 'afmhot'
    fig = plt.figure(figsize=(10.,9.1))
    evals, evect = eigh([[a[i * a._size + j] for j in xrange(a._size)] for i in xrange(a._size)])
    axe = fig.add_subplot(111)
    axe.imshow(np.log2([[a[i * a._size + j] for j in xrange(a._size)]
                        for i in xrange(a._size)]), interpolation='none',
               cmap=cmap)
    fig.subplots_adjust(top=.70, left=0.35)
    data = [a[i * a._size + j] for i in xrange(a._size) for j in xrange(a._size)]
    data = np.log2([(i or 0.1) for i in data])
    gradient = np.linspace(0, max(data), len(data))
    gradient = np.vstack((gradient, gradient))
    axe2 = fig.add_axes([0.1, 0.78, 0.2, 0.1])
    h  = axe2.hist(data, color='grey', bins=20, histtype='step', normed=True)
    cb  = axe2.imshow(gradient, aspect='auto', cmap=cmap, extent=(0, max(data) , 0, max(h[0])))
    axe2.set_xlim((0, max(data)))
    axe2.set_ylim((0, max(h[0])))
    axe3 = fig.add_axes([0.35, 0.71, 0.55, 0.1])
    axe3.vlines(range(a._size), 0, evect[0], color='k')
    axe3.hlines(0, 0, a._size, color='red')
    axe3.xaxis.set_major_formatter(NullFormatter())
    axe.set_xlim((-0.5, a._size - .5))
    axe.set_ylim((-0.5, a._size - .5))
    axe3.set_xlim((-0.5, a._size - .5))
    axe2.set_xlabel('log interaction count')
    normfit = sc_norm.pdf(data, np.mean(data), np.std(data))
    axe2.plot(data, normfit, 'g.', markersize=1, alpha=.1)
    axe2.set_title('skew: %.3f, kurtosis: %.3f' % (skew(data), kurtosis(data)))

    for lenk in range(1, len(max(a.section_sizes.keys(), key=len)) + 1):
        newax = axe.twinx()
        newax.set_frame_on(True)
        newax.patch.set_visible(False)
        newax.yaxis.set_ticks_position('left')
        newax.yaxis.set_label_position('left')
        newax.spines['left'].set_position(('outward', 40 * (1 + len(max(a.section_sizes.keys(), key=len)) - lenk)))
        vals = []
        keys = []
        total = 0
        for i in [s for s in a.ordered_sections if len(s) == lenk]:
            print i
            vals.append(a.section_sizes[i] + total)
            keys.append(i[-1])
            total += a.section_sizes[i]
        newax.set_yticks(vals)
        newax.set_yticklabels(keys)

    plt.show()


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

# d = a.get_sample((10,50))

e = a.get_scaled(90000)
print '\n'.join(['\t'.join([str(a[i, j]) for j in xrange(a._size)]) for i in xrange(a._size)])
print '-'*80
print '\n'.join(['\t'.join([str(e[i, j]) for j in xrange(e._size)]) for i in xrange(e._size)])

view_mat(e)


view_mat(a)
