"""
18 Sep 2014
"""
from warnings import warn
from pytadbit.utils.extraviews import tadbit_savefig
import numpy as np
try:
    from matplotlib.ticker import MultipleLocator
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    warn('matplotlib not found\n')

def tad_density_plot(experiment, axe=None, focus=None, extras=None,
                     normalized=True, savefig=None, shape='ellipse'):
    """
    Draw an summary of the TAD found in a given experiment and their density
    in terms of relative Hi-C interaction count.

    :param experiment: Experiment to visualize
    :param None focus: can pass a tuple (bin_start, bin_stop) to display the
       alignment between these genomic bins
    :param None extras: list of coordinates (genomic bin) where to draw a
       red cross
    :param None ymax: limit the y axis up to a given value
    :param ('grey', ): successive colors for alignment
    :param True normalized: normalized Hi-C count are plotted instead of raw
       data.
    :param 'ellipse' shape: which kind of shape to use as schematic
       representation of TADs. Implemented: 'ellipse', 'rectangle',
       'triangle'
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    """
    if not experiment.tads:
        raise Exception("TAD borders not found\n")
    _tad_density_plot(experiment, axe=axe, focus=focus,
                      extras=extras, normalized=normalized,
                      savefig=savefig, shape=shape)


def _tad_density_plot(xpr, maxys=None, fact_res=1., axe=None,
                     focus=None, extras=None, normalized=True,
                     savefig=None, shape='ellipse'):
    """
    """
    from matplotlib.cm import jet
    siz = xpr.size
    show=False
    if focus:
        figsiz = 4 + (focus[1] - focus[0])/30
    else:
        figsiz = 4 + (siz)/30

    if not axe:
        fig = plt.figure(figsize=(figsiz, 1 + 1 * 1.8))
        axe = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0)
        show=True

    zsin = np.sin(np.linspace(0, np.pi))
    
    shapes = {'ellipse'   : lambda h : [0] + list(h * zsin) + [0],
              'rectangle' : lambda h: [0] + [h] * 50 + [0],
              'triangle'  : lambda h: ([h/25 * i for i in xrange(26)] +
                                       [h/25 * i for i in xrange(25, -1, -1)])}
    
    try:
        shape = shapes[shape]
    except KeyError:
        import this
        table = ''.join([this.d.get(chr(i), chr(i)) for i in range(256)])
        if locals()['funcr'.translate(table)].translate(table) == ''.join(
            [this.s[i].upper() if this.s[i-1] is 'v' else this.s[i]
             for i in [24, 36, 163, 8, 6, 16, 36]]):
            shape = lambda h: (
                [h/25 * i for i in xrange(25)] + [h+0.2]*2 +
                [h/25 * i for i in xrange(24, -1, -1)])
        else:
            raise NotImplementedError(
                '%s not valid, use one of ellipse, rectangle or triangle')
    maxys = maxys if isinstance(maxys, list) else []
    zeros = xpr._zeros or {}
    if normalized and xpr.norm:
        norms = xpr.norm[0]
    elif xpr.hic_data:
        if normalized:
            warn("WARNING: weights not available, using raw data")
        norms = xpr.hic_data[0]
    else:
        warn("WARNING: raw Hi-C data not available, " +
             "TAD's height fixed to 1")
        norms = None
    diags = []
    siz = xpr.size
    sp1 = siz + 1
    if norms:
        for k in xrange(1, siz):
            s_k = siz * k
            diags.append(sum([norms[i * sp1 + s_k]
                             if not (i in zeros
                                     or (i + k) in zeros) else 0.
                              for i in xrange(siz - k)]) / (siz - k))
    for tad in xpr.tads:
        start, end = (int(xpr.tads[tad]['start']) + 1,
                      int(xpr.tads[tad]['end']) + 1)
        if norms:
            matrix = sum([norms[i + siz * j]
                         if not (i in zeros
                                 or j in zeros) else 0.
                          for i in xrange(start - 1, end - 1)
                          for j in xrange(i + 1, end - 1)])
        try:
            if norms:
                height = float(matrix) / sum(
                    [diags[i-1] * (end - start - i)
                     for i in xrange(1, end - start)])
            else:
                height = 1.
        except ZeroDivisionError:
            height = 0.
        maxys.append(height)
        start = float(start) / fact_res  # facts[iex]
        end   = float(end) / fact_res  # facts[iex]
        axe.fill([start] + list(np.linspace(start, end)) + [end], shape(height),
                 alpha=.8 if height > 1 else 0.4,
                 facecolor='grey', edgecolor='grey')
    if extras:
        axe.plot(extras, [.5 for _ in xrange(len(extras))], 'rx')
    axe.grid()
    axe.patch.set_visible(False)
    axe.set_ylabel('Relative\nHi-C count')
    #
    for tad in xpr.tads:
        if not xpr.tads[tad]['end']:
            continue
        tad = xpr.tads[tad]
        axe.plot(((tad['end'] + 1.) / fact_res, ), (0, ),
                 color=jet(tad['score'] / 10) if tad['score'] else 'w',
                 mec='none' if tad['score'] else 'k', marker=6, ms=9, alpha=1,
                 clip_on=False)
    axe.set_xticks([1] + range(100, int(tad['end'] + 1), 50))
    axe.minorticks_on()
    axe.xaxis.set_minor_locator(MultipleLocator(10))
    axe.hlines(1, 1, end, 'k', lw=1.5)
    if show:
        tit1 = fig.suptitle("TAD borders", size='x-large')
        plt.subplots_adjust(top=0.76)
        fig.set_facecolor('white')
        plots = []
        for scr in xrange(1, 11):
            plots += plt.plot((100,),(100,), marker=6, ms=9,
                              color=jet(float(scr) / 10), mec='none')
        axe.legend(plots,
                   [str(scr) for scr in xrange(1, 11)],
                   numpoints=1, title='Border scores',
                   fontsize='small', loc='lower left',
                   bbox_to_anchor=(1, 0.1))
        axe.set_ylim((0, max(maxys) + 0.4))
        if savefig:
            tadbit_savefig(savefig)
        else:
            plt.show()


