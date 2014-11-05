"""
17 Sep 2014


"""

from math                              import sqrt
from pytadbit.utils.extraviews         import tadbit_savefig
from sys                               import stderr

try:
    import matplotlib.pyplot as plt
except ImportError:
    stderr.write('matplotlib not found\n')

def interaction_map(experiments, tad=None, focus=None, paint_tads=False,
                    axe=None, show=True, logarithm=True, normalized=False,
                    relative=True, decorate=True, savefig=None, clim=None,
                    scale=(8, 6)):
    """
    Visualize the matrix of Hi-C interactions of a given experiment

    :param None names: name of the experiment to visualize, or list of
       experiment names. If None, all experiments will be shown
    :param None tad: a given TAD in the form:
       ::

         {'start': start,
          'end'  : end,
          'brk'  : end,
          'score': score}

       **Alternatively** a list of the TADs can be passed (all the TADs
       between the first and last one passed will be showed. Thus, passing
       more than two TADs might be superfluous)
    :param None focus: a tuple with the start and end positions of the 
       region to visualize
    :param False paint_tads: draw a box around the TADs defined for this
       experiment
    :param None axe: an axe object from matplotlib can be passed in order
       to customize the picture
    :param True show: either to pop-up matplotlib image or not
    :param True logarithm: show the logarithm values
    :param True normalized: show the normalized data (weights might have
       been calculated previously). *Note: white rows/columns may appear in
       the matrix displayed; these rows correspond to filtered rows (see*
       :func:`pytadbit.utils.hic_filtering.hic_filtering_for_modelling` *)*
    :param True relative: color scale is relative to the whole matrix of
       data, not only to the region displayed
    :param True decorate: draws color bar, title and axes labels
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param None clim: tuple with minimum and maximum value range for color
       scale. I.e. clim=(-4, 10)
    """
    if not isinstance(experiments, list) and not isinstance(experiments, tuple):
        experiments = [experiments]
        cols = 1
        rows = 1
    else:
        sqrtxpr = sqrt(len(experiments))
        cols = int(round(sqrtxpr + (0.0 if int(sqrtxpr)==sqrtxpr else .5)))
        rows = int(sqrtxpr+.5)
    notaxe = axe == None
    if not scale:
        scale = (8, 6)
    if notaxe and len(experiments) != 1:
        fig = plt.figure(figsize=(scale[0] * cols, scale[1] * rows))
    for i in xrange(rows):
        for j in xrange(cols):
            if i * cols + j >= len(experiments):
                break
            if notaxe and len(experiments) != 1:
                axe = fig.add_subplot(
                    rows, cols, i * cols + j + 1)
            if (isinstance(experiments[i * cols + j], tuple) or
                isinstance(experiments[i * cols + j], list)):
                if not axe:
                    fig = plt.figure(figsize=(scale[0] * cols, scale[1] * rows))
                    axe = fig.add_subplot(
                        rows, cols, i * cols + j + 1)                        
                xpr1 = experiments[i * cols + j][0]
                xpr2 = experiments[i * cols + j][1]
                img = xpr1.view(tad=tad, focus=focus, paint_tads=paint_tads,
                                axe=axe, show=False, logarithm=logarithm,
                                normalized=normalized, relative=relative,
                                decorate=decorate, savefig=False,
                                where='up', clim=clim)
                img = xpr2.view(tad=tad, focus=focus, paint_tads=paint_tads,
                                axe=axe, show=False, logarithm=logarithm,
                                normalized=normalized, relative=relative,
                                decorate=False, savefig=False, where='down',
                                clim=clim or img.get_clim())
                #axe = axe.twinx()
                #axe.set_aspect('equal',adjustable='box-forced',anchor='NE')
                if decorate:
                    plt.text(1.01, .5, 
                             'Experiment %s' % (xpr2.name),
                              rotation=-90, va='center', size='large',
                              ha='left', transform=axe.transAxes)
            else:
                xper = experiments[i * cols + j]
                if not xper.hic_data and not xper.norm:
                    continue
                xper.view(tad=tad, focus=focus, paint_tads=paint_tads,
                          axe=axe, show=False, logarithm=logarithm,
                          normalized=normalized, relative=relative,
                          decorate=decorate, savefig=False, clim=clim)
    if savefig:
        tadbit_savefig(savefig)
    if show:
        plt.show()

