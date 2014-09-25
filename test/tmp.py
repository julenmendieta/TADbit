
from pytadbit                        import Chromosome, Experiment, load_chromosome
from pytadbit                        import tadbit, batch_tadbit, find_tads
from pytadbit                        import align_experiments
from pytadbit.tad_clustering.tad_cmo import optimal_cmo
from pytadbit.imp.structuralmodels   import load_structuralmodels
from pytadbit.imp.impmodel           import load_impmodel_from_cmm
from pytadbit.eqv_rms_drms           import rmsdRMSD_wrapper
from pytadbit.utils.normalize_hic    import iterative
from os                              import system, path, chdir
from warnings                        import warn
from distutils.spawn                 import find_executable
from pytadbit.visual.heatmap         import interaction_map
from pytadbit import find_tads
from pytadbit.visual.tad_plot         import tad_density_plot

PATH = '/Users/fransua/Box/tadbits/tadbit-1.0/test/'
experiments = {}

experiments['exp1'] = Experiment('exp1', resolution=40000,
                                 hic_data=PATH + '/40Kb/chrT/chrT_A.tsv')
experiments['exp2'] = Experiment('exp2', resolution=20000,
                                 hic_data=PATH + '/40Kb/chrT/chrT_B.tsv')
experiments['exp3'] = Experiment('exp3', resolution=20000,
                                 hic_data=PATH + '/40Kb/chrT/chrT_C.tsv')
experiments['exp4'] = Experiment('exp4', resolution=20000,
                                 hic_data=PATH + '/40Kb/chrT/chrT_D.tsv')

interaction_map([(experiments['exp1'], experiments['exp2']), experiments['exp1'], experiments['exp2']])


find_tads(experiments['exp1'])

tad_density_plot(experiments['exp1'])

Experiment('exp4', resolution=20000, hic_data=PATH + '/40Kb/chrT/chrT_D.tsv')

from pytadbit                        import Chromosome, Experiment, load_chromosome
exp = Experiment('hola', hic_data='/Users/fransua/Downloads/visualizingthesematricesintadbit/6711_HFSP-MCF7-E3h-R2-AA.matrix')










