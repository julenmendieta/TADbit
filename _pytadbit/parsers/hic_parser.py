"""
November 7, 2013.

"""

from warnings import warn
from math import sqrt, isnan
from pytadbit.parsers.gzopen import gzopen
from pytadbit.interaction_matrix import InteractionMatrix

import re

HIC_DATA = True


def extract_headers(rownames, colnames, header_re=None):
    if not header_re:
        header_re = '|([^:]+):([0-9]+)-([0-9]+)'
    if isinstance(header_re, str):
        header_re = [header_re]
    for i in xrange(len(rownames)):
        rownames[i] = [re.findall(header_re[j], rownames[i][j])[0]
                       for j, r in enumerate(header_re)]
    for i in xrange(len(colnames)):
        colnames[i] = [re.findall(header_re[j], colnames[i][j])[0]
                       for j, r in enumerate(header_re)]        

# Exception to handle failed autoread.
class AutoReadFail(Exception):
    pass

# Helper functions for the autoreader.
def is_asymmetric(matrix):
    maxn = len(matrix)
    for i in range(maxn):
        maxi = matrix[i] # slightly more efficient
        for j in range(i+1, maxn):
            if maxi[j] != matrix[j][i]:
                if isnan(maxi[j]) and isnan(matrix[j][i]):
                    continue
                return True
    return False

def symmetrize(matrix):
    maxn = len(matrix)
    for i in range(maxn):
        for j in range(i+1, maxn):
            matrix[i][j] = matrix[j][i] = matrix[i][j] + matrix[j][i]


def autoreader(f, header_re=None):
    """
    Auto-detect matrix format of HiC data file.
    
    :param f: an iterable (typically an open file).
    
    :returns: A tuple with integer values and the dimension of
       the matrix.
    """

    # Skip initial comment lines and read in the whole file
    # as a list of lists.
    for line in f:
        if not line.startswith('#'):
            break
    items = [line.split()] + [line.split() for line in f]

    # Count the number of elements per line after the first.
    # Wrapping in a set is a trick to make sure that every line
    # has the same number of elements.
    S = set([len(line) for line in items[1:]])
    ncol = S.pop()
    # If the set 'S' is not empty, at least two lines have a
    # different number of items.
    if S:
        raise AutoReadFail('ERROR: unequal column number')

    nrow = len(items)
    # Auto-detect the format, there are only 4 cases.
    if ncol == nrow:
        print 'SAME NUMBER OF ROW COLUMNS'
        try:
            _ = [float(item) for item in items[0]]
            # Case 1: pure number matrix.
            print '  no header'
            header = False
            trim = 0
        except ValueError:
            # Case 2: matrix with row and column names.
            print '  header'
            header = True
            if header_re:
                trim = len(header_re)
            else:
                trim = 1
    else:
        print 'DIFFERENT NUMBER OF ROW COLUMNS', len(items[0]), len(items[1]), nrow, ncol
        if len(items[0]) == len(items[1]):
            # Case 3: matrix with row information.
            print '  no header'
            header = False
            if header_re:
                trim = len(header_re)
            else:
                trim = ncol - nrow
        else:
            # Case 4: matrix with header and row information.
            print '  header'
            header = True
            if header_re:
                trim = len(header_re)
            else:
                trim = ncol - nrow + 1
    # Keep header line
    if header and not trim:
        rownames = items.pop(0)
        nrow -= 1
    elif not trim:
        rownames = range(1, nrow + 1)
    else:
        del(items[0])
        nrow -= 1
        rownames = [tuple([a for a in line[:trim]]) for line in items]
    if trim:
        colnames = [c for c in items[0] if c]

    extract_headers(colnames, rownames, header_re)
    # print 'HEADER', header
    # Get the numeric values and remove extra columns
    what = int if HIC_DATA else float
    try:
        items = [[what(a) for a in line[trim:]] for line in items]
    except ValueError:
        if not HIC_DATA:
            raise AutoReadFail('ERROR: non numeric values')
        try:
            # Dekker data 2009, uses integer but puts a comma... 
            items = [[int(float(a)+.5) for a in line[trim:]] for line in items]
            warn('WARNING: non integer values')
        except ValueError:
            try:
                # Some data may contain 'NaN' or 'NA'
                items = [
                    [0 if a.lower() in ['na', 'nan']
                     else int(float(a)+.5) for a in line[trim:]]
                for line in items]
                warn('WARNING: NA or NaN founds, set to zero')
            except ValueError:
                raise AutoReadFail('ERROR: non numeric values')

    # Check that the matrix is square.
    ncol -= trim
    if ncol != nrow: raise AutoReadFail('ERROR: non square matrix')

    if is_asymmetric(items):
        warn('WARNING: input matrix not symmetric: symmetrizing')
        symmetrize(items)

    return tuple([a for line in items for a in line]), ncol, (rownames, colnames)


def read_matrix(things, parser=None, hic=True, **kwargs):
    """
    Read and checks a matrix from a file (using
    :func:`pytadbit.parser.hic_parser.autoreader`) or a list.

    :param things: might be either a file name, a file handler or a list of
        list (all with same length)
    :param None parser: a parser function that returns a tuple of lists
       representing the data matrix,
       with this file example.tsv:
       ::
       
         chrT_001	chrT_002	chrT_003	chrT_004
         chrT_001	629	164	88	105
         chrT_002	86	612	175	110
         chrT_003	159	216	437	105
         chrT_004	100	111	146	278

       the output of parser('example.tsv') might be:
       ``([629, 86, 159, 100, 164, 612, 216, 111, 88, 175, 437, 146, 105, 110,
       105, 278])``


    :param True hic: if False, TADbit assumes that files contains normalized
       data
    :returns: the corresponding matrix concatenated into a huge list, also
       returns number or rows

    """
    global HIC_DATA
    HIC_DATA = hic
    parser = parser or autoreader
    if not isinstance(things, list):
        things = [things]
    matrices = []
    for thing in things:
        if isinstance(thing, InteractionMatrix):
            matrices.append(thing)
        elif isinstance(thing, file):
            matrix, size, header = parser(thing, **kwargs)
            thing.close()
            matrices.append(InteractionMatrix([(i, matrix[i]) for i in xrange(size**2)
                                               if matrix[i]], size))
        elif isinstance(thing, str):
            try:
                matrix, size, header = parser(gzopen(thing, **kwargs))
            except IOError:
                if len(thing.split('\n')) > 1:
                    matrix, size, header = parser(thing.split('\n'), **kwargs)
                else:
                    raise IOError('\n   ERROR: file %s not found\n' % thing)
            matrices.append(InteractionMatrix([(i, matrix[i]) for i in xrange(size**2)
                                               if matrix[i]], size, sections=header))
        elif isinstance(thing, list):
            if all([len(thing)==len(l) for l in thing]):
                matrix  = reduce(lambda x, y: x+y, thing)
                size = len(thing)
            else:
                # print thing
                raise Exception('must be list of lists, all with same length.')
            matrices.append(InteractionMatrix([(i, matrix[i]) for i in xrange(size**2)
                                               if matrix[i]], size))
        elif isinstance(thing, tuple):
            # case we know what we are doing and passing directly list of tuples
            matrix = thing
            siz = sqrt(len(thing))
            if int(siz) != siz:
                raise AttributeError('ERROR: matrix should be square.\n')
            size = int(siz)
            matrices.append(InteractionMatrix([(i, matrix[i]) for i in xrange(size**2)
                                               if matrix[i]], size))
        elif 'matrix' in str(type(thing)):
            try:
                row, col = thing.shape
                if row != col:
                    raise Exception('matrix needs to be square.')
                matrix  = thing.reshape(-1).tolist()[0]
                size = row
            except Exception as exc:
                print 'Error found:', exc
            matrices.append(InteractionMatrix([(i, matrix[i]) for i in xrange(size**2)
                                               if matrix[i]], size))
        else:
            raise Exception('Unable to read this file or whatever it is :)')

    return matrices

