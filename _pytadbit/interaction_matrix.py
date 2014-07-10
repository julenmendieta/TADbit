"""
10 juil. 2014
"""

class InteractionMatrix(dict):
    """
    This may also hold the print/write-to-file matrix functions

    :param items: list of key, values
    :param size: size of the matrix
    :param None sections: list of descriptors to assign at each row/columns (two
       lists needed if not symmetric). One descriptors can be for example:
       ('Scer', 'chrIV', '50000') or ('chrX', '980000:1000000')
    :param False normalized: if matrix is already normalized
    :param True symmetric: if matrix is symmetric (the case for Hi-C data but
       not for 5-C)
    :param None scale: can be either a single number if all rows represent the
       same number of nucleotides, or a list if not.
    """
    def __init__(self, items, size, sections=None, normalized=False,
                 symmetric=True, scale=None):
        super(InteractionMatrix, self).__init__(items)
        self._size = size
        self._size2 = size**2
        self.normalized = normalized
        self.bias = None
        self.sections = sections
        self.symmetric = symmetric
        self.scale = scale if (isinstance(scale, list) or
                               isinstance(scale, tuple)) else Scale((scale, ))

    def __len__(self):
        return self._size

    def __getitem__(self, row_col):
        """
        slow one... for user
        for fast item getting, use self.get()
        """
        try:
            row, col = sorted(row_col) if self.symmetric else row_col
            pos = row * self._size + col
            if pos > self._size2:
                raise IndexError(
                    'ERROR: row or column larger than %s' % self._size)
            return self.get(pos, 0)
        except TypeError:
            if row_col > self._size2:
                raise IndexError(
                    'ERROR: position %d larger than %s^2' % (row_col,
                                                             self._size))
            return self.get(row_col, 0)

    def get_section(self, section):
        """
        Get a section of the current interaction matrix, matching with the
        current input word(s)

        :param section: word or list of words to search for in the
           InteractionMatrix

        :retruns: a new Interaction matrix
        """
        # TODO
        pass

    def get_as_tuple(self):
        """
        :returns: tuple corresponding to all values in the matrix
        """
        return tuple([self[i, j] for j  in xrange(len(self))
                      for i in xrange(len(self))])

    def get_matrix(self, focus=None, diagonal=True, normalized=False):
        """
        get the matrix
        """
        siz = len(self)
        if normalized and not self.bias:
            raise Exception('ERROR: experiment not normalized yet')
        if focus:
            start, end = focus
            start -= 1
        else:
            start = 0
            end   = siz
        if normalized:
            if diagonal:
                return [[self[i, j] / self.bias[i] / self.bias[j]
                         for i in xrange(start, end)]
                        for j in xrange(start, end)]
            else:
                mtrx = [[self[i, j] / self.bias[i] / self.bias[j]
                         for i in xrange(start, end)]
                        for j in xrange(start, end)]
                for i in xrange(start, end):
                    mtrx[i][i] = 1 if mtrx[i][i] else 0
                return mtrx
        else:
            if diagonal:
                return [[self[i, j] for i in xrange(start, end)]
                        for j in xrange(start, end)]
            else:
                mtrx = [[self[i, j] for i in xrange(start, end)]
                        for j in xrange(start, end)]
                for i in xrange(start, end):
                    mtrx[i][i] = 1 if mtrx[i][i] else 0
                return mtrx


class Scale(tuple):
    def __init__(self, items):
        super(Scale, self).__init__(items)
        self.__scale = items[0]
    def __getitem__(self, num):
        return self.__scale
