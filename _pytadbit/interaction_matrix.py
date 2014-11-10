"""
10 juil. 2014
"""

from cPickle import load, dump
from gzip import open as gopen

class InteractionMatrix(dict):
    """
    TODO: This may also hold the print/write-to-file matrix functions

    :param items: list of key, values
    :param size: size of the matrix
    :param None sections: list of descriptors to assign at each row/columns.
       One descriptors can be for example:
       ('Scer', 'chrIV', '50000') or ('chrX', '980000:1000000')
    :param None section_sizes: a dictionary containing the size of each section,
       subsection. I.e.: {'Scer': 1200000, ('Scer', 'chrIV'): 150000}
    :param False normalized: if matrix is already normalized
    :param None scale: can be either a single number if all rows represent the
       same number of nucleotides, or a list if not.
    """
    def __init__(self, items, size, name=None, sections=None, normalized=False,
                 normalization=None, scale=1, section_sizes=None):
        super(InteractionMatrix, self).__init__(items)
        self.name = name or 'NoName'
        self._size = size
        self._size2 = size**2
        self.normalized = normalized
        self._normalization = normalization or 'None'
        self.bias = None
        if not sections:
            # create sections...
            self.sections = [(item, ) for item in items]
        else:
            self.sections = sections

        # calculate the size of each section
        self.section_sizes = {}
        if section_sizes:
            self.section_sizes = section_sizes
        else:
            self.__size_sections__()

        self.scale = scale if (isinstance(scale, list) or
                               isinstance(scale, tuple)) else Scale((scale, ))
        self._real_size = sum([self.scale[i] for i in xrange(self._size)])

    def __size_sections__(self):
        for section in self.sections:
            key = tuple((section[0],))
            self.section_sizes.setdefault(key, 0)
            self.section_sizes[key] += 1
            for i in section[1:-1]:
                key += tuple((i,))
                self.section_sizes.setdefault(key, 0)
                self.section_sizes[key] += 1

    def __len__(self):
        return self._size

    def __getitem__(self, row_col):
        """
        slow one... for user
        for fast item getting, use self.get()
        """
        try:
            row, col = sorted(row_col)
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

    def write(self, fname, format='pik', gzip=True, headers=True, focus=None):
        """
        :param fname: save matrix to file
        :param pik format: can be either
            - pik (python pickle object)
            - mtx (matrix, largest file size)
            - abc (3 column matrix)
        :param True header: only used when writting matrices
        :param None focus: a tuple with (start, end) positions to print. These
           positions are not genomic coordinates but bin indexes.
        """
        thisopen = gopen if gzip else open
        out = thisopen(fname, 'w' + ('b' if gzip else ''))
        if focus:
            if focus[0] < 0 or focus[0] > self._size:
                raise IndexError('Focus out of bond')
            thisrange = range(focus[0], focus[1])
        else:
            thisrange = range(self._size)
        if format == 'mtx':
            out.write('\t' +
                      '\t'.join([('_'.join(self.sections[i]))
                                 for i in thisrange]) + '\n')
            for row in thisrange:
                if headers:
                    rowh = '_'.join(self.sections[row])
                    out.write(rowh + '\t' + '\t'.join(
                        [str(self[row * self._size + col])
                         for col in thisrange]) + '\n')
                else:
                    out.write('\t'.join(
                        [str(self[row * self._size + col])
                         for col in thisrange]) + '\n')
        elif format == 'pik':
            if focus:
                dump(self.get_sample(focus), out)
            else:
                dump(self, out)
        elif format == 'abc':
            for row in thisrange:
                rowh = '_'.join(self.sections[row])
                for col in thisrange:
                    colh = '_'.join(self.sections[col])
                    out.write('%s\t%s\t%s\n' % ('_'.join(rowh), '_'.join(colh),
                                                self[row * self._size + col]))
        else:
            out.close()
            raise NotImplementedError('Format should be one of pik, mtx or abc')
        out.close()


    def get_scaled(self, scale):
        """
        Lower the resolution, sections will be unchanged but the last, that will
        be concatenated.
        WARNING: this will produce binned data.
        
        :param resolution: should be multiple of actual resolution
        """
        vals = {}
        maxlen = len(max(self.section_sizes.keys(), key=len))
        sec2scale = dict([(sec, sum([self.scale[i] for i in xrange(self._size)
                                     if sec == self.sections[i][:-1]])
                           / self.section_sizes[sec])
                          for sec in self.section_sizes])
        newsize = sum([self.section_sizes[sec] * sec2scale[sec] / scale
                       for sec in self.section_sizes if len(sec) == maxlen])
        newsections = [sec for sec in self.section_sizes
                       for i in xrange(self._size)
                       if sec == self.sections[i][:-1]]
        print newsections
        # j = 0
        # for supsec in newsections:
        #     for i in xrange(self.section_sizes[supsec]):
        #         newsections[i+j] = supsec
        #     j += i
        for row in xrange(self._size):
            newrow = row  * self.scale[row] / scale
            for col in xrange(self._size):
                newcol = col * self.scale[col] / scale
                cell = (newcol) + (newrow) * newsize
                vals.setdefault(cell, 0)
                vals[cell] += self[row * self._size + col]

        return InteractionMatrix(
            vals.items(), newsize,
            normalized=self.normalized,
            normalization=self._normalization + ' inherited',
            sections = newsections,
            name = self.name + ('_scaled:%s' % (scale)))
       

    def get_sample(self, focus):
        """
        pick portion of the matrix
        :param focus: a tuple with (start, end) positions to print. These
           positions are not genomic coordinates but bin indexes.
        """
        idxs = range(focus[0], focus[1])
        
        return InteractionMatrix(
            [(i, self[j]) for i, j in enumerate(idxs)], focus[1] - focus[0],
            normalized=self.normalized,
            normalization=self._normalization + ' inherited',
            sections = [self.sections[i] for i in idxs],
            name = self.name + ('%s-%s' % focus))


    def get_section(self, section):
        """
        Get a section of the current interaction matrix, matching with the
        current input word(s)

        :param section: word or list of words to search for in the
           InteractionMatrix

        :retruns: a new Interaction matrix
        """
        if isinstance(section, list):
            section_set = set(tuple(section))
        elif isinstance(section, str):
            section_set = set((section,))
        elif isinstance(section, tuple):
            section_set = set(section)
        else:
            raise NotImplementedError('should pass either list tuple or str')
        idxs = [row for row, names in enumerate(self.sections)
                if not section_set.difference(names)]
        if not idxs:
            raise IndexError('Section(s) not found')
        if 1 + idxs[-1] - idxs[0] != len(idxs):
            raise Exception('ERROR: section name should correspond to uniq ' +
                            'section of the matrix')
        new_im = InteractionMatrix(
            [], self.section_sizes.get(section, len(idxs)),
            normalized=self.normalized,
            section_sizes=dict([(tuple(section_set.difference(k)),
                                 self.section_sizes[k])
                                for k in self.section_sizes
                                if section_set.difference(k)]),
            normalization=self._normalization + ' inherited', sections = [],
            name = '_'.join(section))
        minidx = idxs[0]
        if isinstance(section, str):
            pos = self.sections[minidx].index(section)
        else:
            pos = len(section)
        for k in xrange(len(idxs)):
            new_im[idxs[k] - minidx] = self[idxs[k]]
            new_im.sections.append(self.sections[idxs[k]][:pos] +
                                   self.sections[idxs[k]][pos + 1:])
        return new_im

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
