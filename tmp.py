
from pytadbit.parsers.hic_parser import read_matrix

a = read_matrix('/home/fransua/Box/tadbits/tadbit-1.0/chrXX.tsv')[0]

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

