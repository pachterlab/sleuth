import itertools
for x in itertools.product('ACGT',repeat=6):
    print "'{0}'".format( ''.join(x) )
