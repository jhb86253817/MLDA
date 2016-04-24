# statistics for computing the real computational efficiency of vb

with open('en_9000.dat', 'r') as f:
    en = f.read()

en = en.split('\n')
en = en[1:-1]
en = [d.split() for d in en]
en_count = sum([len(d) for d in en])

en_set = [set(d) for d in en]
en_set_count = sum([len(d) for d in en_set])

print 'english, original count: %d\n' % en_count
print 'english, vocabulary count: %d\n' % en_set_count

with open('de_9000.dat', 'r') as f:
    de = f.read()

de = de.split('\n')
de = de[1:-1]
de = [d.split() for d in de]
de_count = sum([len(d) for d in de])

de_set = [set(d) for d in de]
de_set_count = sum([len(d) for d in de_set])

print 'german, original count: %d\n' % de_count
print 'german, vocabulary count: %d\n' % de_set_count
