# statistics for computing the real computational efficiency of vb

with open('en_50000.dat', 'r') as f:
    en = f.read()

en = en.split('\n')
en = en[1:-1]
en = [d.split() for d in en]
en_count = sum([len(d) for d in en])

en_set = [set(d) for d in en]
en_set_count = sum([len(d) for d in en_set])

print 'english, original count: %d\n' % en_count
print 'english, vocabulary count: %d\n' % en_set_count

with open('zh_50000.dat', 'r') as f:
    zh = f.read()

zh = zh.split('\n')
zh = zh[1:-1]
zh = [d.split() for d in zh]
zh_count = sum([len(d) for d in zh])

zh_set = [set(d) for d in zh]
zh_set_count = sum([len(d) for d in zh_set])

print 'chinese, original count: %d\n' % zh_count
print 'chinese, vocabulary count: %d\n' % zh_set_count
