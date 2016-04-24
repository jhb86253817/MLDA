# compute average length of each document

en = []
with open('en_50000.dat', 'r') as f:
    for l in f:
        if l != '50000\n':
            en.append(l.strip())
en = [s.split() for s in en]
en = [w for s in en for w in s]

zh = []
with open('zh_50000.dat', 'r') as f:
    for l in f:
        if l != '50000\n':
            zh.append(l.strip())
zh = [s.split() for s in zh]
zh = [w for s in zh for w in s]
print en[0]
print zh[0]
#print len(en)/50000.0
#print len(zh)/50000.0
