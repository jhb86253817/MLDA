# -*- coding: utf-8 -*-

import re
import json
import sys
import os
from nltk.corpus import stopwords
from collections import Counter

# number of documents for training
pair_num_trn = 9000
# number of documents for inference
pair_num_inf = 667

def readtxt(filename):
    txt = []
    with open(filename, 'r') as f:
        pointer = -1
        for l in f:
            if l.startswith('<DOC de-news'):
                pointer += 1
                txt.append('')
            else:
                txt[pointer] += l
    return txt

en = readtxt('en.txt')
de = readtxt('de.txt')
print len(en)
# remove newline character
en = [s.replace('\n', ' ') for s in en]
de = [s.replace('\n', ' ') for s in de]
# lowercase all the strings
en = [s.lower() for s in en]
de = [s.lower() for s in de]
# split each document to tokens
en = [s.split() for s in en]
de = [s.split() for s in de]
# remove non-alphabet character
en = [[w for w in s if w.isalpha()] for s in en]
de = [[w for w in s if w.isalpha()] for s in de]
# load stopwords for English and German
en_stopwords = set(stopwords.words('english'))
with open('de_stopwords.json', 'r') as f:
    de_stopwords = json.loads(f.read())
de_stopwords = set(de_stopwords)
# remove stopwords
en = [[w for w in s if not w in en_stopwords] for s in en]
de = [[w for w in s if not w in de_stopwords] for s in de]
# convert token list to string
en = [' '.join(s) for s in en]
de = [' '.join(s) for s in de]
# split to training and inference files
en_trn = en[:pair_num_trn]
en_inf = en[pair_num_trn:]
de_trn = de[:pair_num_trn]
de_inf = de[pair_num_trn:]
# save files
# training files for English
with open('en_'+str(pair_num_trn)+'.dat', 'w') as f:
    f.write(str(pair_num_trn)+'\n')
    for doc in en_trn:
        f.write(doc.encode('utf-8'))
        f.write('\n')
# inference files for English
with open('en_'+str(pair_num_inf)+'_inf.dat', 'w') as f:
    f.write(str(pair_num_inf)+'\n')
    for doc in en_inf:
        f.write(doc.encode('utf-8'))
        f.write('\n')
# training files for German
with open('de_'+str(pair_num_trn)+'.dat', 'w') as f:
    f.write(str(pair_num_trn)+'\n')
    for doc in de_trn:
        f.write(doc.encode('utf-8'))
        f.write('\n')
# inference files for German
with open('de_'+str(pair_num_inf)+'_inf.dat', 'w') as f:
    f.write(str(pair_num_inf)+'\n')
    for doc in de_inf:
        f.write(doc.encode('utf-8'))
        f.write('\n')
