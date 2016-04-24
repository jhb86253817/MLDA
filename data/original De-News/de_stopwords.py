# -*- coding: utf-8 -*-

# To map special character to its alternative

import json

de_stopwords = []
with open('stop-words_german_2_de.txt', 'r') as f:
    for l in f:
        de_stopwords.append(l.strip())
de_stopwords = [w.replace('ß', 'ss') for w in de_stopwords]
de_stopwords = [w.replace('ü', 'ue') for w in de_stopwords]
de_stopwords = [w.replace('ä', 'ae') for w in de_stopwords]
de_stopwords = [w.replace('ö', 'oe') for w in de_stopwords]

with open('de_stopwords.json', 'w') as f:
    f.write(json.dumps(de_stopwords))


