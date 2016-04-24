# -*- coding: utf-8 -*-
from __future__ import division
import re
import json
import sys
import os
from nltk.corpus import stopwords
from collections import Counter
import matplotlib.pyplot as plt
from numpy.random import normal

#if len(sys.argv) != 2:
#    sys.exit('Please specify the number of docs for training')
#
## number of doc pair need to extract for training
#pair_num_trn = int(sys.argv[1])
#
#print "Generating %d doc pairs for training\n" % pair_num_trn
#
########################################################################################
#print '#############################################'
#print 'Extracting contents from XML file'
#print '#############################################'
#en = []
#zh = []
## extracting content from xml file
#with open('wikicomp-2014_enzh.xml', 'r') as f:
#    doc_num = 0
#    # indicate whether current line a part of content or not
#    content = False
#    # indicate which language the current content is
#    en_current = True
#    # record previous two lines
#    line_pre = ""
#    line_prepre = ""
#    # have to read 1 line a time because the file is large
#    for line in f:
#        # if current line a part of content 
#        if content:
#            # check if the end of the content
#            found = re.search(r'</content>', line)
#            if found:
#                # if enough doc found, quit
#                if doc_num == pair_num_trn*2:
#                    break
#                content = False
#                en_current = not en_current
#                continue
#            # add current line to the list
#            if en_current == True:
#                en[-1] += line
#            else:
#                zh[-1] += line
#        else:
#            found = re.search(r'<content>', line)
#            # if find the start of content
#            # change state of variables
#            if found:
#                doc_num += 1
#                content = True
#                # also add article name to content
#                article_name = re.search(r'name="(.+)"', line_prepre)
#                if en_current == True: 
#                    en.append('')
#                    en[-1] += article_name.group(1)+' '
#                else:
#                    zh.append('')
#                    zh[-1] += article_name.group(1)+' '
#        line_prepre = line_pre
#        line_pre = line
########################################################################################
## remove xml elements from extracted content
#print '#############################################'
#print 'Removing XML elements from extracted contents'
#print '#############################################'
## remove link 
#en = [re.sub(r'<link.+?>(.*?)</link>', r'\1', s) for s in en]
#zh = [re.sub(r'<link.+?>(.*?)</link>', r'\1', s) for s in zh]
## remove table
#en = [re.sub(r'<table>[\w\W]+?</table>', '', s) for s in en]
#zh = [re.sub(r'<table>[\w\W]+?</table>', '', s) for s in zh]
## remove math
#en = [re.sub(r'<math>[\w\W]+?</math>', '', s) for s in en]
#zh = [re.sub(r'<math>[\w\W]+?</math>', '', s) for s in zh]
## remove &quot
#en = [re.sub(r'&?quot;', ' ', s) for s in en]
#zh = [re.sub(r'&?quot;', ' ', s) for s in zh]
## remove &apos
#en = [re.sub(r'&?apos;', ' ', s) for s in en]
#zh = [re.sub(r'&?apos;', ' ', s) for s in zh]
## remove &amp
#en = [re.sub(r'&?amp;', ' ', s) for s in en]
#zh = [re.sub(r'&?amp;', ' ', s) for s in zh]
## remove reference and the thing after it
#en = [re.sub(r'<p><h>References[\w\W]*</p>$', '', s) for s in en]
#zh = [re.sub(r'<p><h>參考文獻[\w\W]*</p>$', '', s) for s in zh]
#zh = [re.sub(r'<p><h>参考文献[\w\W]*</p>$', '', s) for s in zh]
#zh = [re.sub(r'<p><h>參考資料[\w\W]*</p>$', '', s) for s in zh]
#zh = [re.sub(r'<p><h>資料來源[\w\W]*</p>$', '', s) for s in zh]
#zh = [re.sub(r'<p><h>外部链接[\w\W]*</p>$', '', s) for s in zh]
#zh = [re.sub(r'<p><h>参见[\w\W]*</p>$', '', s) for s in zh]
#zh = [re.sub(r'<p><h>参考[\w\W]*</p>$', '', s) for s in zh]
## remove <p>
#en = [re.sub(r'<p>', '', s) for s in en]
#zh = [re.sub(r'<p>', '', s) for s in zh]
## remove </p>
#en = [re.sub(r'</p>', '', s) for s in en]
#zh = [re.sub(r'</p>', '', s) for s in zh]
## remove <h>
#en = [re.sub(r'<h>', '', s) for s in en]
#zh = [re.sub(r'<h>', '', s) for s in zh]
## remove </h>
#en = [re.sub(r'</h>', '', s) for s in en]
#zh = [re.sub(r'</h>', '', s) for s in zh]
## remove http address
#en = [re.sub(r'http://.*', '', s) for s in en]
#zh = [re.sub(r'http://.*', '', s) for s in zh]

#prop = [len(zh[i])/len(en[i]) for i in range(pair_num_trn)]
#print len(prop)
#with open('doc_len_stats.json', 'w') as f:
#    f.write(json.dumps(prop))

with open('doc_len_stats.json', 'r') as f:
    prop = json.loads(f.read())

#plt.figure(1)
#plt.hist(prop, bins=100, range=[0, 10])
#plt.show()

#for i in range(100):
#    print i, ': ', prop[i]

doc_ids = [i for i in range(len(prop)) if prop[i]>0.5 and prop[i]<2]
doc_ids = doc_ids[:55000]
print len(doc_ids)
with open('doc_ids.json', 'w') as f:
    f.write(json.dumps(doc_ids))
