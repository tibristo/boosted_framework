class scoresObj:
    #import math
    def __init__(self, uepochs=10, log_learn=-10, momentum=0.5, sepochs=20, log_reg=-10, auc=1):
        import math
        #'uepochs': array([20]), u'log_learning': array([-10.]), u'momentum': array([ 0.7]), u'sepochs': array([40]), u'log_regularize': array([-10.])
        self.uepochs = uepochs
        self.log_learn = log_learn
        self.learn = math.exp(log_learn)
        self.momentum = momentum
        self.sepochs = sepochs
        self.log_regularize = log_reg
        self.regularize = math.exp(log_reg)
        self.auc = auc

    def plotuepochs(self):
        return self.uepochs,self.auc

    def plotsepochs(self):
        return self.sepochs,self.auc

    def plotmomentum(self):
        return self.momentum,self.auc

    def plotlearning(self):
        return self.learn,self.auc

    def plotregularize(self):
        return self.regularize,self.auc

import re

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

import os
import numpy as np
files = [f for f in os.listdir('.') if f.endswith('small') == True]
sort_nicely(files)

print files
scores = np.zeros(len(files))
scores_names = []
scores_list = []

for i, f in enumerate(files):
    f_in = open(f,'r')
    params = {'uepochs':0, 'log_learning': 0, 'momentum': 0, 'sepochs':0, 'log_regularize':0}
    params_l = ''
    for l in f_in:
        if l.startswith('{'):
            params_l = l.strip()
    last = l.strip()
    params_spl = params_l.split(',')
    for p in params_spl:
        name = p.split('\'')[1::2][0]
        val = p[p.find('[')+1:p.find(']')]
        params[name] = float(val)
    scores[i] = float(last)        
    scores_list.append(scoresObj(params['uepochs'],params['log_learning'],params['momentum'],params['sepochs'],params['log_regularize'], scores[i]))

    scores_names.append(f)
    f_in.close()

max_idx = np.argmax(scores)
print 'max: ' + str(scores[max_idx])
print 'file name: ' + str(scores_names[max_idx])

max_file = open(scores_names[max_idx], 'r')
for l in max_file:
    if l.startswith('--formula') or l.startswith('{'):
        print l.strip()
max_file.close()


import matplotlib.pyplot as plt

sepochs = np.zeros(len(scores_list))
uepochs = np.zeros(len(scores_list))
momentum = np.zeros(len(scores_list))
learn = np.zeros(len(scores_list))
regul = np.zeros(len(scores_list))
auc = np.zeros(len(scores_list))
for i, s in enumerate(scores_list):
    sepochs[i] = s.plotsepochs()[0]
    uepochs[i] = s.plotuepochs()[0]
    momentum[i] = s.plotmomentum()[0]
    regul[i] = s.plotregularize()[0]
    learn[i] = s.plotlearning()[0]
    auc[i]= s.auc

#plt.scatter(np.arange(1, len(scores_list)+1, dtype=int), auc)
plt.scatter(sepochs, auc)
#plt.scatter(uepochs, auc)
#plt.scatter(momentum, auc)
#plt.scatter(regul, auc)
#plt.scatter(learn, auc)
plt.show()
