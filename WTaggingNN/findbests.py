import WTaggingPerf as ww
import yaml
import os
import operator
import sys

folder = sys.argv[1]
file_id = 'jz5_bkg_v2'

def get_score(fname,directory=''):
    f = open(directory+fname, 'r')

    schema = yaml.load(f)

    taggers = ww.generate_taggers(schema, fname.replace('.yaml',''))

    rejection_power = ww.plot_roc(taggers, schema, 'ROC_plots_full/ROC_plots_' + fname.replace('.yaml','') + '.pdf', save_arr = True, logscale=False)
    return rejection_power

# get all config files to test
files = [f for f in os.listdir(folder+'/') if os.path.isfile(folder+'/'+f) and f.find(file_id)!= -1]#f.startswith('config-trained-paramID')]

# get the scores for each config
scores = {}

# remember that all cross validation folds have to be averaged to get the proper score!
for f in files:
    print f
    scores[f] = get_score(f,folder+'/')
    #raw_input()

# sort on scores
sorted_scores = sorted(scores.items(), key = operator.itemgetter(1) )

print sorted_scores

#save to file
f_out = open('best_rejection_jz5_v2_full.txt','w')
#f_out = open('best_rejection_cv_full_bdt_vars.txt','w')
for s in range(len(sorted_scores)):
    f_out.write(sorted_scores[s][0] + ': ' + str(sorted_scores[s][1])+'\n')

f_out.close()
