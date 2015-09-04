import WTaggingPerf as ww
import yaml
import os
import operator
import sys

folder = sys.argv[1]

def get_score(fname,directory=''):
    f = open(directory+fname, 'r')

    schema = yaml.load(f)

    taggers = ww.generate_taggers(schema, fname.replace('.yaml',''))

    rejection_power = ww.plot_roc(taggers, schema, 'ROC_plots_full/ROC_plots_' + fname.replace('.yaml','') + '.pdf', save_arr = True, logscale=False)
    return rejection_power

# get all config files to test
files = [f for f in os.listdir(folder+'/') if os.path.isfile(folder+'/'+f) and f.startswith('config-trained-')]
#files = ['config-trained-33.yaml']
# get the scores for each config
scores = {}
for f in files:
    print f
    scores[f] = get_score(f,folder+'/')
    #raw_input()

# sort on scores
sorted_scores = sorted(scores.items(), key = operator.itemgetter(1) )

print sorted_scores

#save to file
f_out = open('best_rejection_full.txt','w')
for s in range(len(sorted_scores)):
    f_out.write(sorted_scores[s][0] + ': ' + str(sorted_scores[s][1])+'\n')

f_out.close()
