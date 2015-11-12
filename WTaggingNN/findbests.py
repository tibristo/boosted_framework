import WTaggingPerf as ww
import yaml
import os
import operator
import sys
import argparse
 
def get_score(fname,directory=''):
    f = open(directory+fname, 'r')

    schema = yaml.load(f)

    taggers = ww.generate_taggers(schema, fname.replace('.yaml',''))

    rejection_power = ww.plot_roc(taggers, schema, 'ROC_plots_full/ROC_plots_' + fname.replace('.yaml','') + '.pdf', save_arr = True, logscale=False)
    return rejection_power




def main(args):
    parser = argparse.ArgumentParser('Test the neural network on some new datasets.')
    parser.add_argument('--input-folder', default = 'trained_full/', dest='folder', help = 'Folder containing the config files for the datasets that are being tested. Default is trained_full/')
    parser.add_argument('--output-id', dest = 'output_fname', help = 'A string that is included in the output log file.  The file name will be of the form `best_rejection_[output-id].txt`. Default is mc15_v0.1.')
    parser.add_argument('--file-id', dest='file_id', help = 'Partial string to match when looking for the config files in the folder given. This is used to filter out old config files that we do not want to use. Default is mc15_v0.1.')
    parser.set_defaults(file_id='mc15_v0.1')
    parser.set_defaults(output_fname='mc15_v0.1')
    args = parser.parse_args()

    # get all config files to test
    files = [f for f in os.listdir(args.folder+'/') if os.path.isfile(args.folder+'/'+f) and f.find(args.file_id)!= -1]#f.startswith('config-trained-paramID')]

    # get the scores for each config
    scores = {}

    # remember that all cross validation folds have to be averaged to get the proper score!
    for f in files:
        print f
        scores[f] = get_score(f,args.folder+'/')

    # sort on scores
    sorted_scores = sorted(scores.items(), key = operator.itemgetter(1) )

    print sorted_scores

    #save to file
    f_out = open('best_rejection_'+args.output_fname+'.txt','w')
    for s in range(len(sorted_scores)):
        f_out.write(sorted_scores[s][0] + ': ' + str(sorted_scores[s][1])+'\n')

    f_out.close()
    
if __name__=='__main__':
    main(sys.args)
