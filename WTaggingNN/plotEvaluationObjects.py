import sys
import modelEvaluation as ev
import os
import pickle
import numpy as np
import pandas as pd
from sklearn.externals import joblib
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score
import operator
import re
import matplotlib.pyplot as plt
from math import floor, log10
job_id = 'features_l_2_10'


# set up the job ids
key = 'config-output-paramID'
#key = 'features_l_2_10ID'
#full_dataset = 'persist/data_features_nc_2_10_v2_100.pkl'


print 'finished creating new full objects'
#sys.exit()
#raw_input()
files = [f for f in os.listdir('evaluationObjects/') if f.find(key)!=-1 and f.endswith('.pickle')]


def createDataframe(key, files):
    '''
    create a dataframe containing the info for all of the evaluation objects contained in the files list.
    '''
    data = []
    for f in files:
        print 'plotting file: ' + f
        with open('evaluationObjects/'+f,'r') as p:
            model = pickle.load(p)
        # the filename contains the iteration number and the cv number which together form a unique id.
        #  we have a string like paramID_iternum_key_cvnum[_full].pickle
        # first number in the filename is the iternum
        numbers = re.findall(r'\d+', f)
        iter_num = numbers[0]
        # last number is the cv number, unless we have _full_v2.pickle
        cv_num = numbers[-1]
        # now we want to get max_depth parameter
        params = model.params
        print params
        params_tofind = {'learning':-999.0,'momentum':-999.0, 'regularize':-999.0,'uepochs':-999.0,'sepochs':-999.0}
        for k in params_tofind.keys():
            vals = [value for key, value in params.items() if k in key.lower()]
            if vals[0] is not None:
                params_tofind[k] = float(vals[0]) if k.find('epochs') == -1 else int(vals[0])
                
        # get the f1 scores
        f1_train = model.train_f1
        f1_test = model.test_f1
        # get the recall scores
        recall_train = model.train_recall
        recall_test = model.test_recall
        # get the accuracy scores
        accuracy_train = model.train_accuracy
        accuracy_test = model.test_accuracy
        # get the precision scores
        precision_train = model.train_precision
        precision_test = model.test_precision
        # get the bkg rejection power scores
        bkg_rej_train = model.ROC_rej_power_05_train
        bkg_rej_test = model.ROC_rej_power_05
        thedict  = {'test_id':int(iter_num), 'cv_id': int(cv_num),'momentum':params_tofind['momentum'],'regularize':params_tofind['regularize'], 'learning':params_tofind['learning'], 'uepochs':int(params_tofind['uepochs']), 'sepochs':int(params_tofind['sepochs']),'f1_train':f1_train,'f1_test':f1_test,'recall_train':recall_train,'recall_test':recall_test,'precision_train':precision_train,'precision_test':precision_test,'accuracy_train':accuracy_train,'accuracy_test':accuracy_test,'bkg_rej_train':bkg_rej_train,'bkg_rej_test':bkg_rej_test,'filename': f}
        data.append(thedict)
        
    return data


def plotValidationCurve(key, parameters, train_scores_mean, test_scores_mean, train_scores_std, test_scores_std, param_range):
    '''
    Plot a validation curve. For the moment this does a validation curve for the bdt. A bdt normally depends on the tree depth
    for over and under-fitting.

    Keyword args:
    key -- file identifier
    parameters --- the values for the training parameters: dictionary with learning_rate: lr, n_est: n, max_depth: md
    mean and std for training and testing scores
    param_range --- the range of the max_depth variable
    '''

    # set up the file name
    fname = "validation_curves/"+key+'_lr_'+str(parameters['learning'])+'_n_'+str(parameters['n_estimators'])+'.png'
    plt.clf()
    plt.title("Validation Curve with BDT")
    plt.xlabel("$Max depth$")
    plt.ylabel("Accuracy")
    # find out what the y limits should be
    y_up = max(np.amax(train_scores_mean), np.amax(test_scores_mean))
    y_do = min(np.amin(train_scores_mean), np.amin(test_scores_mean))
    plt.ylim(y_do*0.95, y_up*1.05)
    plt.plot(param_range, train_scores_mean, label="Training score", color="r")
    plt.fill_between(param_range, train_scores_mean - train_scores_std,
                                      train_scores_mean + train_scores_std, alpha=0.2, color="r")
    plt.plot(param_range, test_scores_mean, label="Cross-validation score",
                              color="g")
    plt.fill_between(param_range, test_scores_mean - test_scores_std,
                                      test_scores_mean + test_scores_std, alpha=0.2, color="g")
    plt.legend(loc="best")
    plt.savefig(fname)


recreate_csv = True
columns = ['test_id','cv_id','momentum','regularize','learning','uepochs','sepochs','f1_train','f1_test','recall_train','recall_test','precision_train','precision_test','accuracy_train','accuracy_test','bkg_rej_train','bkg_rej_test','filename']


if recreate_csv:
    #get a dict of the data
    data_dict = createDataframe(key, files)
    # now create a dataframe out of the data
    df = pd.DataFrame(data_dict)
    # save it as a csv so that we don't have to do this all over again every time.
    #if save_csv:
    csv_file = open('data_'+key+'.csv','w')
    # write the column names to the csv file
    csv_file.write('test_id,cv_id,momentum,regularize,learning,uepochs,sepochs,f1_train,f1_test,recall_train,recall_test,precision_train,precision_test,accuracy_train,accuracy_test,bkg_rej_train,bkg_rej_test,filename\n')
    for i,d in enumerate(data_dict):
        # write the keys to the file if it is the first one
        num_keys = len(d.keys())
        #if i != 0:
        for j,k in enumerate(columns):
            csv_file.write(str(d[k]))
            if j != num_keys-1:
                csv_file.write(',')
            else:
                csv_file.write('\n')
    csv_file.close()

else:
    # we can read it in from the csv file
    # check that the file exists
    filefound = False
    fname = 'data_'+key+'.csv'
    while not filefound:
        if not os.path.isfile(fname):
            print "file "+fname+" doesn't exist! Enter another name to try"
            fname = raw_input()
            continue
        filefound = True
        
    df = pd.read_csv(fname)

# now we want to be able to the dataframe sort on a number of criteria and create meaningful stats
# easiest to run in interactive mode first.
grouped = df.groupby(['test_id']) # this combines all of the cv folds
# get the mean scores
gmean = grouped.mean()
gstd = grouped.std()

# sort the groups based on bkg rej power at 50%
srt_bkg = gmean.sort('bkg_rej_test',ascending=False)

# before writing to csv, format it to only use 3 sig digits
to_round = ['momentum','regularize','learning','accuracy_test','bkg_rej_test']
for r in to_round:
    srt_bkg[r] = srt_bkg[r].map(lambda x: round(x, -int(floor(log10(x)))+2) )
# make sure the epochs are ints
srt_bkg[['uepochs','sepochs']] = srt_bkg[['uepochs','sepochs']].astype(int)
    
srt_bkg.to_csv('data_'+key+'_sortedresults.txt',columns=['momentum','regularize','learning','uepochs','sepochs','accuracy_test','bkg_rej_test'])#,float_format='%2.2f')

# now write the top ones to a file
#f = open('data_'+key+'_sortedresults.txt','w')


#f.write()

#f.close()


# get the index
idx = gmean.index
'''
# get the levels - this stores the keys for the different groupby objects
lvls = idx.levels
# lrate is lvls[0], n_est is 1, max_depth is 2
# lvls is a FrozenList with each element being a Frozen64Index
lrates = lvls[0].values
n_est = lvls[1].values
md = lvls[2].values
training_points = np.zeros(len(md))
training_std = np.zeros(len(md))
testing_points = np.zeros(len(md))
testing_std = np.zeros(len(md))
# now we can get the scores for the different max_depths

for lr in lrates:
    for n in n_est:
        # store the values for all possible depths. this is what we want to plot!
        for i,d in enumerate(md):
            #print gmean.loc[lr,n,d]['accuracy_test']
            training_points[i] = gmean.loc[lr,n,d]['accuracy_test']
            training_std[i] = gstd.loc[lr,n,d]['accuracy_test']
            testing_points[i] = gmean.loc[lr,n,d]['accuracy_train']
            testing_std[i] = gstd.loc[lr,n,d]['accuracy_train']
        # now we can plot these
        #print training_points
        plotValidationCurve(key, {'learning_rate':lr, 'n_estimators':n}, training_points, testing_points, training_std, testing_std, md)
     
# we also want to know how all of the max_depth values did, for all parameters.
grouped_md = df.groupby(['max_depth','test_id'])
md_mean = grouped_md.mean()
for m in md:
    md_m = md_mean.loc[m]
    plt.scatter(np.repeat([m],len(md_m)),md_m['accuracy_test'], label=str(m))
    plt.scatter(np.repeat([m],len(md_m)),md_m['accuracy_train'], label=str(m))
plt.savefig('validation_curves/combined_max_depth.png')
'''
def getTaggerScores(files):
    all_taggers_scores = {}
    all_taggers_positions = {}
    scores = {}
    max_score = 0
    max_id = ""
    for f in files:
        print 'plotting file: ' + f
        with open('evaluationObjects/'+f,'r') as p:
            model = pickle.load(p)
        taggers = model.taggers
        for t in taggers:
            if t not in all_taggers_scores.keys():
                all_taggers_scores[t] = []
                all_taggers_positions[t] = []
        feature_importances = 100.0 * (model.feature_importances / model.feature_importances.max())

        #feature_importances = model.feature_importances
        #print feature_importances
        sorted_idx = np.argsort(feature_importances)[::-1]
        for x in range(len(sorted_idx)):
            # add the feature importance score and the position
            all_taggers_scores[taggers[sorted_idx[x]]].append(feature_importances[sorted_idx[x]])
            all_taggers_positions[taggers[sorted_idx[x]]].append(x)
            #print 'feature: ' + taggers[sorted_idx[x]]+' score: ' + str(feature_importances[sorted_idx[x]])

        scores[model.job_id] = model.ROC_rej_power_05
        #print 'sum: ' + str(np.sum(feature_importances))
        if model.score > max_score:
            max_score = model.ROC_rej_power_05
            max_id = model.job_id



    print 'max score: ' + str(max_score)
    print 'id: ' + max_id

    sorted_scores = sorted(scores.items(), key=operator.itemgetter(1))
    print sorted_scores
    print 'press a key to get tagger scores'
    # write to file
    f = open('scores'+key+'.txt','w')
    for s in sorted_scores:
        f.write(s[0] + ': ' + str(s[1]) +'\n')
    f.close()

    # find median of each
    medians = {}
    for t in all_taggers_scores.keys():
        #print t
        #print np.median(all_taggers_scores[t])
        medians[t] = np.median(all_taggers_scores[t])
        #print '***********'

        #print all_taggers_scores


        # sort medians

    sorted_medians = sorted(medians.items(), key=operator.itemgetter(1))
    print sorted_medians
