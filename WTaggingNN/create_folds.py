# Tools for creating folds of data for cross validation
# Persists the folds onto the disk
# Can scale the data
# Can plot the variables
#

import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cross_validation import ShuffleSplit, StratifiedKFold
import os, sys, time
import argparse
import numpy.lib.recfunctions as nf
from pprint import pprint
import math
from root_numpy import array2root
import ROOT
import AtlasStyle as atlas
import copy

class scalerNN:
    def __init__(self, variables, all_means, all_std, weighted_means = None, weighted_std = None):
        '''
        store the scaling for a sample.  Store the means and std deviations. Create standard scaler out of this too.
        
        keywords args:
        variables --- the variables that the means and std correspond to. all_means/std[0] -> variables[0]
        all_means --- np array of means(value) for variables(key)
        all_std --- np array of std dev(value) for variables(key)
        '''
        import numpy as np
        from sklearn.preprocessing import StandardScaler
        self.means = all_means
        self.std = all_std
        self.variables = variables
        
        self.scaler = StandardScaler()
        self.scaler.means_ = all_means
        self.scaler.std_ = all_std

        if weighted_means is not None:
            self.weighted_means = weighted_means
        if weighted_std is not None:
            self.weighted_std = weighted_std
        if weighted_std is not None and weighted_means is not None:
            self.weighted_scaler = StandardScaler()
            self.weighted_scaler.means_ = weighted_means
            self.weighted_scaler.std_ = weighted_std

def weightedMean(X, w):
    '''
    Method to calculate the weighted mean of a sample
    '''
    sum_wx = 0.0
    for i in range(len(X)):
        sum_wx += w[i]*X[i]

    sum_wx /= np.sum(w)
    return sum_wx

def weightedStd(X, w, u):
    '''
    Method to calculate the weighted std of a sample
    '''
    variance = 0.0
    for i in range(len(X)):
        variance += w[i]*(X[i]-u)**2
    variance /= np.sum(w)

    return math.sqrt(variance)

def scaleSample(scaler_filename, filename='/Disk/ecdf-nfs-ppe/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw_merged.csv', prefix='folds/', name='data'):
    '''
    Method to standardise an input file with the mean and std from the training data (this is stored in the scaler file).
    '''
    #import scalerNN
    import pickle
    # check that the scaler pickle file exists
    if not os.path.isfile(scaler_filename):
        print 'Scaler file does not exist'
        return
    
    # open the scaler pickle
    with open(scaler_filename,'r') as p:
        scaler = pickle.load(p)
    # check the type to make sure it's not an error
    if isinstance(scaler, basestring):
        print "Scaler appears to be an error message"
        print scaler
        return

    # need to know how many columns there are in the data
    f = open(filename)
    l = f.readline()
    colcount = l.count(',')
    f.close()
    
    cols = np.linspace(1,colcount,colcount,dtype=int)
    #data = np.recfromcsv(sys.argv[1],usecols=cols)
    X_full = np.recfromcsv(filename, usecols=cols)
    # here we can pull the variables from the scaler and just scale those. the rest are all observers.
    for i,v in enumerate(scaler.variables):
        mean = scaler.means[i]
        std = scaler.std[i]
        X_full[v] = (X_full[v]-mean)/std

    #recfull = nf.append_fields(X_full, names='label', data=y, dtypes=np.int32, usemask=False)#, dtypes=int)#, usemask=False)
    full_file = os.path.abspath(prefix+name + 'full.root')
    array2root(X_full, full_file, 'outputTree', 'recreate')    
    

#import cv_fold
def persist_cv_splits(X, y, w, variables, observers, n_cv_iter=5, test_size=0.25, name='data', prefix='folds/', suffix="_cv_%03d.root", random_state=None, scale=False):
    import pickle
    #import scalerNN
    """Materialize randomized train test splits of a dataset."""
    cv = StratifiedKFold(y,n_folds=n_cv_iter,shuffle=True)
    #cv = ShuffleSplit(X.shape[0], n_iter=n_cv_iter,
    #    test_size=test_size, random_state=random_state)
    cv_split_filenames = []
    # normalise the weights, otherwise agilepack doesnt work
    #observers['weight'] = 1/np.sum(observers['weight'])

    for i, (train, test) in enumerate(cv):
        print ("iteration %03d" % i)
        # can't scale all the variables - we don't want to scale weights!
        Xtrain_weighted = X[train]
        Xtest_weighted = X[test]
        if scale:
            # the standardscaler does not seem to work on numpy structured arrays
            # a quick fix, which is possibly the cleanest way, is to just loop through
            # each variable and scale it manually.  I don't think this should differ from the StandardScaler()
            # it is just more code.
            Xtrain = X[train]#*observers['weight'][train]
            Xtest = X[test]#*observers['weight'][test]
            # set up some empty numpy arrays for the means and std
            curr_means = np.zeros(len(variables))
            curr_std = np.ones(len(variables))
            weighted_means = np.zeros(len(variables))
            weighted_std = np.ones(len(variables))
            # calculate the stats for each variable of interest
            for j,v in enumerate(variables):
                mean = np.mean(X[v][train])
                curr_means[j] = mean
                weighted_mean = weightedMean(X[v][train], w[train])
                weighted_means[j] = weighted_mean
                std = np.std(X[v][train])
                curr_std[j] = std
                weighted_std[j] = weightedStd(X[v][train], w[train], weighted_mean)

                # do the standardisation
                print v
                Xtrain[v] = (Xtrain[v]-mean)/std
                Xtest[v] = (Xtest[v]-mean)/std
                
            # create a scaler object to scale datasets later on
            sc = scalerNN(variables, curr_means, curr_std, weighted_means, weighted_std)
            # pickle the standardscaler and the variables for scaling
            sc_filename = os.path.abspath(prefix+name+'train' + suffix % i)
            sc_filename = sc_filename.replace('.root','_scaler.pkl')
            try:
                with open(sc_filename,'w') as d:
                    pickle.dump(sc,d)
                d.close()
            except:
                msg = 'unable to dump ' + sc_filename
                msg+= str(sys.exc_info()[0])
                with open(sc_filename,'w') as d:
                    pickle.dump(msg,d)
                d.close()
        else:
            # if we're not standardising
            Xtrain = X[train]
            Xtest = X[test]
        ytrain = y[train]
        ytest = y[test]
        # output file names
        cv_split_filename = prefix+name + 'train' + suffix % i
        cv_split_train = os.path.abspath(cv_split_filename)
        cv_split_test = os.path.abspath(cv_split_filename.replace('train','test'))

        # merge all of the training data, observers, labels and weights.
        cv_split_filenames.append([cv_split_train, cv_split_test])
        xtrain_o = [Xtrain,observers[train]]
        xtest_o = [Xtest,observers[test]]
        merged_train = nf.merge_arrays(xtrain_o,flatten=True,usemask=True)
        merged_test = nf.merge_arrays(xtest_o,flatten=True,usemask=True)

        # write this to a root file
        rectrain = nf.append_fields(merged_train, names='label', data=ytrain, dtypes=np.int32, usemask=False)
        array2root(rectrain, cv_split_train, 'outputTree', 'recreate')
        rectest = nf.append_fields(merged_test, names='label', data=ytest, dtypes=np.int32,usemask=False)
        array2root(rectest, cv_split_test, 'outputTree', 'recreate')

        # now do the weighted ones
        '''
        xtrain_w = [Xtrain_weighted,observers[train]]
        xtest_w = [Xtest_weighted,observers[test]]
        merged_train_w = nf.merge_arrays(xtrain_w,flatten=True,usemask=True)
        merged_test_w = nf.merge_arrays(xtest_w,flatten=True,usemask=True)
        rectrain_w = nf.append_fields(merged_train_w, names='label', data=ytrain, dtypes=np.int32, usemask=False)#, dtypes=int)#, usemask=False)
        #rec = nf.append_fields(X[train], 'label', y[train], usemask=False)
        array2root(rectrain_w, cv_split_train.replace('train','train_w'), 'outputTree', 'recreate')
        rectest_w = nf.append_fields(merged_test_w, names='label', data=ytest, dtypes=np.int32,usemask=False)#, dtypes=int)#, usemask=False)
        array2root(rectest_w, cv_split_test.replace('test','test_w'), 'outputTree', 'recreate')
        '''

    return cv_split_filenames


def cross_validation(data,iterations, name='data', scale=True, pt_rw = True):
    # name of all variables in the dataset
    variables = list(data.dtype.names)
    # remove the ones we do not want to standardise
    variables.remove('label')
    # remove the variables that are "observers", ie that do not get used for training, or weights, since those must not be scaled.
    observers = ['mc_event_weight','jet_antikt10truthtrimmedptfrac5smallr20_pt','jet_antikt10truthtrimmedptfrac5smallr20_eta','m','pt','eta','phi','evt_xsec','evt_filtereff','evt_nevts','weight','jet_camkt12truth_pt','jet_camkt12truth_eta','jet_camkt12truth_phi','jet_camkt12truth_m','jet_camkt12lctopo_pt','jet_camkt12lctopo_eta','jet_camkt12lctopo_phi','jet_camkt12lctopo_m','eff','averageintperxing']

    # variables in the mc15 samples
    #mc_event_weight,jet_AntiKt10TruthTrimmedPtFrac5SmallR20_pt,jet_AntiKt10TruthTrimmedPtFrac5SmallR20_eta,EEC_D2_1,Aplanarity,m,FoxWolfram20,ThrustMaj,TauWTA2,EEC_C2_1,Angularity,ThrustMin,pt,PlanarFlow,Dip12,Mu12,YFilt,TauWTA1,phi,EEC_D2_2,SPLIT12,Sphericity,ZCUT12,eta,TauWTA2TauWTA1,EEC_C2_2,averageIntPerXing,evt_xsec,evt_filtereff,evt_nEvts,jet_CamKt12Truth_pt,jet_CamKt12Truth_eta,jet_CamKt12Truth_phi,jet_CamKt12Truth_m,jet_CamKt12LCTopo_pt,jet_CamKt12LCTopo_eta,jet_CamKt12LCTopo_phi,jet_CamKt12LCTopo_m,weight,eff,label
    
    # this can be done in one line
    for o in observers:
        if o.lower() in variables:
            variables.remove(o)
    #print variables
    # Get X, which is all training/ testing variables
    X = data[variables]
    # target variable
    y = data['label']
    # weights
    # sometimes the weights need to be adjusted.

    w = data['weight']
    # variables that are not standardised
    observer_data = nf.append_fields(data[observers], names='weight_train', data=copy.deepcopy(data['weight']), dtypes=np.float, usemask=False)
    #observer_data = data[observers]
    if not pt_rw:
        # for signal all weights set to 1, no pt rw for training the nn
        print 'setting signal weights to 1'
        observer_data['weight_train'][data['label']==1] = 1.0
    # create the folds
    filenames = persist_cv_splits(X, y, w, variables, observer_data, n_cv_iter=iterations, name=name, suffix="_cv_%03d.root", test_size=0.25, random_state=None, scale=scale)

    return filenames


def plotFiles(filenames, variables, key, weight_plots = False, plot_dict = {}, tex_dict = {}):
    # plot a bunch of files and get stats of variables

    ROOT.gROOT.SetBatch(True)
    atlas.SetAtlasStyle()
    print filenames

    # check how many cv folds are represented in this list of files
    cv_folds = filter(lambda x: x.find('cv') != -1 and x.endswith('root'), filenames)
    # get the cv_xxx name and only return unique values
    cv_nums = list(set(map(lambda x: 'cv_'+x.split('_')[-1].replace('.root','') , cv_folds)))
    
    # create a dictionary for the stats of each cv fold
    stats = {}
    # create a dict for the event counts of each cv fold
    event_counts = {}
    for c in cv_nums:
        stats[c] = {'Train':{},'Valid':{}}
        event_counts[c] = {'Train':{},'Valid':{}}
    # dict for the stats of the full dataset
    stats['Full'] = {}
    event_counts['Full'] = {}

    # add this to the end of the filenames to differentiate between weighted and not weighted
    weight_id = '_weighted' if weight_plots else ''
    
    for i, f in enumerate(filenames):
        # is this a test or train file?
        file_type = "Full"
        if f.lower().find('train') != -1:
            file_type = "Train"
        elif f.lower().find('test') != -1:
            file_type = 'Valid'
        # what cv split is it? look for cv_xyz.root
        cv_num = ''
        if file_type != 'Full':
            cv_num = 'cv_'+f.split('_')[-1].replace('.root','')

        # open the file and get the tree
        f_open = ROOT.TFile.Open('folds/'+f)
        tree = f_open.Get('outputTree')
        leg = ROOT.TLegend(0.7,0.55,0.9,0.75);leg.SetFillColor(ROOT.kWhite)
        c = ROOT.TCanvas(f)

        total_events = tree.GetEntries()
        signal_events = tree.GetEntries("label==1")
        bkg_events = tree.GetEntries("label==0")

        # this is the full dataset
        if cv_num != '':
            event_counts[cv_num][file_type] = '{0:15}  {1:10} {2:14}{3:10}'.format(file_type+' ' + cv_num,str(signal_events), str(bkg_events),str(total_events))
        else:
            # this is a cv split
            event_counts['Full'] = '{0:15}  {1:10} {2:14}{3:10}'.format(file_type,str(signal_events), str(bkg_events),str(total_events))
            
        # create histograms for each variable
        hists = {}

        for v in variables:
            # set up the names for the different histograms
            hist_full_name = v
            hist_sig_name = 'sig_'+v
            hist_bkg_name = 'bkg_'+v

            # first get the full histogram to get an idea of the combined mean and rms
            tree.Draw(v+'>>'+v)
            # pull from the global space
            hist_full = ROOT.gDirectory.Get(hist_full_name).Clone()
            mean = '{0:.4f}'.format(float(hist_full.GetMean()))
            std = '{0:.4f}'.format(float(hist_full.GetRMS()))
            
            # set up the variable expression that gets used in the Draw function
            varexp = v+">>"+hist_sig_name
            # cut string to select signal only.  An additional string can be added here to apply the weights.
            cutstring = '(label==1)'#*(weight)
            if weight_plots:
                cutstring += '*weight'
            # create the signal histogram and retrieve it
            tree.Draw(varexp, cutstring)
            hist_sig = ROOT.gDirectory.Get(hist_sig_name).Clone()
            # stats
            sig_mean = '{0:.4f}'.format(float(hist_sig.GetMean()))
            sig_std = '{0:.4f}'.format(float(hist_sig.GetRMS()))
            # normalise
            if hist_sig.Integral()!=0:
                hist_sig.Scale(1/hist_sig.Integral())
            # set some drawing options and titles
            hist_sig.SetLineColor(2)
            hist_sig.SetTitle('Signal')
            hist_sig.GetXaxis().SetTitle(v)
            # now get the background histogram
            tree.Draw(varexp.replace('sig','bkg'),cutstring.replace('1','0'))
            hist_bkg = ROOT.gDirectory.Get(hist_bkg_name).Clone() # pull from global
            # stats
            bkg_mean = '{0:.4f}'.format(float(hist_sig.GetMean()))
            bkg_std = '{0:.4f}'.format(float(hist_sig.GetRMS()))
            # normalise
            if hist_bkg.Integral()!=0:
                hist_bkg.Scale(1/hist_bkg.Integral())
            # drawing options and titles
            hist_bkg.SetLineColor(4)
            hist_bkg.SetTitle('Background')
            hist_sig.GetXaxis().SetTitle(v)
            

            leg.Clear()
            # add the legend entries
            leg.AddEntry(hist_sig, 'Signal','l'); leg.AddEntry(hist_bkg, 'Background','l')            

            # find the maximum for when we draw them together on a single canvas
            sig_max =  hist_sig.GetMaximum()
            bkg_max =  hist_bkg.GetMaximum()
            max_val = max(sig_max, bkg_max)
            hist_sig.SetMaximum(max_val*1.1)
            hist_bkg.SetMaximum(max_val*1.1)

            if v.strip() in plot_dict.keys():
                hist_sig.GetXaxis().SetTitle(plot_dict[v])
                hist_bkg.GetXaxis().SetTitle(plot_dict[v])
                fnum = 'Full' if cv_num == '' else cv_num
                hist_sig.SetTitle('Fold ' +fnum+': ' + plot_dict[v])
                hist_bkg.SetTitle('Fold '+fnum + ': '+plot_dict[v])
            
            hist_sig.Draw()
            hist_bkg.Draw('same')

            leg.Draw('same')
            if cv_num == '':
                c.SaveAs('fold_plots/'+key+'_Full_'+v+weight_id+'.png')
            else:
                c.SaveAs('fold_plots/'+key+'_Full_'+v+weight_id+'.png')
            # write the means and std to the stats file
            result = '{0:15}: {1:10} {2:10} {3:10} {4:10} {5:10} {6:10}'.format(file_type+' '+cv_num,str(mean),str(std),str(sig_mean),str(sig_std),str(bkg_mean),str(bkg_std))
            # check that this variable has a dictionary entry
            #if v not in stats[cv_num][file_type].keys():
            #stats[cv_num][file_type][v] = {}
            if cv_num != '': # full dataset
                stats[cv_num][file_type][v] = result
            else: # cv split
                stats['Full'][v] = result
                
    # now all of the stats can be written to file!
    # first write the combined stats
    combined_stats = open('fold_stats/combined_stats_'+key+weight_id+'.txt','w')
    combined_stats.write('{0:15}  {1:10} {2:14}{3:10}'.format('Sample','Signal','Background','Total')+'\n')
    combined_stats.write(event_counts['Full']+'\n')
    # write out all of the cv splits
    for cv in cv_nums:
        for f in ['Train','Valid']:
            combined_stats.write(event_counts[cv][f]+'\n')
    combined_stats.write('\n'+'\n')
    # now start doing the variables
    combined_stats.write('\n{0:15}: {1:10} {2:10} {3:10} {4:10} {5:10} {6:10}'.format('Variable','Mean','Std','Mean Sig','Std Sig','Mean Bkg','Std Bkg')+'\n\n')

    for v in variables:
        if v.strip() in tex_dict.keys():
            combined_stats.write(tex_dict[v]+'\n')
        else:
            combined_stats.write(v+'\n')
        combined_stats.write(stats['Full'][v]+'\n')
        for c in cv_nums:
            combined_stats.write(stats[c]['Train'][v]+'\n')
            combined_stats.write(stats[c]['Valid'][v]+'\n')
        combined_stats.write('\n')
    combined_stats.close()
    
    # write the stats for each cv fold
    for cv in cv_nums:
        stats_file = open('fold_stats/'+key+'_'+cv+weight_id+'.txt','w')
        stats_file.write('{0:15}  {1:10} {2:14}{3:10}'.format('Sample','Signal','Background','Total')+'\n')
        stats_file.write(event_counts['Full']+'\n')
        stats_file.write(event_counts[cv]['Train']+'\n')
        stats_file.write(event_counts[cv]['Valid']+'\n')

        stats_file.write('\n{0:15}: {1:10} {2:10} {3:10} {4:10} {5:10} {6:10}'.format('Variable','Mean','Std','Mean Sig','Std Sig','Mean Bkg','Std Bkg')+'\n\n')
        # now write each variable
        for v in variables:
            if v.strip() in tex_dict.keys():
                stats_file.write(tex_dict[v]+'\n')
            else:
                stats_file.write(v+'\n')
            stats_file.write(stats['Full'][v]+'\n')
            stats_file.write(stats[cv]['Train'][v]+'\n')
            stats_file.write(stats[cv]['Valid'][v]+'\n')
            stats_file.write('\n')
        stats_file.close()

def main(args):
    # print the current directory
    print os.getcwd()
    from collections import OrderedDict
    # set up the arguments parser
    parser = argparse.ArgumentParser(description='Plot some variables or create cv folds.')
    parser.add_argument('folds', help = 'If the cv folds should be created')
    parser.add_argument('plot', help = 'If the cv folds should be plotted')
    parser.add_argument('--key', help = 'Key to be used for identifying files')
    parser.add_argument('--weight', help = 'If the plots should be weighted')
    parser.add_argument('--algorithm', default = 'AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw_merged', help = 'Name of the algorithm (this is the name of the csv file, without the .csv at the end).')
    parser.add_argument('--fulldataset', default = 'DEFAULT', help = 'Full dataset filename.  This is the NOT Cleaned file.')
    parser.add_argument('--ptrw', dest='ptrw', action='store_true',help = 'If signal should be pt weighted. Default True.')
    parser.add_argument('--no-ptrw', dest='ptrw',action='store_false', help = 'If signal should be pt weighted. Default True.')
    parser.set_defaults(ptrw=True)
    # parse args
    args = parser.parse_args()
    
    if not args.folds or not args.plot:
        print 'need to set if folds should be created! usage: python create_folds.py folds(true/false) plot(true/false) [--key=key] [--weight=true/false] [--algorithm=alg] [--fulldataset=fullset] [--ptrw=true/false]'
        sys.exit(0)

    # this is the default path
    path = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/BosonTagging/csv/'

    name = args.algorithm+'scale' # use this if we want to name the output file something different
    # set up the cols which get used for creating a dataframe from csv
    colcheck = open(path+args.algorithm+'.csv')
    line1 = colcheck.readline()
    colcount = line1.count(',')
    cols = np.linspace(1,colcount, colcount,dtype=int)
    print 'pt rw ' +str(args.ptrw)

    #['mc_event_weight', 'jet_antikt10truthtrimmedptfrac5smallr20_pt', 'jet_antikt10truthtrimmedptfrac5smallr20_eta', 'aplanarity', 'thrustmin', 'tau1', 'sphericity', 'm', 'foxwolfram20', 'tau21', 'thrustmaj', 'eec_c2_1', 'pt', 'eec_c2_2', 'dip12', 'split12', 'phi', 'tauwta2tauwta1', 'eec_d2_1', 'yfilt', 'mu12', 'tauwta2', 'zcut12', 'angularity', 'tau2', 'eec_d2_2', 'eta', 'tauwta1', 'planarflow', 'averageintperxing', 'evt_xsec', 'evt_filtereff', 'evt_nevts', 'jet_camkt12truth_pt', 'jet_camkt12truth_eta', 'jet_camkt12truth_phi', 'jet_camkt12truth_m', 'jet_camkt12lctopo_pt', 'jet_camkt12lctopo_eta', 'jet_camkt12lctopo_phi', 'jet_camkt12lctopo_m', 'weight', 'eff', 'label']

    # if we want to create the folds
    if args.folds.lower() == 'true':
        # read in the data as a numpy recarray
        data = np.recfromcsv(path+args.algorithm+'.csv',usecols=cols)
        # get the names of the variables in the recarray
        variables = list(data.dtype.names)
    
        print 'Creating folds'
        # are we standardising the data?
        scale = True
        # the cross validation method will call the persists_cv method and create the folds
        filenames = cross_validation(data, 4, name, scale, pt_rw = args.ptrw)
        #full_dataset = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw_merged.csv'
        # name of the full dataset which is used for the cv splits
        if args.fulldataset == 'DEFAULT':
            full_dataset = path+algorithm.replace('loose','notcleaned')+'.csv'
            if not os.path.isfile(full_dataset):
                print 'full dataset does not exist!' + full_dataset
        else:
            full_dataset = path + args.fulldataset

        if scale:
            # create the full datasets scaled according to the different training samples:
            for f in filenames:
                # first entry in f is the train sample name, second is test
                # get the scalerNN object for the train sample
                # scalerNN saves the values we used to do the standardisation
                scaler_fname = f[0].replace('.root', '_scaler.pkl')
                # find out which cv split we're on
                cv_split = f[0][f[0].find('cv'):f[0].find('.root')]
                # apply the standardisation to the dataset
                scaleSample(scaler_fname, filename=full_dataset, prefix='folds/', name=full_dataset[full_dataset.rfind('/')+1:].replace('.csv','')+'_full_scaled_'+cv_split)
    else:
        print 'Not creating folds'        
        
    
    # arg because I didn't read this properly, if using plotFiles() after doing cross_validation, the filenames list
    # must be converted to a 1D list
    
    # plot all of the variables in these files
    if args.plot.lower() == 'true':
        # key to be used for finding the correct files
        key = '13tev_matchedM_loose_v2_200_1000_mw'
        if args.key:
            key = args.key

        print 'Plotting files'
        try:
            filenames
        except NameError:
            print 'Filenames were not defined, creating a list from folds directory using key ' + key
            filenames = [f for f in os.listdir('folds') if f.endswith('root') and f.find(key) != -1]
        else:
            print 'Filenames were defined, converting to a 1D list'
            filenames = [item for sublist in filenames for item in sublist]
        # the variables we're interested in                  dc14
        #variables = ['aplanarity','eec_c2_1', 'eec_c2_2', 'split12','eec_d2_1', 'eec_d2_2', 'tauwta2tauwta1','zcut12','sphericity','mu12','planarflow']
        # mc15
        variables = ['aplanarity','eec_c2_1', 'split12','eec_d2_1', 'tauwta2tauwta1','zcut12','sphericity','mu12','planarflow']
        plot_dict = {'tauwta2tauwta1':"#tau^{WTA}_{2}/#tau^{WTA}_{1}",'eec_c2_1':"C^{(#beta=1)}_{2}",'eec_c2_2':"C^{(#beta=2)}_{2}",'eec_d2_1':"D^{(#beta=1)}_{2}",'eec_d2_2':"D^{(#beta=2)}_{2}", 'split12':"#sqrt{d_{12}}",'aplanarity':"#it{A}",'zcut12':"#sqrt{z_{12}}",'sphericity':"#it{S}",'planarflow':"#it{P}"}
        tex_dict = {'tauwta2tauwta1':r"$\tau^{WTA}_{2}/\tau^{WTA}_{1}$",'eec_c2_1':r"$C^{(\beta=1)}_{2}$",'eec_c2_2':r"$C^{(\beta=2)}_{2}$",'eec_d2_1':r"$D^{(\beta=1)}_{2}$",'eec_d2_2':r"$D^{(\beta=2)}_{2}$", 'split12':r"$\sqrt{d_{12}}$",'aplanarity':r"$\textit{A}$",'zcut12':r"$\sqrt{z_{12}}$",'sphericity':r"$\textit{S}$",'planarflow':r"$\textit{P}$"}
        weight = True if args.weight.lower() == 'true' else False

        plotFiles(filenames, variables, key, weight_plots=weight, plot_dict = plot_dict, tex_dict = tex_dict)
        
if __name__ == "__main__":
    print ' running main '
    main(sys.argv)
