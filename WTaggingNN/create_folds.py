# A number of tools for doing MVA analysis.
# Includes:
# Cross validation 
# Grid search
# Learning curves
# Probability distributions
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cross_validation import ShuffleSplit, StratifiedKFold
import os, sys, time
import numpy.lib.recfunctions as nf
from pprint import pprint
from root_numpy import array2root

#import cv_fold
def persist_cv_splits(X, y, observers, n_cv_iter=5, test_size=0.25, name='data', prefix='folds/', suffix="_cv_%03d.root", random_state=None, scale=False):
    """Materialize randomized train test splits of a dataset."""
    cv = StratifiedKFold(y,n_folds=n_cv_iter,shuffle=True)
    #cv = ShuffleSplit(X.shape[0], n_iter=n_cv_iter,
    #    test_size=test_size, random_state=random_state)
    cv_split_filenames = []
    #scale the data
    
    for i, (train, test) in enumerate(cv):
        print ("iteration %03d" % i)
        # can't scale all the variables - we don't want to scale weights!
        if scale:
            scaler = StandardScaler()
            #mean = X[train].mean()
            #std = X[train].std()
            #Xtrain = (X[train] - mean)/std
            #Xtest = (X[test] - mean)/std
            Xtrain = scaler.fit_transform(X[train])
            Xtest = scaler.transform(X[test])
        else:
            Xtrain = X[train]
            Xtest = X[test]
        ytrain = y[train]
        ytest = y[test]
        cv_split_filename = prefix+name + 'train' + suffix % i
        cv_split_train = os.path.abspath(cv_split_filename)
        cv_split_test = os.path.abspath(cv_split_filename.replace('train','test'))

        cv_split_filenames.append([cv_split_train, cv_split_test])
        xtrain_o = [Xtrain,observers[train]]
        xtest_o = [Xtest,observers[test]]
        merged_train = nf.merge_arrays(xtrain_o,flatten=True,usemask=True)
        merged_test = nf.merge_arrays(xtest_o,flatten=True,usemask=True)
        #print merged.shape
        #print merged.dtype.names
        rectrain = nf.append_fields(merged_train, names='label', data=ytrain, dtypes=np.int32, usemask=False)#, dtypes=int)#, usemask=False)
        #rec = nf.append_fields(X[train], 'label', y[train], usemask=False)
        array2root(rectrain, cv_split_train, 'outputTree', 'recreate')
        rectest = nf.append_fields(merged_test, names='label', data=ytest, dtypes=np.int32,usemask=False)#, dtypes=int)#, usemask=False)
        array2root(rectest, cv_split_test, 'outputTree', 'recreate')
    return cv_split_filenames


def cross_validation(data,iterations, name='data'):
    variables = list(data.dtype.names)
    variables.remove('label')
    # remove the variables that are "observers", ie that do not get used for training, or weights, since those must not be scaled.
    observers = ['mc_event_weight','jet_antikt10truthtrimmedptfrac5smallr20_pt','jet_antikt10truthtrimmedptfrac5smallr20_eta','m','pt','eta','phi','evt_xsec','evt_filtereff','evt_nevts','weight','jet_camkt12truth_pt','jet_camkt12truth_eta','jet_camkt12truth_phi','jet_camkt12truth_m','jet_camkt12lctopo_pt','jet_camkt12lctopo_eta','jet_camkt12lctopo_phi','jet_camkt12lctopo_m']
    for o in observers:
        variables.remove(o)
    X = data[variables]
    y = data['label']
    observer_data = data[observers]
    #print observer_data.dtype.names
    #raw_input()
    #w = data['weight']

    filenames = persist_cv_splits(X, y, observer_data, n_cv_iter=iterations, name=name, suffix="_cv_%03d.root", test_size=0.25, random_state=None, scale=False)
    #all_parameters, all_tasks = grid_search(
     #   lb_view, model, filenames, params)
    return filenames
    #return all_parameters, all_tasks

def plotFiles(filenames, variables):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    print filenames
    for f in filenames:
        print f
        f_open = ROOT.TFile.Open(f[0])
        tree = f_open.Get('outputTree')
        leg = ROOT.TLegend(0.8,0.55,0.9,0.65);leg.SetFillColor(ROOT.kWhite)
        c = ROOT.TCanvas(f[0])
        # create histograms for each variable
        hists = {}
        for v in variables:
            hist_sig_name = 'sig_'+v
            hist_bkg_name = 'bkg_'+v
            varexp = v+">>"+hist_sig_name
            cutstring = 'label==0'
            tree.Draw(varexp, cutstring)
            hist_sig = ROOT.gDirectory.Get(hist_sig_name).Clone()
            # normalise
            if hist_sig.Integral()!=0:
                hist_sig.Scale(1/hist_sig.Integral())
            hist_sig.SetLineColor(2)
            hist_sig.SetTitle('Signal')
            hist_sig.GetXaxis().SetTitle(v)
            tree.Draw(varexp.replace('sig','bkg'),cutstring.replace('0','1'))
            hist_bkg = ROOT.gDirectory.Get(hist_bkg_name).Clone()
            if hist_bkg.Integral()!=0:
                hist_bkg.Scale(1/hist_bkg.Integral())
            hist_bkg.SetLineColor(4)
            hist_bkg.SetTitle('Background')
            hist_sig.GetXaxis().SetTitle(v)
            sig_max =  hist_sig.GetBinContent(hist_sig.GetMaximumBin())
            bkg_max =  hist_sig.GetBinContent(hist_bkg.GetMaximumBin())

            leg.Clear()
            leg.AddEntry(hist_sig, 'Signal','l'); leg.AddEntry(hist_bkg, 'Background','l')
            

            if sig_max > bkg_max:
                hist_sig.Draw()
                hist_bkg.Draw('same')
            else:
                hist_bkg.Draw()
                hist_sig.Draw('same')
            leg.Draw('same')
            c.SaveAs(f[0].replace('.root','.png'))
            #maxX = max(hist_sig.GetXaxis().GetXmax(),hist_bkg.GetXaxis().GetXmax())

print os.getcwd()
from collections import OrderedDict

path = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/BosonTagging/csv/'
algorithm = 'AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_500_1000_nomw_merged'
#name = sys.argv[1].replace('.csv','')
#if name.find('/')!=-1:
#    name = name[name.rfind('/')+1:]
#print name
name = algorithm
cols = np.linspace(1,43,43,dtype=int)
#data = np.recfromcsv(sys.argv[1],usecols=cols)
data = np.recfromcsv(path+algorithm+'.csv',usecols=cols)

filenames = cross_validation(data, 4, name)

# plot all of the variables in these files
variables = list(data.dtype.names)
plotFiles(filenames, variables)
