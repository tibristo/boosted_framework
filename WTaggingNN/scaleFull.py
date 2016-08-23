import ROOT as rt
import pandas as pd
import numpy as np
import root_numpy as rn
import os
import sys
import create_folds as cf
import pickle
#import create_folds.scalerNN

#training_files = ['folds/'+f for f in os.listdir('folds') if f.find('train')!=-1]
#training_files = ['/Disk/ds-sopa-group/PPE/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz6_nTrk_v1_1300_1800_mw_merged.root']
training_files = ['/Disk/ds-sopa-group/PPE/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_nonTrk_v3_400_1200_mw_merged.root']
#training_files = ['/Disk/ds-sopa-group/PPE/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_nTrk_v1_800_1200_mw_merged.root']
cols = np.linspace(1,42,42,dtype=int)
#filename = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz6_nTrk_v1_1300_1800_mw_merged.csv'
filename = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_nonTrk_v3_400_1200_mw_merged.csv'
#filename = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/BosonTagging/csv/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_nTrk_v1_800_1200_mw_merged.csv'

for f in training_files:
    #tfile = rt.TFile.Open(f)
    #tree = tfile.Get('outputTree')
    # get the branches and only focus on the ones not listed below
    # remove the variables that are "observers", ie that do not get used for training, or weights, since those must not be scaled.
    #br = tree.GetListOfBranches()
    X = rn.root2rec(f)
    variables = list(X.dtype.names)
    
    #for b in br:
    #    variables.append(br.GetName())
    variables.remove('label')
    observers = ['mc_event_weight','jet_antikt10truthtrimmedptfrac5smallr20_pt','jet_antikt10truthtrimmedptfrac5smallr20_eta','m','pt','eta','phi','evt_xsec','evt_filtereff','evt_nevts','weight','jet_camkt12truth_pt','jet_camkt12truth_eta','jet_camkt12truth_phi','jet_camkt12truth_m','jet_camkt12lctopo_pt','jet_camkt12lctopo_eta','jet_camkt12lctopo_phi','jet_camkt12lctopo_m','eff','averageintperxing']
    for o in observers:
        if o in variables:
            variables.remove(o)
    curr_means = np.zeros(len(variables))
    curr_std = np.ones(len(variables))
    weighted_means = np.zeros(len(variables))
    weighted_std = np.ones(len(variables))
    for j,v in enumerate(variables):
        mean = np.mean(X[v])
        print mean
        curr_means[j] = mean
        weighted_mean = cf.weightedMean(X[v], X['weight'])
        weighted_means[j] = weighted_mean
        std = np.std(X[v])
        curr_std[j] = std
        weighted_std[j] = cf.weightedStd(X[v], X['weight'], weighted_mean)
        # create a scaler object to scale datasets later on
    sc = cf.scalerNN(variables, curr_means, curr_std, weighted_means, weighted_std)
    # pickle the standardscaler and the variables for scaling
    #sc_filename = os.path.abspath(prefix+name+'train' + suffix % i)
    sc_filename = os.path.abspath('test.root')
    #sc_filename = sc_filename.replace('.root','_scaler.pkl')
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

    
        
    X_full = np.recfromcsv(filename,usecols=cols)
    xf = list(X_full.dtype.names)
    print xf
    for i,v in enumerate(variables):
        mean = curr_means[i]
        std = curr_std[i]
        v_low = v.lower()
        X_full[v_low] = (X_full[v_low]-mean)/std
    #full_file = os.path.abspath('folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz6_nTrk_v1_1300_1800_mw_mergedscalefull.root')
    full_file = os.path.abspath('folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_nonTrk_v3_400_1200_mw_mergedscalefull.root')
    #full_file = os.path.abspath('folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_nTrk_v1_800_1200_mw_mergedscalefull.root')
    rn.array2root(X_full, full_file, 'outputTree','recreate')
    #raw_input()

