

<<<<<<< HEAD
<<<<<<< HEAD
algorithms = ['CamKt12LCTopoSplitFilteredMu67SmallR0YCut9','CamKt12LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut9', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut12','CamKt12LCTopoSplitFilteredMu100SmallR30YCut15','AntiKt6LCTopoTrimmedPtFrac5SmallR20','AntiKt4LCTopoTrimmedPtFrac5SmallR20','AntiKt10LCTopoTrimmedPtFrac5SmallR20','AntiKt10LCTopoTrimmedPtFrac5SmallR30','AntiKt8LCTopoTrimmedPtFrac5SmallR20','AntiKt12LCTopoTrimmedPtFrac5SmallR20','CamKt6LCTopoTrimmedPtFrac5SmallR20','CamKt8LCTopoTrimmedPtFrac5SmallR20','CamKt10LCTopoTrimmedPtFrac5SmallR20','CamKt12LCTopoTrimmedPtFrac5SmallR20','CamKt6LCTopoPrunedCaRcutFactor50Zcut15','CamKt6LCTopoPrunedCaRcutFactor50Zcut10','CamKt8LCTopoPrunedCaRcutFactor50Zcut15','CamKt8LCTopoPrunedCaRcutFactor50Zcut10','CamKt10LCTopoPrunedCaRCutFactor50Zcut15','CamKt10LCTopoPrunedCaRCutFactor50Zcut10','CamKt12LCTopoPrunedCaRcutFactor50Zcut15','CamKt12LCTopoPrunedCaRcutFactor50Zcut10', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut9', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut9','CamKt10LCTopoPrunedCaRCutFactor50Zcut15','CamKt10LCTopoPrunedCaRCutFactor50Zcut10']
numbins = [] # give all 20 bins?
axisvalues = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[-2,2],[0,1],[0,100],[0,1],[0,0.3],[0,90],[0,30],[0,15],[-2,5],[-2,5],[-2,10],[-2,5],[0.1,1],[0,1],[0,150],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[0,1]] #going to be 2d for min, max
#var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','Tau21','SPLIT12','Dip12','PlanarFlow','Angularity','constit_index','constit_n','Aplanarity','TauWTA1','TauWTA2','TauWTA2TauWTA1','Sphericity','ThrustMaj', 'ThrustMin', 'ZCUT12','FoxWolfram20','EEC_C1','EEC_C2','EEC_D1','EEC_D2','SoftDrop','config_massFraction','config_ktycut2', 'FoxWolfram_0', 'FoxWolfram_2', 'SoftDropTag']
var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','Tau21','SPLIT12','Dip12','PlanarFlow','Angularity','Aplanarity','TauWTA1','TauWTA2','TauWTA2TauWTA1','Sphericity','ThrustMaj', 'ThrustMin', 'ZCUT12','EEC_C2_1','EEC_C2_2','EEC_D2_1','EEC_D2_2','ECF1','ECF2','ECF3','ECF1_beta2','ECF2_beta2','ECF3_beta2', 'FoxWolfram20', 'FoxWolfram2', 'Mu12','config_massFraction','config_ktycut2','YFilt','evt_sumWeights','evt_kfactor','evt_filtereff','evt_nEvts','evt_sumWeights','evt_xsecs']

subjet_alg = {'AntiKt4LCTopoTrimmedPtFrac5SmallR20':'AntiKt4LCTopoTrimmedSubjetsPtFrac5SmallR20','AntiKt6LCTopoTrimmedPtFrac5SmallR20':'AntiKt6LCTopoTrimmedSubjetsPtFrac5SmallR20','AntiKt8LCTopoTrimmedPtFrac5SmallR20':'AntiKt8LCTopoTrimmedSubjetsPtFrac5SmallR20', 'AntiKt10LCTopoTrimmedPtFrac5SmallR20':'AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR20','AntiKt10LCTopoTrimmedPtFrac5SmallR30':'AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR30', 'AntiKt12LCTopoTrimmedPtFrac5SmallR20':'AntiKt12LCTopoTrimmedSubjetsPtFrac5SmallR20',\
'CamKt12LCTopoSplitFilteredMu67SmallR0YCut9':'CamKt12LCTopoSplitFiltSubjetsMu67SmallR0YCut9','CamKt12LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut9', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut12':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut12', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut15':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut15',\
'CamKt6LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut9',\
'CamKt8LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut9'}


subjet_axes = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[0,1],[0,0.3]]
subjet_vars = ['E','pt','m','eta','phi','constit_index','WIDTH']

#loop through algorithms and add to branches file
for alg in algorithms:
    # open output file
    print alg
    f2 = open('branches/13tev_beta2_'+alg+'_branches.txt','w')
    tcount = 0
    # branch prefix
    prefix = 'jet_'
    # if LC in prefix name, keep it there
    if alg.find('LC') != -1:
        algprefix = alg[:alg.find('LC')+2]
    else: # if there is no LC in the name
        algprefix = alg[:alg.find('Topo')]
        print alg
    # xAODs have Jets appended to the algorithm name
    jetsadd =''
    if alg.find('Jets') != -1:
        jetsadd = 'Jets'
    # loop through all of teh variables
    for t in var_type:        
        #changes numbins[tcount] to 20 for now....        
        # emfrac doesn't exist for truth
        # don't write config_massFraction and config_ktycut2 for topo and truth
        if t.find('config') == -1:
            if t != 'emfrac':
                #f2.write('jet_CamKt12Truth'+jetsadd+'_'+t+'\n')
                f2.write(prefix+alg.replace("LCTopo","Truth")+jetsadd+'_'+t+'\n')
            # topo + groomed

            f2.write(prefix+algprefix+'Topo'+jetsadd+'_'+t+'\n')
        f2.write(prefix+alg+'_'+t+'\n')
        tcount+=1
    ecount = 0
    # now loop through and add the subjet branches
    if alg in subjet_alg.keys():
        for e in subjet_vars:
            f2.write(prefix+subjet_alg[alg]+jetsadd+'_'+e+'\n')
            ecount+=1
    f2.close()

=======
algorithms = ['TopoSplitFilteredMu67SmallR0YCut9','TopoSplitFilteredMu100SmallR30YCut4', 'TopoSplitFilteredMu100SmallR30YCut0', 'TopoSplitFilteredMu100SmallR30YCut9', 'TopoSplitFilteredMu100SmallR30YCut12','TopoSplitFilteredMu100SmallR30YCut15','TopoTrimmedPtFrac5SmallR30','TopoTrimmedPtFrac5SmallR20','TopoPrunedCaRcutFactor50Zcut10','TopoPrunedCaRcutFactor50Zcut20','AntiKt2LCTopo','AntiKt3LCTopo','AntiKt4LCTopo', 'AntiKt10LCTopoPrunedCaRcutFactor50Zcut10Jets','AntiKt10LCTopoSplitFilteredMu100SmallR30YCut4Jets', 'AntiKt10LCTopoTrimmedPtFrac5SmallR30Jets', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut4Jets','AntiKt10LCTopoSplitFilteredMu100SmallR30YCut4Jets']
=======
algorithms = ['CamKt12LCTopoSplitFilteredMu67SmallR0YCut9','CamKt12LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut9', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut12','CamKt12LCTopoSplitFilteredMu100SmallR30YCut15','AntiKt6LCTopoTrimmedPtFrac5SmallR20','AntiKt4LCTopoTrimmedPtFrac5SmallR20','AntiKt10LCTopoTrimmedPtFrac5SmallR20','AntiKt10LCTopoTrimmedPtFrac5SmallR30','AntiKt8LCTopoTrimmedPtFrac5SmallR20','AntiKt12LCTopoTrimmedPtFrac5SmallR20','CamKt6LCTopoTrimmedPtFrac5SmallR20','CamKt8LCTopoTrimmedPtFrac5SmallR20','CamKt10LCTopoTrimmedPtFrac5SmallR20','CamKt12LCTopoTrimmedPtFrac5SmallR20','CamKt6LCTopoPrunedCaRcutFactor50Zcut15','CamKt6LCTopoPrunedCaRcutFactor50Zcut10','CamKt8LCTopoPrunedCaRcutFactor50Zcut15','CamKt8LCTopoPrunedCaRcutFactor50Zcut10','CamKt10LCTopoPrunedCaRCutFactor50Zcut15','CamKt10LCTopoPrunedCaRCutFactor50Zcut10','CamKt12LCTopoPrunedCaRcutFactor50Zcut15','CamKt12LCTopoPrunedCaRcutFactor50Zcut10', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut9', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut9','CamKt10LCTopoPrunedCaRCutFactor50Zcut15','CamKt10LCTopoPrunedCaRCutFactor50Zcut10']
>>>>>>> 393929b3c6507b0be568335e02ce19e45373b4b3
numbins = [] # give all 20 bins?
axisvalues = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[-2,2],[0,1],[0,100],[0,1],[0,0.3],[0,90],[0,30],[0,15],[-2,5],[-2,5],[-2,10],[-2,5],[0.1,1],[0,1],[0,150],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[0,1]] #going to be 2d for min, max
#var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','Tau21','SPLIT12','Dip12','PlanarFlow','Angularity','constit_index','constit_n','Aplanarity','TauWTA1','TauWTA2','TauWTA2TauWTA1','Sphericity','ThrustMaj', 'ThrustMin', 'ZCUT12','FoxWolfram20','EEC_C1','EEC_C2','EEC_D1','EEC_D2','SoftDrop','config_massFraction','config_ktycut2', 'FoxWolfram_0', 'FoxWolfram_2', 'SoftDropTag']
var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','Tau21','SPLIT12','Dip12','PlanarFlow','Angularity','Aplanarity','TauWTA1','TauWTA2','TauWTA2TauWTA1','Sphericity','ThrustMaj', 'ThrustMin', 'ZCUT12','EEC_C2_1','EEC_C2_2','EEC_D2_1','EEC_D2_2','ECF1','ECF2','ECF3','ECF1_beta2','ECF2_beta2','ECF3_beta2', 'FoxWolfram20', 'FoxWolfram2', 'Mu12','config_massFraction','config_ktycut2','YFilt','evt_sumWeights','evt_kfactor','evt_filtereff','evt_nEvts','evt_sumWeights','evt_xsecs']

subjet_alg = {'AntiKt4LCTopoTrimmedPtFrac5SmallR20':'AntiKt4LCTopoTrimmedSubjetsPtFrac5SmallR20','AntiKt6LCTopoTrimmedPtFrac5SmallR20':'AntiKt6LCTopoTrimmedSubjetsPtFrac5SmallR20','AntiKt8LCTopoTrimmedPtFrac5SmallR20':'AntiKt8LCTopoTrimmedSubjetsPtFrac5SmallR20', 'AntiKt10LCTopoTrimmedPtFrac5SmallR20':'AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR20','AntiKt10LCTopoTrimmedPtFrac5SmallR30':'AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR30', 'AntiKt12LCTopoTrimmedPtFrac5SmallR20':'AntiKt12LCTopoTrimmedSubjetsPtFrac5SmallR20',\
'CamKt12LCTopoSplitFilteredMu67SmallR0YCut9':'CamKt12LCTopoSplitFiltSubjetsMu67SmallR0YCut9','CamKt12LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut9', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut12':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut12', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut15':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut15',\
'CamKt6LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut9',\
'CamKt8LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut9'}


subjet_axes = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[0,1],[0,0.3]]
subjet_vars = ['E','pt','m','eta','phi','constit_index','WIDTH']

#loop through algorithms and add to branches file
for alg in algorithms:
    # open output file
    print alg
    f2 = open('branches/13tev_beta2_'+alg+'_branches.txt','w')
    tcount = 0
    # branch prefix
    prefix = 'jet_'
    # if LC in prefix name, keep it there
    if alg.find('LC') != -1:
        algprefix = alg[:alg.find('LC')+2]
    else: # if there is no LC in the name
        algprefix = alg[:alg.find('Topo')]
        print alg
    # xAODs have Jets appended to the algorithm name
    jetsadd =''
    if alg.find('Jets') != -1:
        jetsadd = 'Jets'
    # loop through all of teh variables
    for t in var_type:        
        #changes numbins[tcount] to 20 for now....        
        # emfrac doesn't exist for truth
        # don't write config_massFraction and config_ktycut2 for topo and truth
        if t.find('config') == -1:
            if t != 'emfrac':
                #f2.write('jet_CamKt12Truth'+jetsadd+'_'+t+'\n')
                f2.write(prefix+alg.replace("LCTopo","Truth")+jetsadd+'_'+t+'\n')
            # topo + groomed

            f2.write(prefix+algprefix+'Topo'+jetsadd+'_'+t+'\n')
        f2.write(prefix+alg+'_'+t+'\n')
        tcount+=1
    ecount = 0
    # now loop through and add the subjet branches
    if alg in subjet_alg.keys():
        for e in subjet_vars:
            f2.write(prefix+subjet_alg[alg]+jetsadd+'_'+e+'\n')
            ecount+=1
    f2.close()

<<<<<<< HEAD

'''
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_n
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_E
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_pt
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_m
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_eta
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_phi
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_constit_n
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_constit_index
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_Lepton_n
jet_AntiKt10TruthTrimmedSubjetsPtFrac5SmallR30_WIDTH


jet_AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR30_n

jet_CamKt12TruthSplitFiltSubjetsMu67SmallR0YCut9_n

jet_CamKt12TruthSplitFiltSubjetsMu100SmallR30YCut4_n

jet_CamKt12LCTopoSplitFiltSubjetsMu67SmallR0YCut9_n

jet_CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut4_n
'''
#jet_AntiKt10Truth_TrimmedSubjetsPtFrac5SmallR30_n

#jet_AntiKt10LCTopo_TrimmedSubjetsPtFrac5SmallR30_n

#jet_CamKt12Truth_SplitFiltSubjetsMu67SmallR0YCut9_n

#jet_CamKt12Truth_SplitFiltSubjetsMu100SmallR30YCut4_n
>>>>>>> master
=======
>>>>>>> 393929b3c6507b0be568335e02ce19e45373b4b3
