

algorithms = ['CamKt12LCTopoSplitFilteredMu67SmallR0YCut9','CamKt12LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut9', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut12','CamKt12LCTopoSplitFilteredMu100SmallR30YCut15','AntiKt10LCTopoTrimmedPtFrac5SmallR20','AntiKt8LCTopoTrimmedPtFrac5SmallR20','AntiKt12LCTopoTrimmedPtFrac5SmallR20','CamKt6LCTopoTrimmedPtFrac5SmallR20','CamKt8LCTopoTrimmedPtFrac5SmallR20','CamKt10LCTopoTrimmedPtFrac5SmallR20','CamKt12LCTopoTrimmedPtFrac5SmallR20','CamKt6LCTopoPrunedCaRcutFactor50Zcut15','CamKt6LCTopoPrunedCaRcutFactor50Zcut10','CamKt8LCTopoPrunedCaRcutFactor50Zcut15','CamKt8LCTopoPrunedCaRcutFactor50Zcut10','CamKt10LCTopoPrunedCaRcutFactor50Zcut15','CamKt10LCTopoPrunedCaRcutFactor50Zcut10','CamKt12LCTopoPrunedCaRcutFactor50Zcut15','CamKt12LCTopoPrunedCaRcutFactor50Zcut10', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut9', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut4', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut0', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut9']
numbins = [] # give all 20 bins?
axisvalues = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[-2,2],[0,1],[0,100],[0,1],[0,0.3],[0,90],[0,30],[0,15],[-2,5],[-2,5],[-2,10],[-2,5],[0.1,1],[0,1],[0,150],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[0,1]] #going to be 2d for min, max
var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','SPLIT12','Dip12','PlanarFlow','Angularity','constit_index','constit_n','Aplanarity','TauWTA1','TauWTA2','Sphericity','ThrustMaj', 'ThrustMin', 'ZCUT12','FoxWolfram20','EEC_C1','EEC_C2','EEC_D1','EEC_D2','SoftDrop','QJetsVol','config_massFraction','config_ktycut2']

subjet_alg = {'AntiKt8LCTopoTrimmedPtFrac5SmallR20':'AntiKt8LCTopoTrimmedSubjetsPtFrac5SmallR20', 'AntiKt10LCTopoTrimmedPtFrac5SmallR20':'AntiKt10LCTopoTrimmedSubjetsPtFrac5SmallR20', 'AntiKt12LCTopoTrimmedPtFrac5SmallR20':'AntiKt12LCTopoTrimmedSubjetsPtFrac5SmallR20',\
'CamKt12LCTopoSplitFilteredMu67SmallR0YCut9':'CamKt12LCTopoSplitFiltSubjetsMu67SmallR0YCut9','CamKt12LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut9', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut12':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut12', 'CamKt12LCTopoSplitFilteredMu100SmallR30YCut15':'CamKt12LCTopoSplitFiltSubjetsMu100SmallR30YCut15',\
'CamKt6LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt6LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt6LCTopoSplitFiltSubjetsMu100SmallR30YCut9',\
'CamKt8LCTopoSplitFilteredMu100SmallR30YCut0':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut0', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut4':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut4', 'CamKt8LCTopoSplitFilteredMu100SmallR30YCut9':'CamKt8LCTopoSplitFiltSubjetsMu100SmallR30YCut9'}


subjet_axes = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[0,1],[0,0.3]]
subjet_vars = ['E','pt','m','eta','phi','constit_index','WIDTH']

#loop through algorithms and add to branches file
for alg in algorithms:
    # open output file
    print alg
    f2 = open('branches/'+alg+'_branches.txt','w')
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
                f2.write('jet_CamKt12Truth'+jetsadd+'_'+t+'\n')
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

