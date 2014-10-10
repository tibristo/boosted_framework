

algorithms = ['TopoSplitFilteredMu67SmallR0YCut9','TopoSplitFilteredMu100SmallR30YCut4','TopoTrimmedPtFrac5SmallR30','TopoTrimmedPtFrac5SmallR20','TopoPrunedCaRcutFactor50Zcut10','TopoPrunedCaRcutFactor50Zcut20','AntiKt2LCTopo','AntiKt3LCTopo','AntiKt4LCTopo']
numbins = [] # give all 20 bins?
axisvalues = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[-2,2],[0,1],[0,100],[0,1],[0,0.3],[0,90],[0,30],[0,15],[-2,5],[-2,5],[-2,10],[-2,5],[0.1,1],[0,1],[0,150],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[0,1]] #going to be 2d for min, max
var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','Tau3','WIDTH','SPLIT12','SPLIT23','SPLIT34','Dip12','Dip13','Dip23','DipExcl12','PlanarFlow','Angularity','QW','PullMag','PullPhi','Pull_C00','Pull_C01','Pull_C10','Pull_C11','constit_index']

subjet_alg = {'TopoTrimmedPtFrac5SmallR30':'TopoTrimmedSubjetsPtFrac5SmallR30','TopoSplitFilteredMu67SmallR0YCut9':'TopoSplitFiltSubjetsMu67SmallR0YCut9','TopoSplitFilteredMu100SmallR0YCut4':'TopoSplitFiltSubjetsMu100SmallR0YCut4'}
alg_prefix = {'TopoSplitFilteredMu67SmallR0YCut9':'CamKt12','TopoSplitFilteredMu100SmallR30YCut4':'CamKt12','TopoTrimmedPtFrac5SmallR30':'AntiKt10','TopoTrimmedPtFrac5SmallR20':'AntiKt10','TopoPrunedCaRcutFactor50Zcut10':'AntiKt10','TopoPrunedCaRcutFactor50Zcut20':'AntiKt10','AntiKt2LCTopo':'','AntiKt3LCTopo':'','AntiKt4LCTopo':''}
subjet_axes = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[0,1],[0,0.3]]
subjet_vars = ['E','pt','m','eta','phi','constit_index','WIDTH']

for alg in algorithms:
    f = open(alg+'.config','w')
    f2 = open(alg+'_branches.txt','w')
    tcount = 0
    prefix = 'jet_' + alg_prefix[alg]
    if alg.find('LCTopo') == -1:
        prefix+='LC'
    for t in var_type:
        #changes numbins[tcount] to 20 for now....
        
        if t != 'emfrac':
            f.write('truth_'+t+','+'jet_CamKt12Truth_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')
            #f.write('truth_'+t+','+'jet_CamKt12Truth'+alg[4:]+'_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')
            f2.write('jet_CamKt12Truth_'+t+'\n')
            #f2.write('jet_CamKt12Truth'+alg[4:]+'_'+t+'\n')
        #f.write('topo_'+t+','+'jet_CamKt12LCTopo_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')        
        f.write('topo_'+t+','+prefix+'Topo_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')        
        f.write('groomed_'+t+','+prefix+alg+'_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')
        


        #f2.write('jet_CamKt12LCTopo_'+t+'\n')
        f2.write(prefix+'Topo_'+t+'\n')
        #f2.write(prefix+'_'+t+'\n')
        f2.write(prefix+alg+'_'+t+'\n')
        tcount+=1
    ecount = 0
    if alg in subjet_alg:
        for e in subjet_vars:
            f.write('subjet_'+e+','+prefix+subjet_alg[alg]+'_'+e+','+str(20)+','+str(subjet_axes[ecount][0]) + ',' + str(subjet_axes[ecount][1])+',\n')    
            f2.write(prefix+subjet_alg[alg]+'_'+e+'\n')
            ecount+=1

    f.close()
    f2.close()


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
