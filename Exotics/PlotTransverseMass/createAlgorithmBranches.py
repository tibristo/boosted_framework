

algorithms = ['TopoSplitFilteredMu67SmallR0YCut9','TopoSplitFilteredMu100SmallR30YCut4','TopoTrimmedPtFrac5SmallR30','TopoTrimmedPtFrac5SmallR20','TopoPrunedCaRcutFactor50Zcut10','TopoPrunedCaRcutFactor50Zcut20','AntiKt2LCTopo','AntiKt3LCTopo','AntiKt4LCTopo']
numbins = [] # give all 20 bins?
axisvalues = [[0,1700],[0,1200],[0,120],[-5,5],[-4,4],[-2,2],[0,1],[0,100],[0,1],[0,0.3],[0,90],[0,30],[0,15],[-2,5],[-2,5],[-2,10],[-2,5],[-1500,100],[-1500,100],[0,150],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100]] #going to be 2d for min, max
var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','Tau3','WIDTH','SPLIT12','SPLIT23','SPLIT34','Dip12','Dip13','Dip23','DipExcl12','PlanarFlow','Angularity','QW','PullMag','PullPhi','Pull_C00','Pull_C01','Pull_C10','Pull_C11']

for alg in algorithms:
    f = open(alg+'.config','w')
    f2 = open(alg+'_branches.txt','w')
    tcount = 0
    prefix = 'jet_CamKt12'
    if alg.find('LCTopo') == -1:
        prefix+='LC'
    for t in var_type:
        #changes numbins[tcount] to 20 for now....
        f.write(t+','+'jet_CamKt12Truth_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+'\n')
        f.write(t+','+'jet_CamKt12LCTopo_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+'\n')        
        f.write(t+','+prefix+alg+'_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+'\n')

        f2.write('jet_CamKt12Truth_'+t+'\n')
        f2.write('jet_CamKt12Topo_'+t+'\n')
        f2.write(prefix+alg+'_'+t+'\n')
        tcount+=1
    f.close()
    f2.close()
