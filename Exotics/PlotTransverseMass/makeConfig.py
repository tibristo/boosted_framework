import sys
import subprocess
import re
basecfg = sys.argv[2]
inputalgorithms = sys.argv[1]
name = sys.argv[3]
version = sys.argv[4]

algs = open(inputalgorithms)
algorithms = []
alg_prefix = {}

def stripPrefix(a):
    nums = re.findall(r'\d+', a)
    offset = 4 # remove the jet_ at the beginning
    if a.find('LC') != -1:
        offset+=2 #removes the LC
    if a.find('CamKt') != -1:
        return a[offset+len('CamKt')+len(nums[0]):], str('CamKt'+nums[0])
    elif a.find('AntiKt') != -1:
        return a[offset+len('AntiKt')+len(nums[0]):], str('AntiKt'+nums[0])

for a in algs:
    algostrip, prefix = stripPrefix(a)
    algorithms.append(algostrip.strip())
    alg_prefix[algostrip.strip()] = prefix
    a = a[4:]
    outname = 'config/'+name+'_'+algostrip.strip()+'.cfg'
    #print 'cp ' + basecfg + ' ' + outname
    subprocess.call('cp ' + basecfg + ' ' + outname,shell=True)
    repstr = 's/algorithm\ =/algorithm\ =\ '+str(a.strip())+'/g'
    #print 'grep -rl algorithm ' + outname+' | xargs sed -i \''+repstr+'\''
    subprocess.call('grep -rl algorithm ' + outname+' | xargs sed -i \''+repstr+'\'', shell=True)
    fileid = '_'+name + '_'+version 
    repfileid = 's/fileid\ =/fileid\ =\ '+fileid+'/g'
    #print 'grep -rl fileid ' + outname+' | xargs sed -i \''+ repfileid+'\''
    subprocess.call('grep -rl fileid ' + outname+' | xargs sed -i \''+ repfileid+'\'', shell=True)
    # add branches.txt and config file names to cfg
    br = name+'_'+algostrip.strip()+'_branches.txt'
    conf = name+'_'+algostrip.strip()+'.config'
    #print 'echo \'branches-file = '+br+ '\' >> '+outname
    subprocess.call('echo \'branches-file = branches/'+br+ '\' >> '+outname,shell = True)

# this does what createAlgorithmBranches.py does
numbins = [] # give all 20 bins?
axisvalues = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[-2,2],[0,1],[0,100],[0,1],[0,0.3],[0,90],[0,30],[0,15],[-2,5],[-2,5],[-2,10],[-2,5],[0.1,1],[0,1],[0,150],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[-1500,100],[0,1]] #going to be 2d for min, max
var_type = ['E','pt','m','eta','phi','emfrac','Tau1','Tau2','Tau3','WIDTH','SPLIT12','SPLIT23','SPLIT34','Dip12','Dip13','Dip23','DipExcl12','PlanarFlow','Angularity','QW','PullMag','PullPhi','Pull_C00','Pull_C01','Pull_C10','Pull_C11','constit_index']


subjet_alg = {'TopoTrimmedPtFrac5SmallR30':'TopoTrimmedSubjetsPtFrac5SmallR30','TopoSplitFilteredMu67SmallR0YCut9':'TopoSplitFiltSubjetsMu67SmallR0YCut9','TopoSplitFilteredMu100SmallR30YCut0':'TopoSplitFiltSubjetsMu100SmallR30YCut0', 'TopoSplitFilteredMu100SmallR30YCut4':'TopoSplitFiltSubjetsMu100SmallR30YCut4', 'TopoSplitFilteredMu100SmallR30YCut9':'TopoSplitFiltSubjetsMu100SmallR30YCut9', 'TopoSplitFilteredMu100SmallR30YCut12':'TopoSplitFiltSubjetsMu100SmallR30YCut12', 'TopoSplitFilteredMu100SmallR30YCut15':'TopoSplitFiltSubjetsMu100SmallR30YCut15'}
subjet_axes = [[0,1700],[0,1200],[0,300],[-5,5],[-4,4],[0,1],[0,0.3]]
subjet_vars = ['E','pt','m','eta','phi','constit_index','WIDTH']


for alg in algorithms:
    f = open(name+'_'+alg+'.config','w')
    f2 = open('branches/'+name+'_'+alg+'_branches.txt','w')
    tcount = 0
    prefix = 'jet_' + alg_prefix[alg]
    if alg.find('LCTopo') == -1:
        prefix+='LC'
    for t in var_type:
        #changes numbins[tcount] to 20 for now....
        
        if t != 'emfrac':
            f.write('truth_'+t+','+'jet_CamKt12Truth_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')
            f2.write('jet_CamKt12Truth_'+t+'\n')

        f.write('topo_'+t+','+prefix+'Topo_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')        
        f.write('groomed_'+t+','+prefix+alg+'_'+t+','+str(20) + ',' + str(axisvalues[tcount][0]) + ','+str(axisvalues[tcount][1])+',\n')
        
        f2.write(prefix+'Topo_'+t+'\n')
        f2.write(prefix+alg+'_'+t+'\n')
        tcount+=1
    ecount = 0
    if alg in subjet_alg.keys():
        for e in subjet_vars:
            f.write('subjet_'+e+','+prefix+subjet_alg[alg]+'_'+e+','+str(20)+','+str(subjet_axes[ecount][0]) + ',' + str(subjet_axes[ecount][1])+',\n')    
            f2.write(prefix+subjet_alg[alg]+'_'+e+'\n')
            ecount+=1

    f.close()
    f2.close()

