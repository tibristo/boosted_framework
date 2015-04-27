import sys
import subprocess
import re

# usage: python makeConfig.py list_of_algorithms baseconfig_file output_name version_number_in_output prefix_of_branches_file

inputalgorithms = sys.argv[1]
basecfg = sys.argv[2]
name = sys.argv[3]
version = sys.argv[4]
branches_prefix = sys.argv[5]

algs = open(inputalgorithms)
algorithms = []
alg_prefix = {}

# get the prefix of the file, like AntiKtX
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
    #algostrip, prefix = stripPrefix(a)
    # remove the jet_ at the beginning
    a = a.strip()
    if a.find('jet_') != -1:
        a = a.replace('jet_','')
    algorithms.append(a)
    #alg_prefix[algostrip.strip()] = prefix
    outname = 'config/'+name+'_'+version+'_'+a+'.cfg'
    #print 'cp ' + basecfg + ' ' + outname
    subprocess.call('cp ' + basecfg + ' ' + outname,shell=True)
    repstr = 's/algorithm\ =/algorithm\ =\ '+str(a)+'/g'
    #print 'grep -rl algorithm ' + outname+' | xargs sed -i \''+repstr+'\''
    subprocess.call('grep -rl algorithm ' + outname+' | xargs sed -i \''+repstr+'\'', shell=True)
    fileid = '_'+name + '_'+version 
    repfileid = 's/fileid\ =/fileid\ =\ '+fileid+'/g'
    #print 'grep -rl fileid ' + outname+' | xargs sed -i \''+ repfileid+'\''
    subprocess.call('grep -rl fileid ' + outname+' | xargs sed -i \''+ repfileid+'\'', shell=True)
    # add branches.txt and config file names to cfg
    if branches_prefix == '':
        br = name+'_'+a+'_branches.txt'
    else:
        br = branches_prefix+'_'+a+'_branches.txt'
    conf = name+'_'+a+'.config'
    #print 'echo \'branches-file = '+br+ '\' >> '+outname
    subprocess.call('echo \'branches-file = branches/'+br+ '\' >> '+outname,shell = True)
