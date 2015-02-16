import ROOT
import sys
file_in = ROOT.TFile(sys.argv[1])
physics = file_in.Get("physics")
ob = physics.GetListOfBranches()

if (len(sys.argv) <= 2):
    f = open(sys.argv[1][:-5]+"_algorithms.txt",'w')
else:
    f = open(sys.argv[2]+"_algorithms.txt",'w')
for tb in ob:
    if tb.GetName().startswith('jet_') and tb.GetName().endswith('_m') and tb.GetName().find('emscale') == -1 and tb.GetName().find('constscale') == -1 and (tb.GetName().find('Trim')!=-1 or tb.GetName().find('Prune')!= -1 or tb.GetName().find('Split') != -1) and tb.GetName().find('Sub') == -1:
        f.write(tb.GetName()[:-2]+'\n')
f.close()
