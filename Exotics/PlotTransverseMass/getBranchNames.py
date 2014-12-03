import ROOT
import sys
file_in = ROOT.TFile(sys.argv[1])
physics = file_in.Get("physics")
ob = physics.GetListOfBranches()
if (len(sys.argv) <= 2):
    f = open(sys.argv[1][:-5]+"_branches.txt",'w')
else:
    f = open(sys.argv[2]+"_branches.txt",'w')
for tb in ob:
    f.write(tb.GetName()+'\n')
f.close()
