from ROOT import *
import sys
from array import array
import copy

f = TFile.Open(sys.argv[1])

t = f.Get("outputTree")

fo = TFile.Open(sys.argv[1]+".2.root","recreate")

to = t.CloneTree(0)
'''
elx = array('f', [0.0,0.0])
ely = array('f', [0.0,0.0])
elz = array('f', [0.0,0.0])
elt = array('f', [0.0,0.0])
mux = array('f', [0.0,0.0])
muy = array('f', [0.0,0.0])
muz = array('f', [0.0,0.0])
mut = array('f', [0.0,0.0])
'''

elx = std.vector(float)()
ely = std.vector(float)()
elz = std.vector(float)()
elt = std.vector(float)()
mux = std.vector(float)()
muy = std.vector(float)()
muz = std.vector(float)()
mut = std.vector(float)()


to.Branch("electronX",elx)#,"electronX/F")
to.Branch("electronY",ely)#,"electronY/F")
to.Branch("electronZ",elz)#,"electronZ/F")
to.Branch("electronT",elt)#,"electronT/F")
to.Branch("muonX",mux)#,"muonX/F")
to.Branch("muonY",muy)#,"muonY/F")
to.Branch("muonZ",muz)#,"muonZ/F")
to.Branch("muonT",mut)#,"muonT/F")

for x in range(0,t.GetEntries()):
    t.GetEntry(x)
    if x %1000 == 0:
        print str(x) + '/' + str(t.GetEntries())
    #for i in range(len(elx)):
        #elx.pop();
        #elx[i] = -999.0
        #ely[i] = -999.0;
        #elz[i] = -999.0;
        #elt[i] = -999.0;
        #mux[i] = -999.0;
        #muy[i] = -999.0;
        #muz[i] = -999.0;
        #mut[i] = -999.0;
    elx.clear()
    ely.clear()
    elz.clear()
    elt.clear()
    mux.clear()
    muy.clear()
    muz.clear()
    mut.clear()

    for e in range(t.electrons.size()):
        #elx.append(copy.deepcopy(t.electrons[e].X()))
        if e > 1:
            continue
        #elx[e] = t.electrons[e].X()
        #ely[e] = t.electrons[e].Y()
        #elz[e] = t.electrons[e].Z()
        #elt[e] = t.electrons[e].T()
        elx.push_back(t.electrons[e].X())
        ely.push_back(t.electrons[e].Y())
        elz.push_back(t.electrons[e].Z())
        elt.push_back(t.electrons[e].T())

    for m in range(t.muons.size()):
        if m > 1:
            continue
        #mux[m] = t.muons[m].X()
        #muy[m] = t.muons[m].Y()
        #muz[m] = t.muons[m].Z()
        #mut[m] = t.muons[m].T()
        mux.push_back(t.muons[m].X())
        muy.push_back(t.muons[m].Y())
        muz.push_back(t.muons[m].Z())
        mut.push_back(t.muons[m].T())

    to.Fill()

to.Write()
fo.Close()
