import ROOT as rt
import math
import sys
rt.gROOT.LoadMacro("include/TLorentzVectorDict.h+")

# calculate dR
def deltaR(eta1, phi1, eta2, phi2):
    deltaPhi = abs(phi1-phi2)
    deltaEta = eta1-eta2
    if(deltaPhi > math.pi):
        deltaPhi = 2*math.pi - deltaPhi
    return math.sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi)


file_in_name = sys.argv[1]
file_out_name = sys.argv[2]
signal = False
if len(sys.argv) == 4:
    if sys.argv[3].lower() == 'signal':
        signal = True

# define the branches we're going to write out
b_topo_pt = rt.std.vector(float)()
b_topo_eta = rt.std.vector(float)()
b_topo_phi = rt.std.vector(float)()
b_topo_m = rt.std.vector(float)()

b_truth_pt = rt.std.vector(float)()
b_truth_eta = rt.std.vector(float)()
b_truth_phi = rt.std.vector(float)()
b_truth_m = rt.std.vector(float)()
b_truth_id = rt.std.vector(int)()

b_nTracks = rt.std.vector(int)()
b_mc_channel_number = rt.std.vector(int)()


file_in = rt.TFile.Open(file_in_name,'READ')

tree = file_in.Get('outputTree')

entries = tree.GetEntries()

file_out = rt.TFile(file_out_name+'.root','RECREATE')
tree_out = rt.TTree('outputTree','outputTree')

tree_out.Branch('topo_pt', b_topo_pt)
tree_out.Branch('topo_eta', b_topo_eta)
tree_out.Branch('topo_phi', b_topo_phi)
tree_out.Branch('topo_m', b_topo_m)
tree_out.Branch('truth_pt', b_truth_pt)
tree_out.Branch('truth_eta', b_truth_eta)
tree_out.Branch('truth_phi', b_truth_phi)
tree_out.Branch('truth_m', b_truth_m)
tree_out.Branch('truth_id', b_truth_id)
tree_out.Branch('nTracks', b_nTracks)
tree_out.Branch('mc_channel_number', b_mc_channel_number)

csv_out = open(file_out_name+'.csv','write')
csv_out.write('topo_pt,topo_eta,topo_phi,topo_m,truth_pt,truth_eta,truth_phi,truth_m,truth_id,nTracks,mc_channel_number\n')

for i in range(entries):
    tree.GetEntry(i)
    if i % 10000 == 0:
        print 'entry: ' + str(i)
    # clear all of the vectors
    b_topo_pt.clear()
    b_topo_eta.clear()
    b_topo_phi.clear()
    b_topo_m.clear()
    b_truth_pt.clear()
    b_truth_eta.clear()
    b_truth_phi.clear()
    b_truth_m.clear()
    b_truth_id.clear()
    b_nTracks.clear()
    b_mc_channel_number.clear()

    # find the leading topo jet
    max_pt = 0
    max_idx = 0
    if tree.jet_AntiKt10LCTopo_pt.size() == 0:
        continue
    for x, p in enumerate(tree.jet_AntiKt10LCTopo_pt):
        if p > max_pt:
            max_pt = p
            max_idx = x
    topo_pt = tree.jet_AntiKt10LCTopo_pt[max_idx]/1000.0
    topo_eta = tree.jet_AntiKt10LCTopo_eta[max_idx]
    topo_phi = tree.jet_AntiKt10LCTopo_phi[max_idx]
    topo_m = tree.jet_AntiKt10LCTopo_m[max_idx]/1000.0
    nTracks = int(tree.jet_AntiKt10LCTopo_nTracks[max_idx])
    # match this with a truth boson parent
    # truth boson info is a 4vec
    truth_id = 99
    truth_pt = -99
    truth_eta = -99
    truth_phi = -99
    truth_m = -99
    if signal:
        truth_matched = False
        truth_idx = 0
        for x, t in enumerate(tree.truthBosons):
            pt = t.Pt()
            eta = t.Eta()
            phi = t.Phi()
            if not truth_matched and deltaR(topo_eta, topo_phi, eta, phi) and pt > 5000 and abs(eta) < 6:
                # check the id
                truth_id = abs(tree.truthBoson_ID[x])
                truth_pt = tree.truthBosons[x].Pt()/1000.0
                truth_eta = tree.truthBosons[x].Eta()
                truth_phi = tree.truthBosons[x].Phi()
                truth_m = tree.truthBosons[x].M()/1000.0
                                                
                truth_matched = True
            
    if signal and not truth_matched:
        continue

    b_truth_pt.push_back(truth_pt)
    b_truth_eta.push_back(truth_eta)
    b_truth_phi.push_back(truth_phi)
    b_truth_m.push_back(truth_m)
    b_truth_id.push_back(truth_id)

    b_topo_pt.push_back(topo_pt)
    b_topo_eta.push_back(topo_eta)
    b_topo_phi.push_back(topo_phi)
    b_topo_m.push_back(topo_m)
    b_nTracks.push_back(nTracks)

    b_mc_channel_number.push_back(tree.mc_channel_number)
    
    tree_out.Fill()
    csv_out.write(str(topo_pt)+','+str(topo_eta)+','+str(topo_phi)+','+str(topo_m)+','+str(truth_pt)+','+str(truth_eta)+','+str(truth_phi)+','+str(truth_m)+','+str(truth_id)+','+str(nTracks)+','+str(tree.mc_channel_number)+'\n')
    
csv_out.close()
tree_out.Write()
file_out.Close()
