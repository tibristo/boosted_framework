import ROOT as rt
import math
import numpy as np
rt.gROOT.LoadMacro("include/TLorentzVectorDict.h+")
def DeltaR(eta1, phi1, eta2, phi2):
    dphi = abs(phi1-phi2)
    deta = eta1-eta2
    if dphi > rt.TMath.Pi():
        dphi = rt.TMath.TwoPi() - dphi
    return rt.TMath.Sqrt(deta*deta + dphi*dphi)

tc = rt.TChain('outputTree')
tc.Add('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/VVXAOD-00-00-22/mc15_13TeV.301257.Pythia8EvtGen_A14NNPDF23LO_Wprime_WZqqqq_m1000.merge.DAOD_EXOT3.e3749_s2616_s2183_r6869_r6282_p2411/output.root')
#tc.Add('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/VVXAOD-00-00-22/mc15_13TeV.301257.Pythia8EvtGen_A14NNPDF23LO_Wprime_WZqqqq_m1000.merge.DAOD_EXOT3.e3749_s2616_s2183_r6869_r6282_p2411/output.root')
tc.Add('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/VVXAOD-00-00-22/mc15_13TeV.301262.Pythia8EvtGen_A14NNPDF23LO_Wprime_WZqqqq_m1500.merge.DAOD_EXOT3.e3743_s2608_s2183_r6869_r6282_p2411/output.root')
tc.Add('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/VVXAOD-00-00-22/mc15_13TeV.301267.Pythia8EvtGen_A14NNPDF23LO_Wprime_WZqqqq_m2000.merge.DAOD_EXOT3.e3749_s2616_s2183_r6869_r6282_p2411/output.root')
tc.Add('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/VVXAOD-00-00-22/mc15_13TeV.301272.Pythia8EvtGen_A14NNPDF23LO_Wprime_WZqqqq_m2500.merge.DAOD_EXOT3.e3743_s2608_s2183_r6869_r6282_p2411/output.root')
tc.Add('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/VVXAOD-00-00-22/mc15_13TeV.301277.Pythia8EvtGen_A14NNPDF23LO_Wprime_WZqqqq_m3000.merge.DAOD_EXOT3.e3743_s2608_s2183_r6869_r6282_p2411/output.root')
#fil = rt.TFile('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/VVXAOD-00-00-22/mc15_13TeV.301257.Pythia8EvtGen_A14NNPDF23LO_Wprime_WZqqqq_m1000.merge.DAOD_EXOT3.e3749_s2616_s2183_r6869_r6282_p2411/output.root','read')
#tree = fil.Get('outputTree')

topomass = np.zeros(1,dtype=float)#rt.std.vector(float)()
groomedmass = np.zeros(1,dtype=float)#rt.std.vector(float)()
truthmass = np.zeros(1,dtype=float)#rt.std.vector(float)()

outfile = rt.TFile('output.root','recreate')
otree = rt.TTree('outputTree','outputTree')
#otree.Branch('topomass',topomass,'topomass/D')
otree.Branch('groomedmass',groomedmass,'groomedmass/D')
otree.Branch('truthmass',truthmass,'truthmass/D')

entries = tc.GetEntries()



for i in range(entries):
    tc.GetEntry(i)
    if i%10000 == 0:
        print i

    # find leading truth
    maxjet = -1
    maxjet_pt = 0.0
    for j,x in enumerate(tc.jet_CamKt12Truth_m):
        if tc.jet_CamKt12Truth_pt[j] > maxjet_pt:
            maxjet = j
            maxjet_pt = tc.jet_CamKt12Truth_pt[j]
    if maxjet == -1 or tc.jet_CamKt12Truth_pt[maxjet] < 50*1000 or abs(tc.jet_CamKt12Truth_eta[maxjet]) > 1.2:
        continue

    # find leading groomed
    maxgr = -1
    maxgr_pt = 0.0
    for j,x in enumerate(tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_m):
        if tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_pt[j] > maxgr_pt and abs(tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_eta[j]) < 1.2:
            maxgr = j
            maxgr_pt = tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_pt[j]
            
    if maxgr == -1:
        continue

    # match to truth boson
    
    count = 0
    truthBosonIndex = -1
    for jet_i,x in enumerate(tc.truthBoson_ID):
        
        #// delta R matching of groomed jet with truth boson
        if DeltaR(tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_eta[maxgr],tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_phi[maxgr],tc.truthBosons[jet_i].Eta(),tc.truthBosons[jet_i].Phi())<(0.75*1.2) and tc.truthBosons[jet_i].Pt() > 5*1000.0 and abs(tc.truthBosons[jet_i].Eta()) < 6:


            #// if there is a Z boson within this radius, regardless of whether or not there is also a W, we veto the event
            if abs(tc.truthBoson_ID[jet_i]) == 24:
                truthBosonIndex = count;

        count+=1
    if truthBosonIndex == -1:
        continue
    

    # is there a matching truth to groomed?
    chosenLeadTruthJetIndex = -1
    
    for jet_i,x in enumerate(tc.jet_CamKt12Truth_pt):
        dr = DeltaR(tc.jet_CamKt12Truth_eta[jet_i],tc.jet_CamKt12Truth_phi[jet_i],tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_eta[maxgr],tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_phi[maxgr])
        if dr<(0.75*1.2) and abs(tc.jet_CamKt12Truth_eta[jet_i])<4.5:
            chosenLeadTruthJetIndex=jet_i;
 
    if chosenLeadTruthJetIndex == -1:
        continue

    # fill the branches
    #topomass[0] = tc.jet_CamKt12Truth_m[]
    groomedmass[0] = tc.jet_CamKt12LCTopoBDRSFilteredMU100Y4_m[maxgr]
    truthmass[0] = tc.jet_CamKt12Truth_m[chosenLeadTruthJetIndex]
    otree.Fill()
    
#tc.Close()

outfile.Write()
outfile.Close()
    
