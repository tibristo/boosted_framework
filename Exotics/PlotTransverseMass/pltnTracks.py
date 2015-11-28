import ROOT as rt
import sys

rt.gROOT.SetBatch(True)

file_in_name = sys.argv[1]
signal = False

if sys.argv[2].lower() == 'signal':
    signal = True
    
f = rt.TFile.Open(file_in_name)

tree = f.Get('outputTree')

# draw the nTracks vs parentID plots
if signal:
    tcanv = rt.TCanvas('canv1')
    histw = 'hist_w'
    tree.Draw("nTracks>>"+histw,"(truth_id==24)")
    histw = rt.gDirectory.Get(histw)
    histw.SetTitle("nTracks for W boson parents")
    histw.GetXaxis().SetTitle("nTracks")
    histw.Draw()
    tcanv.SaveAs(file_in_name.replace('.root','_w.png'))
    
    histz = 'hist_z'
    tree.Draw("nTracks>>"+histz,"(truth_id==23)")
    histz = rt.gDirectory.Get(histz)
    histz.SetTitle("nTracks for Z boson parents")
    histz.GetXaxis().SetTitle("nTracks")
    histz.Draw()
    tcanv.SaveAs(file_in_name.replace('.root','_z.png'))

else:
    histjz = 'hist_qcd'
    tree.Draw("nTracks>>"+histjz,"(truth_id==99)")
    histjz = rt.gDirectory.Get(histjz)
    histjz.SetTitle("nTracks for QCD boson parents")
    histjz.GetXaxis().SetTitle("nTracks")
    histjz.Draw()
    tcanv.SaveAs(file_in_name.replace('.root','_qcd.png'))

# draw 2D plot with nTracks on x-axis, pT on y-axis
entries = tree.GetEntries()
nTrk_pt = rt.TH2F("nTrk_pT","nTrk_pT",100,0,150,100,150,2000)
nTrk_pt.GetXaxis().SetTitle("nTracks")
nTrk_pt.GetYaxis().SetTitle("Topo Jet pT")
pt_nTrk = rt.TH2F("pT_nTrk","pT_nTrk",100,150,2000,100,0,150)
pt_nTrk.GetYaxis().SetTitle("nTracks")
pt_nTrk.GetXaxis().SetTitle("Topo Jet pT")
nTrk_m = rt.TH2F("nTrk_m","nTrk_m",100,0,150,100,0,1000)
nTrk_m.GetXaxis().SetTitle("nTracks")
nTrk_m.GetYaxis().SetTitle("Topo Jet Mass")
m_nTrk = rt.TH2F("m_nTrk","m_nTrk",100,0,1000,100,0,150)
m_nTrk.GetYaxis().SetTitle("nTracks")
m_nTrk.GetXaxis().SetTitle("Topo Jet Mass")

for n in range(entries):
    tree.GetEntry(n)
    nTrk_pt.Fill(tree.nTracks[0], tree.topo_pt[0])
    nTrk_m.Fill(tree.nTracks[0], tree.topo_m[0])
    pt_nTrk.Fill(tree.topo_pt[0], tree.nTracks[0])
    m_nTrk.Fill(tree.topo_m[0], tree.nTracks[0])

# draw and save
nTrk_pt.Draw()
tcanv.SaveAs(file_in_name.replace('.root','_ntrk_pt.png'))
nTrk_m.Draw()
tcanv.SaveAs(file_in_name.replace('.root','_ntrk_m.png'))

pt_nTrk.Draw()
tcanv.SaveAs(file_in_name.replace('.root','_pt_ntrk.png'))
m_nTrk.Draw()
tcanv.SaveAs(file_in_name.replace('.root','_m_ntrk.png'))

f.Close()
