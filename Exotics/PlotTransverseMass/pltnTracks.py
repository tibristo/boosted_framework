import ROOT as rt
import sys

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptFit(1111)


file_in_name = sys.argv[1]
signal = False

if sys.argv[2].lower() == 'signal':
    signal = True
    
f = rt.TFile.Open(file_in_name)

tree = f.Get('outputTree')

# draw the nTracks vs parentID plots
tcanv = rt.TCanvas('canv1')
if signal:

    histw = 'hist_w'
    tree.Draw("nTracks>>"+histw,"(truth_id==24)")
    histw = rt.gDirectory.Get(histw)
    histw.SetTitle("nTracks for W boson parents")
    histw.GetXaxis().SetTitle("nTracks")
    histw.Draw()
    tcanv.SaveAs(file_in_name.replace('.root','_w_ca12.png'))
    
    histz = 'hist_z'
    tree.Draw("nTracks>>"+histz,"(truth_id==23)")
    histz = rt.gDirectory.Get(histz)
    histz.SetTitle("nTracks for Z boson parents")
    histz.GetXaxis().SetTitle("nTracks")
    histz.Draw()
    tcanv.SaveAs(file_in_name.replace('.root','_z_ca12.png'))

else:
    histjz = 'hist_qcd'
    tree.Draw("nTracks>>"+histjz,"(truth_id==99)")
    histjz = rt.gDirectory.Get(histjz)
    histjz.SetTitle("nTracks for QCD boson parents")
    histjz.GetXaxis().SetTitle("nTracks")
    histjz.Draw()
    tcanv.SaveAs(file_in_name.replace('.root','_qcd_ca12.png'))

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

gr_topo = rt.TH2F("groomed_vs_topo","groomed_vs_topo",100,150,2000,100,150,2000)
gr_topo.GetYaxis().SetTitle("Groomed pT")
gr_topo.GetXaxis().SetTitle("Topo pT")

mc_channel_number = 0

for n in range(entries):
    tree.GetEntry(n)
    mc_channel_number = tree.mc_channel_number[0]
    nTrk_pt.Fill(tree.nTracks[0], tree.topo_pt[0])
    nTrk_m.Fill(tree.nTracks[0], tree.topo_m[0])
    pt_nTrk.Fill(tree.topo_pt[0], tree.nTracks[0])
    m_nTrk.Fill(tree.topo_m[0], tree.nTracks[0])
    gr_topo.Fill(tree.topo_pt[0], tree.groomed_pt[0])




#create a linear function with set parameters - always want it to start at origin


f2 = rt.TF1("f2","[0]+[1]*x",0,150)#,0,1000)
f2.FixParameter(0,0)

#mc_channel_number

# draw and save
#nTrk_pt.Fit("pol1")
#nTrk_pt.Draw()
#tcanv.SaveAs(file_in_name.replace('.root','_ntrk_pt_ca12.png'))

nTrk_m.Fit("f2")
with open('nTrk_m_fitparams_ca12.csv','a') as fits:
    fits.write(str(mc_channel_number) + ',' + str(f2.GetParameter(1))+'\n')
nTrk_m.Draw()
tcanv.SaveAs(file_in_name.replace('.root','_ntrk_m_ca12.png'))

#pt_nTrk.Fit("pol1")
#pt_nTrk.Draw()
#tcanv.SaveAs(file_in_name.replace('.root','_pt_ntrk_ca12.png'))
m_nTrk.Fit("f2")
with open('m_nTrk_fitparams_ca12.csv','a') as fits:
    fits.write(str(mc_channel_number) + ',' + str(f2.GetParameter(1))+'\n')
m_nTrk.Draw()
tcanv.SaveAs(file_in_name.replace('.root','_m_ntrk_ca12.png'))

gr_topo.Fit("pol1")
with open('gr_topo_pt_fitparams_ca12.csv','a') as fits:
    fits.write(str(mc_channel_number) + ',' + str(gr_topo.GetFunction("pol1").GetParameter(1))+'\n')
gr_topo.Draw()
tcanv.SaveAs(file_in_name.replace('.root','_pt_comp_ca12.png'))
#fits.close()
f.Close()
